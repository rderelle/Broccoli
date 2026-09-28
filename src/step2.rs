//! Step 2: phylomes. DIAMOND searches of every reduced proteome against all of them
//! (N best targets per species); then, for every query with duplicated species among
//! its hits, a trimmed alignment built from the pairwise HSPs is turned into a
//! FastTree BioNJ tree, midpoint-rooted. Queries without duplications are kept as
//! similarity groups.
//!
//! Default: one DIAMOND run per pair of species, as in Broccoli v1.
//! -combined_search: one DIAMOND run (all threads) of all proteins against one combined
//! database, per-species limit with --taxon-k (fake taxonomy: one taxid per species)
//! and --dbsize set to the mean proteome size so e-values match per-species searches.
//! DIAMOND output is streamed; batches of queries go to single-threaded workers
//! (alignments + one FastTree run per batch), so memory does not grow with the data.
//!
//! Outputs in dir_step2/:
//!   phylomes.tsv  `T<tab>query<tab>rooted newick`  or  `S<tab>query<tab>id id ...`
//!   hits.tsv      `query<tab>target<tab>min qstart<tab>max qend` (self hits excluded)

use crate::{create, fresh_dir, run_cmd, step1::read_fasta, Opts, Proteins, R};
use crate::tree::Tree;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{mpsc, Arc, Mutex};

const MAX_POSITIONS: usize = 4988; // FastTree limitation (5000 characters per line)
const BATCH: usize = 100; // queries per worker job (= alignments per FastTree run, at most)
const OUTFMT: [&str; 7] = ["6", "qseqid", "sseqid", "qstart", "qend", "sstart", "cigar"];

struct Hit { t: u32, qs: u32, qe: u32, ss: u32, cigar: Box<[u8]> } // 32 bytes (+ cigar)

type Batch = Vec<(u32, Vec<Hit>)>;

/// Trimmed alignment: (protein id, row) pairs.
type Ali = Vec<(u32, Vec<u8>)>;

enum Query { Tree(Ali), Sim(Vec<u32>), Empty }

/// Everything the workers share.
struct Ctx<'a> {
    sp: &'a [u32],
    seqs: &'a [Vec<u8>],
    nsp: usize,
    nb_hits: usize,
    max_gap: f64,
    fasttree: &'a str,
    method: &'a [&'a str],
    counts: Vec<[AtomicUsize; 4]>, // per species: phylo, no phylo, empty alignment, tree problem
    out: Mutex<(BufWriter<File>, BufWriter<File>)>,
}

pub fn run(o: &Opts) -> R<()> {
    let (evalue, nb_hits, max_gap): (f64, usize, f64) = (o.get("e_value")?, o.get("nb_hits")?, o.get("max_gap")?);
    let threads: usize = o.get("threads")?;
    let exact = o.str("combined_search") != "true";
    println!(" --- STEP 2: phylomes\n e_value: {evalue}\n nb_hits: {nb_hits}\n gaps: {max_gap}\n phylogenies: {}\n search: {}",
        o.str("phylogenies"), if exact { "one DIAMOND run per pair of species" } else { "combined database" });
    // FastTree options, or ["nj"] for the built-in neighbor-joining
    let method: &[&str] = match o.str("phylogenies") { "bionj" => &["-noml", "-nome"], "me" => &["-noml"], "nj" => &["nj"], _ => &[] };

    let mut prots = Proteins::load()?;
    prots.name = vec![]; // protein names are not needed in this step
    let nsp = prots.files.len();
    fresh_dir("dir_step2")?;
    std::fs::create_dir("dir_step2/db")?;

    // reduced sequences, indexed by protein id
    let mut seqs: Vec<Vec<u8>> = vec![vec![]; prots.sp.len()];
    for i in 0..nsp {
        for (n, s) in read_fasta(Path::new(&format!("dir_step1/{i}.fas")))? { seqs[n.parse::<usize>()?] = s }
    }

    let ctx = Ctx {
        sp: &prots.sp, seqs: &seqs, nsp, nb_hits, max_gap, fasttree: o.str("path_fasttree"), method,
        counts: (0..nsp).map(|_| Default::default()).collect(),
        out: Mutex::new((create("dir_step2/phylomes.tsv")?, create("dir_step2/hits.tsv")?)),
    };
    let search = Search { o, nsp, nb_hits, evalue, threads, letters: seqs.iter().map(Vec::len).sum() };

    // bounded queue: the DIAMOND reader blocks when workers lag behind
    let (tx, rx) = mpsc::sync_channel::<Batch>(2 * threads);
    let rx = Arc::new(Mutex::new(rx));
    let failed = AtomicBool::new(false);
    std::thread::scope(|s| -> R<()> {
        let workers: Vec<_> = (0..threads).map(|_| {
            let (rx, ctx, failed) = (rx.clone(), &ctx, &failed);
            s.spawn(move || -> R<()> {
                loop {
                    // the lock is released at the end of this statement, not held while processing
                    let Ok(b) = rx.lock().unwrap().recv() else { break };
                    if failed.load(Relaxed) { break }
                    if let Err(e) = process(b, ctx) { failed.store(true, Relaxed); return Err(e) }
                }
                Ok(())
            })
        }).collect();
        drop(rx); // workers own the receiver: if they all stop, sending fails
        let fed = if exact { search.exact(&tx, &prots.files) } else { search.combined(&tx, &prots) };
        drop(tx);
        for w in workers { w.join().expect("worker panicked")? }
        fed
    })?;

    let (mut a, mut b) = ctx.out.into_inner().unwrap();
    a.flush()?;
    b.flush()?;
    let mut log = create("dir_step2/log_step2.txt")?;
    writeln!(log, "#species_file\tnb_phylo\tnb_NO_phylo\tnb_empty_ali\tnb_pbm_tree")?;
    for (f, c) in prots.files.iter().zip(&ctx.counts) {
        writeln!(log, "{f}\t{}\t{}\t{}\t{}", c[0].load(Relaxed), c[1].load(Relaxed), c[2].load(Relaxed), c[3].load(Relaxed))?;
    }
    std::fs::remove_dir_all("dir_step2/db")?;
    println!(" done\n");
    Ok(())
}

struct Search<'a> { o: &'a Opts, nsp: usize, nb_hits: usize, evalue: f64, threads: usize, letters: usize }

impl Search<'_> {
    fn blastp(&self, db: &str, query: &str, threads: usize) -> Command {
        let mut c = Command::new(self.o.str("path_diamond"));
        c.args(["blastp", "--quiet", "--more-sensitive", "--threads", &threads.to_string(), "--db", db, "--query", query,
            "-e", &self.evalue.to_string(), "--outfmt"]).args(OUTFMT);
        c
    }

    /// All proteins against one database holding every proteome.
    fn combined(&self, tx: &mpsc::SyncSender<Batch>, prots: &Proteins) -> R<()> {
        let mut all = create("dir_step2/db/all.fas")?;
        for i in 0..self.nsp { std::io::copy(&mut File::open(format!("dir_step1/{i}.fas"))?, &mut all)?; }
        all.flush()?;
        // fake taxonomy: species i = taxid i + 2, child of the root (taxid 1)
        let mut map = create("dir_step2/db/taxonmap.tsv")?;
        writeln!(map, "accession\taccession.version\ttaxid\tgi")?;
        for (id, (&s, &r)) in prots.sp.iter().zip(&prots.rep).enumerate() {
            if r as usize == id { writeln!(map, "{id}\t{id}\t{}\t0", s + 2)? }
        }
        map.flush()?;
        let (mut nodes, mut names) = ("1\t|\t1\t|\tno rank\t|\t\t|\n".to_string(), "1\t|\troot\t|\t\t|\tscientific name\t|\n".to_string());
        for s in 0..self.nsp {
            nodes += &format!("{}\t|\t1\t|\tspecies\t|\t\t|\n", s + 2);
            names += &format!("{}\t|\tspecies{s}\t|\t\t|\tscientific name\t|\n", s + 2);
        }
        std::fs::write("dir_step2/db/nodes.dmp", nodes)?;
        std::fs::write("dir_step2/db/names.dmp", names)?;
        println!(" build DIAMOND database");
        run_cmd(Command::new(self.o.str("path_diamond")).args(["makedb", "--quiet", "--in", "dir_step2/db/all.fas", "--db", "dir_step2/db/all",
            "--taxonmap", "dir_step2/db/taxonmap.tsv", "--taxonnodes", "dir_step2/db/nodes.dmp", "--taxonnames", "dir_step2/db/names.dmp"]), None)?;

        println!(" similarity search and phylomes ... be patient");
        let mut cmd = self.blastp("dir_step2/db/all", "dir_step2/db/all.fas", self.threads);
        cmd.args(["-k", "0", "--taxon-k", &self.nb_hits.to_string(), "--dbsize", &(self.letters / self.nsp).to_string()])
            .stdin(Stdio::null()).stdout(Stdio::piped()).stderr(File::create("dir_step2/diamond.log")?);
        let mut child = cmd.spawn().map_err(|e| format!("cannot run {}: {e}", self.o.str("path_diamond")))?;
        let res = stream(BufReader::new(child.stdout.take().unwrap()), tx);
        if res.is_err() { let _ = child.kill(); }
        let status = child.wait()?;
        res?;
        if !status.success() {
            return Err(format!("DIAMOND failed: {}", std::fs::read_to_string("dir_step2/diamond.log").unwrap_or_default()).into());
        }
        Ok(())
    }

    /// One search per pair of species (Broccoli v1).
    fn exact(&self, tx: &mpsc::SyncSender<Batch>, files: &[String]) -> R<()> {
        (0..self.nsp).into_par_iter().try_for_each(|i| {
            run_cmd(Command::new(self.o.str("path_diamond")).args(["makedb", "--quiet", "--in", &format!("dir_step1/{i}.fas"), "--db", &format!("dir_step2/db/{i}")]), None).map(drop)
        })?;
        let t = (self.threads / self.nsp).max(1);
        for (q, file) in files.iter().enumerate() {
            println!(" phylome {}/{}: {file}", q + 1, self.nsp);
            // each DIAMOND output is parsed as soon as it is produced (the raw text is not kept)
            let parsed: Vec<Vec<(u32, Hit)>> = (0..self.nsp).into_par_iter().map(|db| {
                let mut c = self.blastp(&format!("dir_step2/db/{db}"), &format!("dir_step1/{q}.fas"), t);
                let out = run_cmd(c.args(["-k", &self.nb_hits.to_string()]), None)?;
                out.split(|&c| c == b'\n').filter(|l| !l.is_empty()).map(parse_hit).collect()
            }).collect::<R<_>>()?;
            let mut by_query: FxHashMap<u32, Vec<Hit>> = FxHashMap::default();
            for (q, h) in parsed.into_iter().flatten() { by_query.entry(q).or_default().push(h) }
            let mut queries: Vec<(u32, Vec<Hit>)> = by_query.into_iter().collect();
            queries.sort_by_key(|x| x.0);
            while !queries.is_empty() {
                let rest = queries.split_off(queries.len().min(BATCH));
                tx.send(std::mem::replace(&mut queries, rest)).map_err(|_| "phylome workers stopped")?;
            }
        }
        Ok(())
    }
}

fn parse_hit(line: &[u8]) -> R<(u32, Hit)> {
    let f: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
    if f.len() < 6 { return Err(format!("unexpected DIAMOND output: {}", String::from_utf8_lossy(line)).into()) }
    let num = |x: &[u8]| -> R<usize> { Ok(std::str::from_utf8(x)?.trim().parse()?) };
    Ok((num(f[0])? as u32, Hit { t: num(f[1])? as u32, qs: num(f[2])? as u32, qe: num(f[3])? as u32, ss: num(f[4])? as u32, cigar: f[5].trim_ascii().into() }))
}

/// Group DIAMOND output lines by query (they come grouped, in query order) and send batches.
fn stream(mut r: impl BufRead, tx: &mpsc::SyncSender<Batch>) -> R<()> {
    let (mut batch, mut line): (Batch, Vec<u8>) = (vec![], vec![]);
    while r.read_until(b'\n', &mut line)? > 0 {
        let (q, h) = parse_hit(&line)?;
        line.clear();
        match batch.last_mut() {
            Some((last, hits)) if *last == q => { hits.push(h); continue }
            Some((last, _)) if *last > q => return Err("DIAMOND output is not grouped by query".into()),
            _ => {}
        }
        if batch.len() == BATCH { tx.send(std::mem::take(&mut batch)).map_err(|_| "phylome workers stopped")? }
        batch.push((q, vec![h]));
    }
    if !batch.is_empty() { tx.send(batch).map_err(|_| "phylome workers stopped")? }
    Ok(())
}

/// Worker job: alignments, trees and output lines for a batch of queries (single thread).
fn process(batch: Batch, c: &Ctx) -> R<()> {
    let (mut phy, mut hit_lines, mut alis) = (String::new(), String::new(), vec![]);
    for (q, mut hits) in batch {
        // targets in species order (as in per-species searches), at most nb_hits per species
        // (DIAMOND --taxon-k may report one more)
        hits.sort_by_key(|h| c.sp[h.t as usize]);
        let (mut cur, mut kept): (u32, Vec<u32>) = (u32::MAX, vec![]);
        hits.retain(|h| {
            let s = c.sp[h.t as usize];
            if s != cur { cur = s; kept.clear() }
            if kept.contains(&h.t) { return true }
            if kept.len() == c.nb_hits { return false }
            kept.push(h.t);
            true
        });

        // hit coordinates, one line per (query, target) pair
        let mut agg: Vec<(u32, u32, u32)> = vec![];
        for h in hits.iter().filter(|h| h.t != q) {
            match agg.iter_mut().find(|a| a.0 == h.t) {
                Some(a) => { a.1 = a.1.min(h.qs); a.2 = a.2.max(h.qe) }
                None => agg.push((h.t, h.qs, h.qe)),
            }
        }
        for (t, s, e) in agg { hit_lines += &format!("{q}\t{t}\t{s}\t{e}\n") }

        let count = &c.counts[c.sp[q as usize] as usize];
        match analyse(q, &hits, c.sp, c.seqs, c.nsp, c.max_gap) {
            Query::Sim(ids) => {
                count[1].fetch_add(1, Relaxed);
                phy += &format!("S\t{q}\t{}\n", ids.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "));
            }
            Query::Empty => { count[0].fetch_add(1, Relaxed); count[2].fetch_add(1, Relaxed); }
            Query::Tree(a) => { count[0].fetch_add(1, Relaxed); alis.push((q, a)) }
        }
    }

    // unrooted trees, in the order of `alis`
    let trees: Vec<String> = if c.method == ["nj"] {
        alis.iter().map(|(_, a)| {
            let (ids, rows): (Vec<u32>, Vec<&[u8]>) = a.iter().map(|(id, r)| (*id, r.as_slice())).unzip();
            crate::nj::tree(&ids, &rows)
        }).collect()
    } else if alis.is_empty() {
        vec![]
    } else {
        let mut input = String::new();
        for (_, a) in &alis {
            input += &format!("{}\t{}\n", a.len(), a[0].1.len());
            for (id, row) in a { input += &format!("{id:<10}{}\n", String::from_utf8_lossy(row)) }
        }
        let stdout = run_cmd(Command::new(c.fasttree).args(["-quiet", "-nosupport", "-fastest", "-bionj", "-pseudo"])
            .args(c.method).args(["-n", &alis.len().to_string()]), Some(input))?;
        // tree lines, in input order (skipping FastTree warnings)
        String::from_utf8(stdout)?.lines().map(str::trim)
            .filter(|l| !l.starts_with("Ign") && !l.starts_with("WARNING")).map(String::from).collect()
    };
    for (l, (q, _)) in trees.iter().zip(&alis) {
        if l.starts_with('(') {
            phy += &format!("T\t{q}\t{}\n", Tree::parse(l)?.midpoint_newick());
        } else if c.counts[c.sp[*q as usize] as usize][3].fetch_add(1, Relaxed) >= 100 {
            return Err("too many errors in phylogenetic analyses -> stopped".into());
        }
    }
    let mut out = c.out.lock().unwrap();
    out.0.write_all(phy.as_bytes())?;
    out.1.write_all(hit_lines.as_bytes())?;
    Ok(())
}

/// Similarity group if no species is duplicated among the hits, otherwise the
/// trimmed alignment (phylip) of the query and its HSPs.
fn analyse(q: u32, hits: &[Hit], sp: &[u32], seqs: &[Vec<u8>], nsp: usize, max_gap: f64) -> Query {
    let mut ids: Vec<u32> = hits.iter().map(|h| h.t).chain([q]).collect();
    ids.sort_unstable();
    ids.dedup();
    let mut per_sp = vec![0u32; nsp];
    for &i in &ids { per_sp[sp[i as usize] as usize] += 1 }
    if per_sp.iter().filter(|&&c| c > 0).count() < 2 || per_sp.iter().all(|&c| c < 2) {
        return Query::Sim(ids);
    }

    // one row per target, gaps outside its HSPs (query coordinates)
    let qlen = seqs[q as usize].len();
    let mut rows: Vec<(u32, Vec<u8>)> = vec![];
    for h in hits {
        let i = rows.iter().position(|r| r.0 == h.t).unwrap_or_else(|| { rows.push((h.t, vec![b'-'; qlen])); rows.len() - 1 });
        place_hsp(&mut rows[i].1, &seqs[h.t as usize], h);
    }
    if !rows.iter().any(|r| r.0 == q) { rows.push((q, seqs[q as usize].clone())) }

    let n = rows.len() as f64;
    let keep: Vec<usize> = (0..qlen)
        .filter(|&c| (rows.iter().filter(|r| r.1[c] == b'-').count() as f64 / n) < max_gap)
        .take(MAX_POSITIONS).collect();
    if keep.is_empty() { return Query::Empty }
    Query::Tree(rows.into_iter().map(|(id, row)| (id, keep.iter().map(|&c| row[c]).collect())).collect())
}

/// Write the target residues of an HSP into `row` (query coordinates), following the CIGAR.
fn place_hsp(row: &mut [u8], target: &[u8], h: &Hit) {
    let (mut pos, mut out, mut n) = (h.ss as usize - 1, h.qs as usize - 1, 0usize);
    for &c in h.cigar.iter() {
        if c.is_ascii_digit() { n = n * 10 + (c - b'0') as usize; continue }
        for _ in 0..n {
            match c {
                b'M' => { if out < row.len() { row[out] = *target.get(pos).unwrap_or(&b'-') } out += 1; pos += 1 }
                b'I' => { if out < row.len() { row[out] = b'-' } out += 1 }
                b'D' => pos += 1,
                _ => {}
            }
        }
        n = 0;
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn hsp_follows_cigar() {
        let mut row = vec![b'-'; 8];
        let h = super::Hit { t: 0, qs: 2, qe: 7, ss: 1, cigar: b"2M1I1D3M".as_slice().into() };
        super::place_hsp(&mut row, b"ABCDEFG", &h);
        assert_eq!(&row, b"-AB-DEF-");
    }
}
