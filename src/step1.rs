//! Step 1: per-proteome k-mer simplification. Proteins sharing a k-mer are grouped and
//! only the longest one of each group is kept for the next steps.

use crate::{create, fresh_dir, Opts, R};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use rayon::prelude::*;
use rustc_hash::{FxBuildHasher, FxHashMap};
use std::hash::BuildHasher;
use std::io::Write;
use std::path::Path;

pub fn run(o: &Opts) -> R<()> {
    let dir = Path::new(o.str("dir"));
    let (k, min_aa): (usize, usize) = (o.get("kmer_size")?, o.get("kmer_min_aa")?);
    let (dir_arg, out) = (o.str("dir_arg"), o.str("output")); // as given by the user
    println!(" --- STEP 1: k-mer simplification\n input dir: {dir_arg}\n output dir: {}\n kmer size: {k}\n kmer nb aa: {min_aa}",
        if out.is_empty() { "." } else { out });
    if !dir.is_dir() { return Err(format!("the directory '{dir_arg}' does not exist").into()) }
    fresh_dir("dir_step1")?;

    // input files, by name (species index and protein ids do not depend on compression)
    let mut files: Vec<String> = std::fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().is_file())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| is_proteome(n))
        .collect();
    if files.is_empty() { return Err(format!("no proteome file (.fas/.fasta/.faa, optionally .gz) in {dir_arg}").into()) }
    files.sort();

    // pass 1: protein names only (ids, duplicate check); sequences are re-read per file below
    let names: Vec<Vec<String>> = files.par_iter()
        .map(|f| Ok(read_fasta(&dir.join(f))?.into_iter().map(|(n, _)| n).collect()))
        .collect::<R<_>>()?;
    let mut offset = vec![0u32];
    for n in &names { offset.push(offset.last().unwrap() + n.len() as u32) }
    println!(" {} input files\n {} sequences", files.len(), offset.last().unwrap());

    let mut seen: FxHashMap<&str, usize> = FxHashMap::default();
    let mut dupli = vec![];
    for (i, ns) in names.iter().enumerate() {
        for n in ns {
            if let Some(&j) = seen.get(n.as_str()) { dupli.push(format!("{n}\t{}\t{}", files[j], files[i])) }
            else { seen.insert(n, i); }
        }
    }
    if !dupli.is_empty() {
        println!(" WARNING: some protein names are present multiple times, see {}", crate::shown("dir_step1/duplicate_names.txt").display());
        std::fs::write("dir_step1/duplicate_names.txt", format!("#protein_name\tfile1\tfile2\n{}\n", dupli.join("\n")))?;
    }

    // pass 2: simplify each proteome and write its reduced version
    let reps: Vec<Vec<u32>> = files.par_iter().enumerate().map(|(i, f)| {
        let seqs: Vec<Vec<u8>> = read_fasta(&dir.join(f))?.into_iter().map(|(_, s)| s).collect();
        let rep = simplify(&seqs, k, min_aa);
        let mut out = create(format!("dir_step1/{i}.fas"))?;
        for (j, s) in seqs.iter().enumerate() {
            if rep[j] as usize == j { writeln!(out, ">{}\n{}", offset[i] + j as u32, String::from_utf8_lossy(s))? }
        }
        Ok(rep.into_iter().map(|r| r + offset[i]).collect())
    }).collect::<R<_>>()?;

    let mut species = create("dir_step1/species.tsv")?;
    let mut prots = create("dir_step1/proteins.tsv")?;
    let mut log = create("dir_step1/log_step1.txt")?;
    writeln!(log, "#index\tfile_name\tnb_initial\tnb_final")?;
    let mut nb_final = 0;
    for (i, f) in files.iter().enumerate() {
        writeln!(species, "{i}\t{f}")?;
        let kept = reps[i].iter().enumerate().filter(|&(j, &r)| r == offset[i] + j as u32).count();
        writeln!(log, "{i}\t{f}\t{}\t{kept}", names[i].len())?;
        nb_final += kept;
        for (j, n) in names[i].iter().enumerate() {
            writeln!(prots, "{}\t{i}\t{}\t{n}", offset[i] + j as u32, reps[i][j])?;
        }
    }
    println!(" -> {nb_final} proteins saved for the next step\n");
    Ok(())
}

/// Proteome files: .fas, .fasta or .faa (any case), optionally gzipped.
fn is_proteome(name: &str) -> bool {
    let n = name.to_lowercase();
    let n = n.strip_suffix(".gz").unwrap_or(&n);
    [".fas", ".fasta", ".faa"].iter().any(|e| n.ends_with(e))
}

/// Returns (name = first word of header, sequence) records. Reads gzipped files too.
pub fn read_fasta(path: &Path) -> R<Vec<(String, Vec<u8>)>> {
    let f = File::open(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.to_string_lossy().to_lowercase().ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    };
    let mut v: Vec<(String, Vec<u8>)> = vec![];
    for l in reader.lines() {
        let l = l?;
        let l = l.trim_end();
        if let Some(h) = l.strip_prefix('>') {
            v.push((h.split_whitespace().next().unwrap_or("").to_string(), vec![]));
        } else if let Some(last) = v.last_mut() {
            last.1.extend_from_slice(l.as_bytes());
        }
    }
    Ok(v)
}

/// Groups sequences sharing at least one informative k-mer (>= min_aa distinct
/// residues, no 'X' or '*'), transitively. Returns, for each sequence, the index of
/// the longest sequence of its group (ties: lowest index).
fn simplify(seqs: &[Vec<u8>], k: usize, min_aa: usize) -> Vec<u32> {
    // every informative kmer as one u64: top bits of its hash | sequence | start (8 bytes
    // per kmer); sorting by hash then kmer brings equal kmers together
    let bits = |x: usize| usize::BITS - x.leading_zeros(); // bits needed for 0..=x
    let sb = bits(seqs.len().saturating_sub(1));
    let pb = bits(seqs.iter().map(|s| s.len().saturating_sub(k)).max().unwrap_or(0));
    // the hash keeps the remaining top bits (fewer bits only mean more exact comparisons)
    let top = u64::MAX.checked_shl(sb + pb).unwrap_or(0);
    let hash = |x: u64| x & top;
    let seq = |x: u64| ((x >> pb) & ((1u64 << sb) - 1)) as usize;
    let start = |x: u64| (x & ((1u64 << pb) - 1)) as usize;
    let mut kmers: Vec<u64> = Vec::with_capacity(seqs.iter().map(|s| (s.len() + 1).saturating_sub(k)).sum());
    for (i, s) in seqs.iter().enumerate() {
        // sliding window of residue counts
        let (mut cnt, mut distinct, mut bad) = ([0u32; 256], 0, 0);
        for j in 0..s.len() {
            let c = s[j] as usize;
            if cnt[c] == 0 { distinct += 1 }
            cnt[c] += 1;
            if c == b'X' as usize || c == b'*' as usize { bad += 1 }
            if j >= k {
                let c = s[j - k] as usize;
                cnt[c] -= 1;
                if cnt[c] == 0 { distinct -= 1 }
                if c == b'X' as usize || c == b'*' as usize { bad -= 1 }
            }
            if j + 1 >= k && distinct >= min_aa && bad == 0 {
                let start = j + 1 - k;
                kmers.push(FxBuildHasher.hash_one(&s[start..=j]) & top | (i as u64) << pb | start as u64);
            }
        }
    }
    let kmer = |x: u64| &seqs[seq(x)][start(x)..start(x) + k];
    kmers.sort_unstable_by(|&a, &b| hash(a).cmp(&hash(b)).then_with(|| kmer(a).cmp(kmer(b)))); // exact despite hash collisions
    let mut uf: Vec<u32> = (0..seqs.len() as u32).collect();
    for w in kmers.windows(2) {
        if hash(w[0]) == hash(w[1]) && kmer(w[0]) == kmer(w[1]) {
            let (a, b) = (find(&mut uf, seq(w[0]) as u32), find(&mut uf, seq(w[1]) as u32));
            uf[a.max(b) as usize] = a.min(b); // root = smallest index
        }
    }
    let mut best: Vec<u32> = (0..seqs.len() as u32).collect();
    for i in 0..seqs.len() {
        let r = find(&mut uf, i as u32) as usize;
        if seqs[i].len() > seqs[best[r] as usize].len() { best[r] = i as u32 }
    }
    (0..seqs.len()).map(|i| best[find(&mut uf, i as u32) as usize]).collect()
}

fn find(uf: &mut [u32], mut x: u32) -> u32 {
    while uf[x as usize] != x {
        uf[x as usize] = uf[uf[x as usize] as usize];
        x = uf[x as usize];
    }
    x
}

#[cfg(test)]
mod tests {
    #[test]
    fn proteome_extensions() {
        for n in ["a.fas", "a.FASTA", "a.faa.gz", "a.Faa.GZ"] { assert!(super::is_proteome(n), "{n}") }
        for n in ["a.fa", "a.txt", "a.gz", "fasta"] { assert!(!super::is_proteome(n), "{n}") }
    }

    #[test]
    fn simplify_keeps_longest() {
        let a = b"ACDEFGHIKLMNPQ".to_vec();
        let mut b = a.clone();
        b.extend(b"RSTVW");
        let c = b"WWWWWWWWWWWWWW".to_vec(); // low complexity: never grouped
        let d = b"ACDEXGHIKLMNPQ".to_vec(); // contains X
        let rep = super::simplify(&[a, b, c.clone(), c, d], 10, 5);
        assert_eq!(rep, vec![1, 1, 2, 3, 4]);
    }
}
