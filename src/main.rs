//! Broccoli: orthology inference combining phylogenies and network analysis.
//! Rust rewrite of https://github.com/rderelle/Broccoli (Derelle et al. 2020, MBE).

mod nj;
mod step1;
mod step2;
mod step3;
mod step4;
mod tree;

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};

pub type R<T> = Result<T, Box<dyn std::error::Error + Send + Sync>>;

const HELP: &str = concat!("
            Broccoli v", env!("CARGO_PKG_VERSION"), "

 general options:
  -steps            steps to be performed, comma separated [default = 1,2,3,4]
  -threads          number of threads [default = 1]

 STEP 1  kmer clustering:
  -dir              directory containing the proteome files [required]
                    (.fas, .fasta or .faa, any case, optionally gzipped)
  -kmer_size        length of kmers [default = 100]
  -kmer_min_aa      minimum nb of different aa a kmer should have [default = 15]

 STEP 2  phylomes:
  -path_diamond     path of DIAMOND [default = diamond]
  -path_fasttree    path of FastTree [default = fasttree]
  -e_value          e-value for similarity search [default = 0.001]
  -nb_hits          max nb of hits per species [default = 6]
  -max_gap          max fraction of gap per position [default = 0.7]
  -phylogenies      FastTree 'bionj', 'me' or 'ml', or 'nj' (built-in NJ, much faster) [default = bionj]
  -combined_search  one DIAMOND search against a combined database (faster, near-identical hits)
                    [default: one DIAMOND search per pair of species, as Broccoli v1]

 STEP 3  network analysis:
  -sp_overlap       max ratio of overlapping species in phylogenetic trees [default = 0.5]
  -min_weight       min weight for an edge to be kept in the orthology network [default = 0.1]
  -min_nb_hits      spurious hits: min number of hits belonging to the OG [default = 2]
  -chimeric_shared  chimeric prot: min fraction of connected nodes in each OG [default = 0.5]
  -chimeric_nb_sp   chimeric prot: min nb of species in OGs involved in gene-fusions [default = 3]

 STEP 4  orthologous pairs:
  -ratio_ortho      limit ratio ortho/total [default = 0.5]
  -not_same_sp      ignore ortho relationships between proteins of the same species
");

const DEFAULTS: &[(&str, &str)] = &[
    ("steps", "1,2,3,4"), ("threads", "1"),
    ("dir", ""), ("kmer_size", "100"), ("kmer_min_aa", "15"),
    ("path_diamond", "diamond"), ("path_fasttree", "fasttree"), ("e_value", "0.001"),
    ("nb_hits", "6"), ("max_gap", "0.7"), ("phylogenies", "bionj"),
    ("sp_overlap", "0.5"), ("min_weight", "0.1"), ("min_nb_hits", "2"),
    ("chimeric_shared", "0.5"), ("chimeric_nb_sp", "3"),
    ("ratio_ortho", "0.5"), ("not_same_sp", "false"), ("combined_search", "false"),
];

/// Command-line options, keyed by name (without the leading dash).
pub struct Opts(HashMap<String, String>);

impl Opts {
    fn parse() -> R<Self> {
        let mut m: HashMap<String, String> =
            DEFAULTS.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect();
        let mut args = std::env::args().skip(1);
        while let Some(a) = args.next() {
            let k = a.trim_start_matches('-').to_string();
            match k.as_str() {
                "h" | "help" => { print!("{HELP}"); std::process::exit(0) }
                "not_same_sp" | "combined_search" => { m.insert(k, "true".into()); }
                _ if m.contains_key(&k) => {
                    let v = args.next().ok_or(format!("missing value for -{k}"))?;
                    m.insert(k, v);
                }
                _ => return Err(format!("unknown option '{a}' (see -help)").into()),
            }
        }
        Ok(Opts(m))
    }

    pub fn str(&self, k: &str) -> &str { &self.0[k] }

    pub fn get<T: std::str::FromStr>(&self, k: &str) -> R<T> {
        self.0[k].parse().map_err(|_| format!("invalid value '{}' for -{k}", self.0[k]).into())
    }
}

/// Protein table produced by step 1. Protein ids are 0..n, contiguous per species.
pub struct Proteins {
    pub files: Vec<String>, // species index -> proteome file name
    pub sp: Vec<u32>,       // protein id -> species index
    pub rep: Vec<u32>,      // protein id -> id of its kmer-cluster representative
    pub name: Vec<String>,  // protein id -> original name
}

impl Proteins {
    pub fn load() -> R<Self> {
        let files = lines("dir_step1/species.tsv")?
            .map(|l| Ok(l?.split('\t').nth(1).unwrap_or("").to_string()))
            .collect::<R<Vec<_>>>()?;
        let (mut sp, mut rep, mut name) = (vec![], vec![], vec![]);
        for l in lines("dir_step1/proteins.tsv")? {
            let l = l?;
            let f: Vec<&str> = l.splitn(4, '\t').collect();
            sp.push(f[1].parse()?);
            rep.push(f[2].parse()?);
            name.push(f[3].to_string());
        }
        Ok(Proteins { files, sp, rep, name })
    }

    /// Representative id -> all proteins it stands for (itself included).
    pub fn members(&self) -> Vec<Vec<u32>> {
        let mut m = vec![vec![]; self.rep.len()];
        for (i, &r) in self.rep.iter().enumerate() { m[r as usize].push(i as u32) }
        m
    }
}

pub fn lines(path: impl AsRef<Path>) -> R<std::io::Lines<BufReader<File>>> {
    let p = path.as_ref();
    let f = File::open(p).map_err(|e| format!("cannot open {}: {e}", p.display()))?;
    Ok(BufReader::new(f).lines())
}

pub fn create(path: impl AsRef<Path>) -> R<BufWriter<File>> {
    Ok(BufWriter::new(File::create(path)?))
}

/// Recreate an empty output directory.
pub fn fresh_dir(path: &str) -> R<()> {
    if Path::new(path).exists() { fs::remove_dir_all(path)? }
    Ok(fs::create_dir_all(path)?)
}

/// Run an external program, feeding `input` on stdin, and return its stdout.
pub fn run_cmd(cmd: &mut Command, input: Option<String>) -> R<Vec<u8>> {
    let prog = cmd.get_program().to_string_lossy().to_string();
    let mut child = cmd
        .stdin(if input.is_some() { Stdio::piped() } else { Stdio::null() })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("cannot run {prog}: {e}"))?;
    // write stdin from another thread so a full stdout pipe cannot deadlock us
    let writer = input.map(|s| {
        let mut stdin = child.stdin.take().unwrap();
        std::thread::spawn(move || stdin.write_all(s.as_bytes()))
    });
    let out = child.wait_with_output()?;
    if !out.status.success() {
        return Err(format!("{prog} failed: {}", String::from_utf8_lossy(&out.stderr)).into());
    }
    if let Some(w) = writer { w.join().expect("stdin writer panicked")? }
    Ok(out.stdout)
}

fn main() {
    if let Err(e) = run() {
        eprintln!("\n            ERROR: {e}\n");
        std::process::exit(1);
    }
}

fn run() -> R<()> {
    let o = Opts::parse()?;
    println!("\n            Broccoli v{}\n", env!("CARGO_PKG_VERSION"));
    let start = std::time::Instant::now();

    let mut steps: Vec<u32> = o.str("steps").split(',').map(|s| s.trim().parse())
        .collect::<Result<_, _>>().map_err(|_| "the steps should be integers (-steps)")?;
    steps.sort();
    steps.dedup();
    if steps.windows(2).any(|w| w[1] != w[0] + 1) || steps.iter().any(|s| !(1..=4).contains(s)) {
        return Err("the steps should be consecutive, between 1 and 4 (-steps)".into());
    }
    if steps.contains(&1) && o.str("dir").is_empty() {
        return Err("you need to specify an input directory with -dir (see -help)".into());
    }
    if !["bionj", "me", "ml", "nj"].contains(&o.str("phylogenies")) {
        return Err("-phylogenies should be 'bionj', 'me', 'ml' or 'nj'".into());
    }
    if steps.contains(&2) {
        let builtin_nj = o.str("phylogenies") == "nj";
        for (k, arg) in [("path_diamond", "version"), ("path_fasttree", "-expert")] {
            if builtin_nj && k == "path_fasttree" { continue }
            let ok = Command::new(o.str(k)).arg(arg)
                .stdin(Stdio::null()).stdout(Stdio::null()).stderr(Stdio::null()).status().is_ok();
            if !ok { return Err(format!("cannot execute '{}' (-{k})", o.str(k)).into()) }
        }
    }
    rayon::ThreadPoolBuilder::new().num_threads(o.get("threads")?).build_global()?;

    for s in steps {
        match s {
            1 => step1::run(&o)?,
            2 => step2::run(&o)?,
            3 => step3::run(&o)?,
            _ => step4::run(&o)?,
        }
    }
    let t = start.elapsed().as_secs();
    println!("\n            Total runtime: {} min {} sec\n", t / 60, t % 60);
    Ok(())
}
