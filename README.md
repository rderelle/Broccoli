[![Cargo Build & Test](https://github.com/rderelle/Broccoli/actions/workflows/build.yml/badge.svg)](https://github.com/rderelle/Broccoli/actions/workflows/build.yml)
[![Clippy check](https://github.com/rderelle/Broccoli/actions/workflows/clippy.yml/badge.svg)](https://github.com/rderelle/Broccoli/actions/workflows/clippy.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/broccoli/README.html)

<p align="center">
  <img width="300" height="auto" src="./images/logo_broccoli.png">
</p>


## Overview

Broccoli is designed to infer orthologous groups and orthologous pairs using a mixed phylogeny-network approach. It also detects chimeric proteins resulting from gene-fusion events and assigns these proteins to the corresponding orthologous groups.

<p align="center">
  <img width="650" height="auto" src="./images/overview_broccoli.png">
</p>

## What's new in v2

Broccoli v2 is a re-implementation of Broccoli v1 (Python) in Rust, written with Claude Opus 5.5. It runs the same four steps with the same parameters and output files, as a single binary, and no longer needs Python or ete3.

Outputs are near-identical to v1: the only differences come from ties, which v1 breaks arbitrarily and v2 deterministically (lowest id or first species); when v1 was modified to break them the same way, both versions gave identical trees, orthologous groups and pairs.

v2 also adds two options that speed up the two most time-consuming parts of the analysis, similarity searches and phylogenetic analyses:
- `-combined_search`: runs a single DIAMOND search of all proteins against one combined database, instead of one search per pair of proteomes, which avoids repeating the database and query setup and makes step 2 faster with near-identical results. Caveats: it uses more memory, and since e-values are computed with the mean proteome size, it should not be used with proteomes of very different sizes.
- `-phylogenies nj`: builds gene trees with a <a href="https://github.com/rderelle/kamino">built-in neighbor-joining</a> instead of FastTree BioNJ (the default, as in v1), which is much faster per tree and removes the need for FastTree. Results are only marginally affected (slightly less accurate).

Efficiency metrics v1 vs v2 with 60 fungal proteomes using 8 CPUs (only steps 1-3; Intel Xeon Platinum 8358 CPU @ 2.60 GHz):


| Version | Options | Runtime (mn) | Peak memory (GB) |
|---------|---------|--------------|------------------|
| v1.4 (Python)   | default | 571             | 36            |
| v2     | default | 419        | 4 |
| v2 | `-combined_search` | 196 | 18 |
| v2 | `-combined_search -phylogenies nj` | 106 | 18 |


## Installation

```bash
conda install bioconda::broccoli
```

or from source:

```bash
git clone https://github.com/rderelle/Broccoli.git
cd Broccoli
cargo build --release           # binary: target/release/broccoli
```

Unless installed via Bioconda, you will also need <a href="https://github.com/bbuchfink/diamond">DIAMOND</a> v0.9.30 or above and <a href="http://www.microbesonline.org/fasttree/">FastTree</a> v2.1.11 or above (FastTree is not needed with `-phylogenies nj`).


## Running Broccoli

All parameters and options are available using the `-help` argument. See the [**documentation**](https://docs.rs/broccoli-rs) for a full description of the input and output files, options and tips.

```bash
# display help menu
broccoli -help

# run Broccoli with 8 threads
broccoli -dir <input_dir> -threads 8

# use one combined DIAMOND search and NJ trees at step 2 (fastest)
broccoli -dir <input_dir> -threads 8 -combined_search -phylogenies nj
```

Broccoli will store the temporary and output files in 4 directories named `dir_step1` to `dir_step4` (one for each step) located in the current directory, or in the directory given with `-output` (e.g. `-output output_broccoli`).

## Citation

If you use Broccoli, please cite:

> Romain Derelle, Hervé Philippe, John K Colbourne. 2020.
> Broccoli: combining phylogenetic and network analyses for orthology assignment.
> [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/37/11/3389/5865275)
