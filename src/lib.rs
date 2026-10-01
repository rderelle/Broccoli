//! <p align="center">
//!   <img width="300" src="https://raw.githubusercontent.com/rderelle/Broccoli/master/images/logo_broccoli.png" alt="Broccoli">
//! </p>
//!
//! Broccoli is designed to infer orthologous groups and orthologous pairs using a mixed phylogeny-network approach. It also detects 
//! chimeric proteins resulting from gene-fusion events and assigns these proteins to the corresponding orthologous groups.
//!
//! This page is the user manual of the `broccoli` command-line program.
//!
//! <p align="center">
//!   <img width="650" src="https://raw.githubusercontent.com/rderelle/Broccoli/master/images/overview_broccoli.png" alt="Overview of the 4 steps of Broccoli">
//! </p>
//!
//! # Installation
//!
//! The easiest way to install Broccoli is via [Bioconda](https://bioconda.github.io/),
//! which also installs DIAMOND and FastTree:
//!
//! ```bash
//! conda install bioconda::broccoli
//! ```
//!
//! Alternatively, with a [Rust toolchain](https://www.rust-lang.org/tools/install)
//! installed, it can be built from the GitHub repository:
//!
//! ```bash
//! git clone https://github.com/rderelle/Broccoli.git
//! cd Broccoli
//! cargo build --release     # executable: target/release/broccoli
//! ```
//!
//! Broccoli uses two external programs during step 2 (both are installed by Bioconda):
//!
//! - [DIAMOND](https://github.com/bbuchfink/diamond) v0.9.30 or above (similarity
//!   searches). The `-combined_search` option needs a DIAMOND release that supports
//!   `--taxon-k`.
//! - [FastTree](http://www.microbesonline.org/fasttree/) v2.1.11 or above, in its
//!   single-threaded version (phylogenetic analyses; Broccoli runs many FastTree
//!   instances in parallel). FastTree is not needed with `-phylogenies nj`.
//!
//! By default they should be named `diamond` and `fasttree` and be in your `PATH`;
//! otherwise, give their location with `-path_diamond` and `-path_fasttree`.
//!
//! # Quick start
//!
//! The only mandatory argument is the directory containing the proteome files:
//!
//! ```bash
//! # display the help menu
//! broccoli -help
//!
//! # run the whole pipeline on 8 threads
//! broccoli -dir proteomes -threads 8
//!
//! # fastest settings: one combined DIAMOND search and built-in NJ trees
//! broccoli -dir proteomes -threads 8 -combined_search -phylogenies nj
//!
//! # orthologous groups only (steps 1 to 3), with maximum-likelihood trees
//! broccoli -dir proteomes -threads 8 -steps 1,2,3 -phylogenies ml
//!
//! # write the dir_step* directories in output_broccoli/ instead of the current directory
//! broccoli -dir proteomes -threads 8 -output output_broccoli
//! ```
//!
//! Broccoli writes its temporary and output files in 4 directories, `dir_step1` to
//! `dir_step4` (one per step), created in the **current directory**, or in the
//! directory given with `-output` (e.g. `-output output_broccoli`). The main
//! results are:
//!
//! - `dir_step3/orthologous_groups.txt`
//! - `dir_step3/chimeric_proteins.txt`
//! - `dir_step4/orthologous_pairs.txt`
//!
//! # Input data
//!
//! The input is a directory containing one proteome file per species:
//!
//! - sequences in FASTA format;
//! - file extension `.fas`, `.fasta` or `.faa` (any case), optionally gzipped
//!   (e.g. `human.faa.gz`); other files in the directory are ignored;
//! - the protein name is the part of the FASTA header before the first space.
//!   Protein names should be unique across all files: if not, Broccoli prints a
//!   warning and lists the duplicates in `dir_step1/duplicate_names.txt`.
//!
//! Species are identified by their file names in all output tables, so naming the
//! files after the species (e.g. `Homo_sapiens.fasta`) makes the results easier to read.
//!
//! # How Broccoli works
//!
//! The analysis is divided into 4 steps, each writing its files in its own directory:
//!
//! 1. **k-mer simplification** (`dir_step1`): within each proteome, proteins sharing a
//!    k-mer (i.e. an identical fragment) are grouped and only the longest protein of
//!    each group is kept for the next steps. This removes redundancy (isoforms,
//!    identical copies) and speeds up the analysis. The removed proteins are added
//!    back, alongside their representative, in the final outputs.
//! 2. **phylomes** (`dir_step2`): all-against-all similarity searches with DIAMOND
//!    (keeping the best hits in each species), then, for each protein whose hits
//!    contain several proteins of the same species, a phylogenetic tree is built from
//!    the trimmed alignment of its hits.
//! 3. **network analysis** (`dir_step3`): orthologous relationships are extracted from
//!    the trees and combined into an orthology network, in which orthologous groups
//!    are identified by label propagation. Chimeric proteins are then detected and
//!    assigned to all their orthologous groups, and spurious hits are removed.
//! 4. **orthologous pairs** (`dir_step4`): within each orthologous group, pairs of
//!    proteins are classified as orthologs or paralogs from their relationships in
//!    all the trees they appear in.
//!
//! Two things to keep in mind:
//!
//! - each step requires the directories generated by the previous steps, so steps
//!   can be run separately (e.g. `-steps 1,2,3` then later `-steps 4`), but always
//!   from the same directory, or with the same `-output`;
//! - each step first deletes its own directory if it exists, to avoid mixing files
//!   from different runs. Rename or move the results you want to keep before
//!   re-running a step, or use a different `-output` for each analysis.
//!
//! # Output files
//!
//! Orthologous groups are named `OG_1`, `OG_2`, ... and contain proteins from at
//! least 2 species. Every step also writes a log file (`log_step1.txt` to
//! `log_step3.txt`) in its directory.
//!
//! ## Orthologous groups and chimeric proteins (`dir_step3`)
//!
//! | File | Content |
//! |------|---------|
//! | `orthologous_groups.txt` | one orthologous group per line: its name, then the names of its proteins separated by spaces |
//! | `table_OGs_protein_counts.txt` | number of proteins of each species (columns) in each orthologous group (rows) |
//! | `table_OGs_protein_names.txt` | same table with the protein names instead of their number |
//! | `chimeric_proteins.txt` | chimeric proteins: species file, protein name, number and list of the orthologous groups it belongs to |
//! | `unclassified_proteins.txt` | proteins not assigned to any orthologous group |
//! | `statistics_per_OG.txt` | for each orthologous group: number of species, number of proteins after the k-mer simplification and in total, clustering coefficient |
//! | `statistics_per_species.txt` | percentage and number of proteins of each species assigned to an orthologous group |
//! | `statistics_nb_OGs_VS_nb_species.txt` | number of orthologous groups containing 2, 3, ... species |
//! | `log_step3.txt` | size of the orthology network, distribution of the edge weights, numbers of communities and chimeric proteins |
//!
//! `OGs_in_network.txt` is an intermediate file used by step 4.
//!
//! ## Orthologous pairs (`dir_step4`)
//!
//! `orthologous_pairs.txt` lists one pair of orthologous proteins per line (two
//! protein names separated by a tab).
//!
//! ## Other files
//!
//! - `dir_step1/log_step1.txt`: number of proteins in each proteome before and after
//!   the k-mer simplification.
//! - `dir_step2/log_step2.txt`: for each species, number of proteins with a
//!   phylogeny, without a phylogeny (no duplicated species among their hits), with an
//!   empty alignment after trimming, and with a failed tree.
//! - the other files in `dir_step1` and `dir_step2` (reduced proteomes, similarity
//!   hits, trees) are intermediate files used by the next steps.
//!
//! # Options
//!
//! Options are given with a single dash (e.g. `-threads 8`); `-help` (or `-h`)
//! prints the list of options with their default values.
//!
//! ## General options
//!
//! | Option | Default | Description |
//! |--------|---------|-------------|
//! | `-steps` | `1,2,3,4` | steps to perform, comma separated, consecutive |
//! | `-threads` | `1` | number of threads |
//! | `-output`, `-o` | current directory | directory where the `dir_step*` directories are written (created if needed) |
//!
//! `-steps` selects which steps are performed. For instance, if you only need
//! orthologous groups, use `-steps 1,2,3`. The steps must be consecutive
//! (`-steps 1,3` is not valid).
//!
//! ## Step 1: k-mer simplification
//!
//! | Option | Default | Description |
//! |--------|---------|-------------|
//! | `-dir` | *required* | directory containing the proteome files (full or relative path) |
//! | `-kmer_size` | `100` | length of the k-mers |
//! | `-kmer_min_aa` | `15` | minimum number of different amino acids a k-mer should contain |
//!
//! Proteins of the same proteome sharing at least one k-mer of length `-kmer_size` are
//! grouped, and only the longest one is kept. K-mers made of fewer than `-kmer_min_aa`
//! different amino acids, or containing `X` or `*`, are ignored to avoid grouping
//! proteins through low-complexity or repetitive regions.
//!
//! ## Step 2: phylomes
//!
//! | Option | Default | Description |
//! |--------|---------|-------------|
//! | `-path_diamond` | `diamond` | path of the DIAMOND executable |
//! | `-path_fasttree` | `fasttree` | path of the FastTree executable |
//! | `-e_value` | `0.001` | e-value threshold of the similarity searches |
//! | `-nb_hits` | `6` | maximum number of hits per species |
//! | `-max_gap` | `0.7` | maximum fraction of gaps per alignment position |
//! | `-phylogenies` | `bionj` | tree-building method: `bionj`, `me`, `ml` or `nj` |
//! | `-combined_search` | off | one DIAMOND search against a combined database |
//!
//! `-e_value` and `-nb_hits` correspond to the DIAMOND options `--evalue` and
//! `--max-target-seqs` (applied per species). `-max_gap` controls the trimming of the
//! alignments: positions with a fraction of gaps equal to or above this value are removed.
//!
//! `-phylogenies` sets how gene trees are built:
//!
//! - `bionj`: FastTree BioNJ (fast);
//! - `me`: FastTree minimum evolution (slower, more precise);
//! - `ml`: FastTree maximum likelihood (much slower, most precise);
//! - `nj`: built-in neighbor-joining, much faster per tree than FastTree and does not
//!   require FastTree; results are only marginally less accurate than with `bionj`.
//!
//! By default, Broccoli runs one DIAMOND search per pair of proteomes. With
//! `-combined_search`, it runs a single search of all proteins against one database
//! containing all proteomes, which avoids repeating the database and query setup and
//! makes step 2 much faster, with near-identical hits. Two caveats: it uses more
//! memory, and since e-values are computed with the mean proteome size, it should not
//! be used with proteomes of very different sizes.
//!
//! ## Step 3: network analysis
//!
//! | Option | Default | Description |
//! |--------|---------|-------------|
//! | `-sp_overlap` | `0.5` | maximum ratio of overlapping species in phylogenetic trees |
//! | `-min_weight` | `0.1` | minimum weight for an edge to be kept in the orthology network |
//! | `-min_nb_hits` | `2` | spurious hits: minimum number of hits within the orthologous group |
//! | `-chimeric_shared` | `0.5` | chimeric proteins: minimum fraction of connected proteins in each orthologous group |
//! | `-chimeric_nb_sp` | `3` | chimeric proteins: minimum number of species connected in each orthologous group |
//!
//! `-sp_overlap` controls how orthologous groups are delineated in the trees. At each
//! node, the ratio is the maximum fraction of species shared between the leaves
//! already browsed and the sister leaves. Low values (e.g. 0.3) produce small, tightly
//! defined orthologous groups, while high values (e.g. 0.7) produce larger groups,
//! similar to the orthogroups of OrthoFinder.
//!
//! `-min_weight` cleans the orthology network by removing weak edges (the weight of an
//! edge is relative to the strongest edge of each protein). Besides removing spurious
//! connections, this speeds up step 3.
//!
//! The 3 remaining options control the two corrections applied after the
//! identification of orthologous groups:
//!
//! - **spurious hits**: a protein is removed from its orthologous group if it has
//!   fewer than `-min_nb_hits` DIAMOND hits with other proteins of the group (only
//!   applied to groups containing more than `-min_nb_hits` + 1 proteins);
//! - **chimeric proteins**: a protein connected to several orthologous groups through
//!   non-overlapping regions of its sequence is considered chimeric (gene fusion) if,
//!   for each group, it is connected to at least a fraction `-chimeric_shared` of the
//!   group's proteins, belonging to at least `-chimeric_nb_sp` species. It is then
//!   assigned to all these groups.
//!
//! ## Step 4: orthologous pairs
//!
//! | Option | Default | Description |
//! |--------|---------|-------------|
//! | `-ratio_ortho` | `0.5` | minimum ratio ortho / (ortho + para) to report a pair |
//! | `-not_same_sp` | off | ignore orthologous relationships between proteins of the same species |
//!
//! For each pair of proteins of an orthologous group, Broccoli counts how many times
//! they appear as orthologs or as paralogs across the trees. `-ratio_ortho` sets the
//! precision/recall trade-off: increasing it improves the precision of the
//! orthologous pairs, decreasing it improves their recall.
//!
//! `-not_same_sp` discards pairs of proteins from the same species (e.g. for the
//! Quest for Orthologs benchmark).
//!
//! # Tips and FAQ
//!
//! ### How to make Broccoli faster?
//!
//! The most time-consuming parts of the analysis are the similarity searches and
//! the phylogenetic analyses (step 2):
//!
//! - `-combined_search` speeds up the similarity searches (see the caveats in
//!   [step 2](#step-2-phylomes));
//! - `-phylogenies nj` speeds up the phylogenetic analyses;
//! - reducing `-nb_hits` reduces the number of sequences in alignments and trees, at
//!   a potential cost in precision;
//! - reducing `-kmer_size` simplifies the proteomes further at step 1, at the risk of
//!   grouping paralogs together. This has nearly no effect on prokaryotic
//!   proteomes, and in some species many unrelated proteins share identical fragments
//!   (e.g. viral domains in *Trichomonas vaginalis*).
//!
//! As an example, with 60 fungal proteomes on 8 threads (steps 1 to 3; Intel Xeon
//! Platinum 8358 @ 2.60 GHz):
//!
//! | Options | Runtime (min) | Peak memory (GB) |
//! |---------|---------------|------------------|
//! | default | 419 | 4 |
//! | `-combined_search` | 196 | 18 |
//! | `-combined_search -phylogenies nj` | 106 | 18 |
//!
//! ### How to make Broccoli more precise?
//!
//! Minimum-evolution or maximum-likelihood trees (`-phylogenies me` or
//! `-phylogenies ml`) improve the delineation of closely related orthologous groups,
//! but increase the runtime of the phylogenetic analyses by a factor of about 30 and
//! 100 respectively.
//!
//! Lowering `-sp_overlap` (e.g. 0.4) improves precision but lowers recall (small
//! orthologous groups missing some divergent in-paralogs), which can be suitable for
//! the identification of phylogenomic markers.
//!