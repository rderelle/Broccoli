//! Built-in neighbor-joining (`-phylogenies nj`), adapted from R. Derelle's phylo.rs.
//!
//! Distances: amino-acid p-distances ignoring gaps and unknown residues, corrected for
//! multiple substitutions (LG frequencies). The alignment is bit-sliced (64 columns per word) and
//! constant / all-gap columns are folded away. Single-threaded: Broccoli runs one tree
//! per worker.

const LG_PI: [f64; 20] = [
    0.079_066, 0.055_941, 0.041_977, 0.053_052, 0.012_937, 0.040_767, 0.071_586, 0.057_337,
    0.022_355, 0.062_157, 0.099_081, 0.064_600, 0.022_951, 0.042_302, 0.044_040, 0.061_197,
    0.053_287, 0.012_777, 0.027_843, 0.070_200,
];

const C: f64 = {
    let (mut sum, mut i) = (0.0, 0);
    while i < LG_PI.len() { sum += LG_PI[i] * LG_PI[i]; i += 1 }
    1.0 - sum
};

const GAP: u8 = 20;
const PLANES: usize = 6; // 5 bit planes for the residue code + 1 validity plane
const VALID: usize = 5;

const AA_CODE: [u8; 256] = {
    let mut table = [GAP; 256];
    let order = *b"ARNDCQEGHILKMFPSTWYV";
    let mut i = 0;
    while i < order.len() { table[order[i] as usize] = i as u8; i += 1 }
    table
};

fn corrected_distance(valid: u64, mismatches: u64) -> f64 {
    if valid == 0 { return 0.0 }
    let p = mismatches as f64 / valid as f64;
    -(1.0 - (p / C).min(0.999_999_999)).ln()
}

/// Bit-sliced alignment: for each sequence, `words` blocks of PLANES u64.
struct BitAlignment { words: usize, planes: Vec<u64>, const_cols: u64 }

impl BitAlignment {
    fn new(seqs: &[&[u8]]) -> Self {
        let len = seqs[0].len();
        // columns: all gaps, constant without gaps (valid, never different), or informative
        let mut kept = vec![];
        let mut const_cols = 0;
        for c in 0..len {
            let first = AA_CODE[seqs[0][c] as usize];
            let (mut gap, mut diff) = (first == GAP, false);
            let mut residue = first != GAP;
            for s in &seqs[1..] {
                let x = AA_CODE[s[c] as usize];
                if x == GAP { gap = true } else { residue = true; diff |= x != first }
            }
            if !residue { continue }
            if !gap && !diff { const_cols += 1 } else { kept.push(c) }
        }
        let words = kept.len().div_ceil(64);
        let mut planes = vec![0u64; seqs.len() * words * PLANES];
        for (s, out) in seqs.iter().zip(planes.chunks_mut((words * PLANES).max(1))) {
            for (block, cols) in kept.chunks(64).enumerate() {
                let acc = &mut out[block * PLANES..(block + 1) * PLANES];
                for (bit, &c) in cols.iter().enumerate() {
                    let code = AA_CODE[s[c] as usize];
                    let valid = (code != GAP) as u64;
                    let value = code as u64 * valid;
                    for (plane, slot) in acc[..VALID].iter_mut().enumerate() { *slot |= ((value >> plane) & 1) << bit }
                    acc[VALID] |= valid << bit;
                }
            }
        }
        BitAlignment { words, planes, const_cols }
    }

    /// (sites valid in both sequences, mismatches among them)
    fn counts(&self, i: usize, j: usize) -> (u64, u64) {
        let stride = self.words * PLANES;
        let (a, b) = (&self.planes[i * stride..(i + 1) * stride], &self.planes[j * stride..(j + 1) * stride]);
        let (mut valid, mut mismatches) = (self.const_cols, 0);
        for (x, y) in a.chunks_exact(PLANES).zip(b.chunks_exact(PLANES)) {
            let both = x[VALID] & y[VALID];
            let diff = (x[0] ^ y[0]) | (x[1] ^ y[1]) | (x[2] ^ y[2]) | (x[3] ^ y[3]) | (x[4] ^ y[4]);
            valid += both.count_ones() as u64;
            mismatches += (diff & both).count_ones() as u64;
        }
        (valid, mismatches)
    }
}

/// Full n x n distance matrix (row-major).
fn distances(seqs: &[&[u8]]) -> Vec<f64> {
    let n = seqs.len();
    let ali = BitAlignment::new(seqs);
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        for j in i + 1..n {
            let (v, m) = ali.counts(i, j);
            d[i * n + j] = corrected_distance(v, m);
            d[j * n + i] = d[i * n + j];
        }
    }
    d
}

/// Neighbor-joining tree of the aligned sequences, as newick with branch lengths.
pub fn tree(names: &[u32], seqs: &[&[u8]]) -> String {
    let n = seqs.len();
    let mut dist = distances(seqs);
    let mut sums: Vec<f64> = (0..n).map(|i| dist[i * n..(i + 1) * n].iter().sum()).collect();
    let mut slot: Vec<usize> = (0..n).collect(); // matrix row -> node id (leaves 0..n, then joins)
    let mut joins: Vec<[(usize, f64); 2]> = vec![];
    let mut new_row = vec![0.0; n];
    let mut m = n;
    while m > 2 {
        let factor = (m - 2) as f64;
        // pair minimising Q; ties: lowest pair of node ids (deterministic)
        let (mut best_q, mut best, mut best_key) = (f64::INFINITY, (0, 1), (usize::MAX, usize::MAX));
        for a in 0..m - 1 {
            for b in a + 1..m {
                let q = factor * dist[a * n + b] - sums[a] - sums[b];
                let key = (slot[a].min(slot[b]), slot[a].max(slot[b]));
                if q < best_q || (q == best_q && key < best_key) { (best_q, best, best_key) = (q, (a, b), key) }
            }
        }
        let (a, b) = best;
        let dab = dist[a * n + b];
        let la = (0.5 * dab + (sums[a] - sums[b]) / (2.0 * factor)).max(0.0);
        let lb = (dab - la).max(0.0);
        joins.push([(slot[a], la), (slot[b], lb)]);
        slot[a] = n + joins.len() - 1;

        // row a becomes the new node, row b is replaced by the last row
        let mut new_sum = 0.0;
        for t in (0..m).filter(|&t| t != a && t != b) {
            let (dat, dbt) = (dist[a * n + t], dist[b * n + t]);
            let dut = 0.5 * (dat + dbt - dab);
            new_row[t] = dut;
            sums[t] += dut - dat - dbt;
            new_sum += dut;
        }
        for t in (0..m).filter(|&t| t != a && t != b) {
            dist[a * n + t] = new_row[t];
            dist[t * n + a] = new_row[t];
        }
        dist[a * n + a] = 0.0;
        sums[a] = new_sum;
        let last = m - 1;
        if b != last {
            dist.copy_within(last * n..last * n + m, b * n);
            dist[b * n + b] = 0.0;
            for t in 0..m { dist[t * n + b] = dist[b * n + t] }
            sums[b] = sums[last];
            slot[b] = slot[last];
        }
        m -= 1;
    }
    let len = (dist[1] * 0.5).max(0.0);
    joins.push([(slot[0], len), (slot[1], len)]);

    // newick, from the last join
    fn write(node: usize, n: usize, names: &[u32], joins: &[[(usize, f64); 2]], out: &mut String) {
        if node < n { out.push_str(&names[node].to_string()); return }
        out.push('(');
        for (k, &(child, len)) in joins[node - n].iter().enumerate() {
            if k > 0 { out.push(',') }
            write(child, n, names, joins, out);
            out.push_str(&format!(":{len:.6}"));
        }
        out.push(')');
    }
    let mut out = String::new();
    write(n + joins.len() - 1, n, names, &joins, &mut out);
    out.push(';');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bit_counts_match_naive() {
        let seqs: Vec<&[u8]> = vec![b"AARN-XVA-KKMAAV", b"ARRNAXVA-KKLAAV", b"NNRNA-VR-KKMAAV", b"NNRDA-VR-KKMAAC"];
        let ali = BitAlignment::new(&seqs);
        for i in 0..4 {
            for j in i + 1..4 {
                let (mut v, mut m) = (0, 0);
                for (&a, &b) in seqs[i].iter().zip(seqs[j]) {
                    let (a, b) = (AA_CODE[a as usize], AA_CODE[b as usize]);
                    if a < GAP && b < GAP { v += 1; m += (a != b) as u64 }
                }
                assert_eq!(ali.counts(i, j), (v, m));
            }
        }
    }

    #[test]
    fn groups_close_sequences() {
        let seqs: Vec<&[u8]> = vec![b"AAAAAAAARND", b"AAAAAAARRND", b"RRRRRRRRRNE", b"RRRRRRRNRNE", b"---AAAAARND"];
        let t = crate::tree::Tree::parse(&tree(&[1, 2, 3, 4, 5], &seqs)).unwrap();
        let nwk = t.midpoint_newick();
        assert!(nwk.contains("(3,4)") || nwk.contains("(4,3)"), "{nwk}");
    }
}
