//! Step 3: orthology network and orthologous groups. Orthologous/paralogous pairs
//! are extracted from the trees (species overlap), turned into a weighted network,
//! and communities are found with an asynchronous weighted label propagation,
//! followed by gene-fusion (chimeric protein) detection and spurious-hit removal.

use crate::{create, fresh_dir, lines, Opts, Proteins, R};
use crate::tree::Tree;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io::Write;
use std::sync::atomic::{AtomicU32, Ordering::Relaxed};
use std::sync::Mutex;

struct Edge { to: u32, w: f64, start: u32, end: u32 } // end == 0: no DIAMOND hit from this node

pub fn key(a: u32, b: u32) -> u64 { (a.min(b) as u64) << 32 | a.max(b) as u64 }

const SHARDS: usize = 64; // pair-count maps, each holding a fixed subset of the pairs
fn shard(k: u64) -> usize { ((k >> 32) ^ k) as usize % SHARDS }

/// Leaves of one phylomes.tsv line (tree or similarity group) and the size of the
/// ortho part: leaves of the deepest node passing the species-overlap test. The rest is para.
fn group(line: &str, sp: &[u32], nsp: usize, sp_overlap: f64) -> R<(Vec<u32>, usize)> {
    let mut it = line.splitn(3, '\t');
    let (kind, q, data) = (it.next().unwrap_or(""), it.next().unwrap_or("").parse::<u32>()?, it.next().unwrap_or(""));
    if kind == "S" {
        let ids = data.split(' ').map(|x| x.parse()).collect::<Result<Vec<u32>, _>>()?;
        let n = ids.len();
        return Ok((ids, n));
    }
    let mut good = 1;
    let order = Tree::parse(data)?.walk(q, sp, nsp, |_| true, |s| {
        let (ov, ns, nb) = (s.ov as f64, s.ns as f64, s.nb as f64);
        if s.nb + s.ns - s.ov == 1 || (ov / ns <= sp_overlap && ov / nb <= sp_overlap) {
            good = s.browsed.len() + s.sister.len();
        }
    });
    Ok((order, good))
}

fn merge(mut a: FxHashMap<u64, u32>, mut b: FxHashMap<u64, u32>) -> FxHashMap<u64, u32> {
    if a.len() < b.len() { std::mem::swap(&mut a, &mut b) }
    for (k, v) in b { *a.entry(k).or_default() += v }
    a
}

/// Folds `f(acc, query, target, qstart, qend)` over hits.tsv, in parallel (one accumulator per thread).
fn fold_hits<T: Default + Send>(f: impl Fn(&mut T, u32, u32, u32, u32) + Sync) -> R<Vec<T>> {
    lines("dir_step2/hits.tsv")?.par_bridge().try_fold(T::default, |mut acc, l| -> R<T> {
        let v = l?.split('\t').map(|x| x.parse()).collect::<Result<Vec<u32>, _>>()?;
        f(&mut acc, v[0], v[1], v[2], v[3]);
        Ok(acc)
    }).collect()
}

pub fn run(o: &Opts) -> R<()> {
    let (sp_overlap, min_weight): (f64, f64) = (o.get("sp_overlap")?, o.get("min_weight")?);
    let (min_nb_hits, chim_shared, chim_nb_sp): (usize, f64, usize) = (o.get("min_nb_hits")?, o.get("chimeric_shared")?, o.get("chimeric_nb_sp")?);
    println!(" --- STEP 3: network analysis\n species overlap: {sp_overlap}\n min edge weight: {min_weight}\n min nb hits: {min_nb_hits}\n chimeric edges: {chim_shared}\n chimeric species: {chim_nb_sp}");
    let prots = Proteins::load()?;
    let (sp, nsp, n) = (&prots.sp, prots.files.len(), prots.sp.len());
    fresh_dir("dir_step3")?;
    let mut log = create("dir_step3/log_step3.txt")?;
    let phylomes = "dir_step2/phylomes.tsv";

    // ortho pairs, kept if found at least twice. Counted in shared shards (buffered
    // inserts) rather than one map per thread, so that each pair is stored only once.
    println!(" extract ortho and para");
    let shards: Vec<Mutex<FxHashMap<u64, u32>>> = (0..SHARDS).map(|_| Mutex::default()).collect();
    let flush = |buf: &mut Vec<u64>, s: usize| {
        let mut m = shards[s].lock().unwrap();
        for &k in buf.iter() { *m.entry(k).or_insert(0) += 1 }
        buf.clear();
    };
    lines(phylomes)?.par_bridge().try_fold(|| vec![Vec::new(); SHARDS], |mut bufs: Vec<Vec<u64>>, l| -> R<Vec<Vec<u64>>> {
        let (leaves, good) = group(&l?, sp, nsp, sp_overlap)?;
        let ortho = &leaves[..good];
        for (i, &a) in ortho.iter().enumerate() {
            for &b in &ortho[i + 1..] {
                let (k, s) = (key(a, b), shard(key(a, b)));
                bufs[s].push(k);
                if bufs[s].len() == 1024 { flush(&mut bufs[s], s) }
            }
        }
        Ok(bufs)
    }).try_for_each(|bufs| -> R<()> {
        for (s, mut b) in bufs?.into_iter().enumerate() { flush(&mut b, s) }
        Ok(())
    })?;
    let shards: Vec<FxHashMap<u64, u32>> = shards.into_iter().map(|m| m.into_inner().unwrap()).collect();
    let nb = shards.iter().map(|m| m.values().filter(|&&v| v > 1).count()).sum();
    let mut pairs: FxHashMap<u64, (u32, AtomicU32)> = FxHashMap::with_capacity_and_hasher(nb, Default::default());
    for m in shards { pairs.extend(m.into_iter().filter(|x| x.1 > 1).map(|(k, v)| (k, (v, AtomicU32::new(0))))) }
    // para counts of those pairs (second pass over the trees)
    lines(phylomes)?.par_bridge().try_for_each(|l| -> R<()> {
        let (leaves, good) = group(&l?, sp, nsp, sp_overlap)?;
        let (ortho, para) = leaves.split_at(good);
        for &a in para { for &b in ortho { if let Some(x) = pairs.get(&key(a, b)) { x.1.fetch_add(1, Relaxed); } } }
        Ok(())
    })?;

    // network: edge if ortho > para, weight = ortho / max ortho of the node
    let mut deg = vec![0usize; n];
    for (&k, (nb_o, nb_p)) in &pairs {
        if *nb_o > nb_p.load(Relaxed) { deg[(k >> 32) as usize] += 1; deg[k as u32 as usize] += 1 }
    }
    let mut adj: Vec<Vec<Edge>> = deg.into_iter().map(Vec::with_capacity).collect();
    for (&k, (nb_o, nb_p)) in &pairs {
        if *nb_o > nb_p.load(Relaxed) {
            let (a, b) = ((k >> 32) as u32, k as u32);
            adj[a as usize].push(Edge { to: b, w: *nb_o as f64, start: 0, end: 0 });
            adj[b as usize].push(Edge { to: a, w: *nb_o as f64, start: 0, end: 0 });
        }
    }
    drop(pairs);
    let is_node: Vec<bool> = adj.iter().map(|e| !e.is_empty()).collect();
    let nb_nodes = is_node.iter().filter(|&&x| x).count();
    let nb_edges = adj.iter().map(Vec::len).sum::<usize>() / 2;
    println!(" network: {nb_nodes} nodes, {nb_edges} edges");
    writeln!(log, "#network size:\n{nb_nodes} nodes\n{nb_edges} edges\n\n\n#edge_weight\tnb_edges")?;
    let max: Vec<f64> = adj.iter().map(|e| e.iter().map(|x| x.w).fold(0.0, f64::max)).collect();
    let mut hist = [0usize; 21];
    for (u, edges) in adj.iter_mut().enumerate() {
        edges.retain_mut(|e| {
            let w = e.w / max[u];
            hist[(20.0 * w).round_ties_even() as usize] += 1;
            let keep = w >= min_weight && e.w / max[e.to as usize] >= min_weight;
            e.w = w;
            keep
        });
        edges.sort_by(|a, b| b.w.total_cmp(&a.w).then(a.to.cmp(&b.to)));
    }
    for (i, h) in hist.iter().enumerate() { writeln!(log, "{:?}\t{h}", i as f64 / 20.0)? }
    let removed = nb_edges - adj.iter().map(Vec::len).sum::<usize>() / 2;
    writeln!(log, "\n-> {removed} edges removed\n")?;

    // hit coordinates on network edges (query side)
    println!(" load similarity search outputs");
    let coords = fold_hits(|acc: &mut Vec<_>, q, t, s, e| {
        if adj[q as usize].iter().any(|x| x.to == t) { acc.push((q, t, s, e)) }
    })?;
    for (q, t, s, e) in coords.into_iter().flatten() {
        let x = adj[q as usize].iter_mut().find(|x| x.to == t).unwrap();
        (x.start, x.end) = (s, e);
    }

    // connected components, then LPA + chimeric proteins per component
    let mut comp = vec![false; n];
    let mut ccs: Vec<Vec<u32>> = vec![];
    for s in (0..n).filter(|&i| is_node[i]) {
        if comp[s] { continue }
        comp[s] = true;
        let mut cc = vec![s as u32];
        let mut i = 0;
        while i < cc.len() {
            for e in &adj[cc[i] as usize] { if !comp[e.to as usize] { comp[e.to as usize] = true; cc.push(e.to) } }
            i += 1;
        }
        cc.sort_unstable();
        ccs.push(cc);
    }
    let limit = if 2 * nsp < 10 { 10 } else { 2 * nsp };
    let res: Vec<(Vec<Vec<u32>>, Vec<u32>)> = ccs.par_iter().map(|cc| {
        let nb_sp = cc.iter().map(|&x| sp[x as usize]).collect::<std::collections::HashSet<_>>().len();
        if cc.len() < 4 || nb_sp == 1 { return (vec![cc.clone()], vec![]) }
        let coms = label_propagation(cc, &adj, limit);
        if coms.len() == 1 { return (coms, vec![]) }
        chimeric(coms, &adj, sp, chim_shared, chim_nb_sp)
    }).collect();
    let mut coms: Vec<Vec<u32>> = vec![];
    let mut is_chim = vec![false; n];
    for (c, ch) in res {
        coms.extend(c);
        for x in ch { is_chim[x as usize] = true }
    }
    let nb_chim = is_chim.iter().filter(|&&x| x).count();
    println!(" {} connected components\n {} communities\n {nb_chim} chimeric proteins", ccs.len(), coms.len());
    writeln!(log, "#network analysis:\n{} connected components\n{} communities\n{nb_chim} chimeric proteins", ccs.len(), coms.len())?;

    // spurious hits: proteins with < min_nb_hits DIAMOND hits inside their community
    let mut memb: Vec<Vec<u32>> = vec![vec![]; n];
    for (c, com) in coms.iter().enumerate().filter(|(_, c)| c.len() > min_nb_hits + 1) {
        for &x in com { memb[x as usize].push(c as u32) }
    }
    let hits_in: FxHashMap<u64, u32> = fold_hits(|acc: &mut FxHashMap<u64, u32>, q, t, _, _| {
        for &c in &memb[q as usize] {
            if memb[t as usize].contains(&c) { *acc.entry((q as u64) << 32 | c as u64).or_default() += 1 }
        }
    })?.into_iter().reduce(merge).unwrap_or_default();
    let mut nb_removed = 0;
    for (c, com) in coms.iter_mut().enumerate().filter(|(_, c)| c.len() > min_nb_hits + 1) {
        let before = com.len();
        com.retain(|&x| hits_in.get(&((x as u64) << 32 | c as u64)).copied().unwrap_or(0) as usize >= min_nb_hits);
        nb_removed += before - com.len();
    }
    println!(" {nb_removed} spurious hits removed");

    save_outputs(&prots, &coms, &is_chim, &adj)?;
    println!();
    Ok(())
}

/// Local clustering coefficient, computed on at most `limit` neighbours (lowest ids).
fn lcc(u: u32, adj: &[Vec<Edge>], limit: usize) -> f64 {
    let d = adj[u as usize].len();
    if d < 4 { return 0.0 }
    let mut nb: Vec<u32> = adj[u as usize].iter().map(|e| e.to).collect();
    nb.sort_unstable();
    nb.truncate(limit);
    let max = (nb.len() * (nb.len() - 1)) as f64;
    let found: usize = nb.iter().map(|&v| adj[v as usize].iter().filter(|e| nb.binary_search(&e.to).is_ok()).count()).sum();
    (found as f64 / max * 1e5).round() / 1e5 // found counts each linked pair twice
}

/// Asynchronous weighted label propagation, nodes visited by decreasing lcc.
/// Deterministic: ties go to the label met first (neighbours sorted by weight).
fn label_propagation(cc: &[u32], adj: &[Vec<Edge>], limit: usize) -> Vec<Vec<u32>> {
    let lccs: Vec<f64> = cc.iter().map(|&u| lcc(u, adj, limit)).collect();
    let mut order: Vec<usize> = (0..cc.len()).collect();
    order.sort_by(|&a, &b| lccs[b].total_cmp(&lccs[a]));
    let nodes: Vec<u32> = order.iter().map(|&i| cc[i]).collect();
    let pos: FxHashMap<u32, usize> = nodes.iter().enumerate().map(|(i, &u)| (u, i)).collect();
    let mut label: Vec<usize> = (0..nodes.len()).collect();
    let (mut sum, mut met) = (vec![0.0f64; nodes.len()], vec![]);
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..nodes.len() {
            for e in &adj[nodes[i] as usize] {
                let l = label[pos[&e.to]];
                if sum[l] == 0.0 { met.push(l) }
                sum[l] += e.w;
            }
            let mut best = met[0];
            for &l in &met { if sum[l] > sum[best] { best = l } }
            for &l in &met { sum[l] = 0.0 }
            met.clear();
            if label[i] != best { label[i] = best; changed = true }
        }
    }
    let mut groups: Vec<(usize, Vec<u32>)> = vec![];
    for (i, &u) in nodes.iter().enumerate() {
        match groups.iter_mut().find(|g| g.0 == label[i]) {
            Some(g) => g.1.push(u),
            None => groups.push((label[i], vec![u])),
        }
    }
    groups.into_iter().map(|g| g.1).collect()
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let m = v.len() / 2;
    if v.len() % 2 == 1 { v[m] } else { (v[m - 1] + v[m]) / 2.0 }
}

/// Gene fusions: a protein linked to several communities through non-overlapping
/// regions of its sequence is added to all of them.
fn chimeric(mut coms: Vec<Vec<u32>>, adj: &[Vec<Edge>], sp: &[u32], shared: f64, nb_sp: usize) -> (Vec<Vec<u32>>, Vec<u32>) {
    let limit: Vec<f64> = coms.iter().map(|c| c.len() as f64 * shared).collect();
    let og: FxHashMap<u32, usize> = coms.iter().enumerate().flat_map(|(i, c)| c.iter().map(move |&u| (u, i))).collect();
    let nodes: Vec<(u32, usize)> = coms.iter().enumerate().flat_map(|(i, c)| c.iter().map(move |&u| (u, i))).collect();
    let mut chim = vec![];
    for (u, own) in nodes {
        // neighbours with a DIAMOND hit, grouped by community (in order of appearance)
        let mut found: Vec<(usize, Vec<&Edge>)> = vec![];
        for e in adj[u as usize].iter().filter(|e| e.end != 0) {
            let Some(&c) = og.get(&e.to) else { continue };
            match found.iter_mut().find(|f| f.0 == c) { Some(f) => f.1.push(e), None => found.push((c, vec![e])) }
        }
        if found.len() < 2 { continue }
        found.retain(|(c, es)| {
            let s: std::collections::HashSet<u32> = es.iter().map(|e| sp[e.to as usize]).collect();
            s.len() >= nb_sp && es.len() as f64 >= limit[*c]
        });
        let Some(i) = found.iter().position(|f| f.0 == own) else { continue };
        if found.len() < 2 { continue }
        let span = |es: &[&Edge]| (median(es.iter().map(|e| e.start as f64).collect()), median(es.iter().map(|e| e.end as f64).collect()));
        let (rs, re) = span(&found[i].1);
        for (c, es) in found.iter().filter(|f| f.0 != own) {
            let (ts, te) = span(es);
            let overlap = if ts > rs { re - ts } else { te - rs };
            if overlap <= 0.0 {
                if chim.last() != Some(&u) { chim.push(u) }
                coms[*c].push(u);
            }
        }
    }
    (coms, chim)
}

fn save_outputs(prots: &Proteins, coms: &[Vec<u32>], is_chim: &[bool], adj: &[Vec<Edge>]) -> R<()> {
    let (sp, nsp) = (&prots.sp, prots.files.len());
    let members = prots.members();
    let header: String = prots.files.join("\t");
    let mut f_ogs = create("dir_step3/orthologous_groups.txt")?;
    let mut f_net = create("dir_step3/OGs_in_network.txt")?;
    let mut f_stats = create("dir_step3/statistics_per_OG.txt")?;
    let mut f_counts = create("dir_step3/table_OGs_protein_counts.txt")?;
    let mut f_names = create("dir_step3/table_OGs_protein_names.txt")?;
    writeln!(f_ogs, "#OG_name\tprotein_names")?;
    writeln!(f_stats, "#OG_name\tnb_species\tnb_reduced_prot\tnb_all_prot\tclustering_coefficient")?;
    writeln!(f_counts, "#OG_name\t{header}")?;
    writeln!(f_names, "#OG_name\t{header}")?;

    let mut nb_ogs_per_nb_sp = vec![0usize; nsp + 1];
    let mut assigned = vec![0usize; nsp];
    let mut classified = vec![false; sp.len()];
    let mut chim_ogs: FxHashMap<u32, Vec<String>> = FxHashMap::default();
    let mut c = 0;
    for com in coms {
        let mut per_sp: Vec<Vec<&str>> = vec![vec![]; nsp];
        for &u in com {
            for &m in &members[u as usize] { per_sp[sp[u as usize] as usize].push(&prots.name[m as usize]) }
        }
        let nb_sp = per_sp.iter().filter(|v| !v.is_empty()).count();
        if nb_sp < 2 { continue }
        for &u in com { for &m in &members[u as usize] { classified[m as usize] = true } }
        c += 1;
        let og = format!("OG_{c}");
        nb_ogs_per_nb_sp[nb_sp] += 1;
        for (s, v) in per_sp.iter().enumerate() { assigned[s] += v.len() }
        let all: Vec<u32> = com.iter().flat_map(|&u| members[u as usize].iter().copied()).collect();
        writeln!(f_ogs, "{og}\t{}", all.iter().map(|&m| prots.name[m as usize].as_str()).collect::<Vec<_>>().join(" "))?;
        writeln!(f_net, "{og}\t{}", com.iter().map(|u| u.to_string()).collect::<Vec<_>>().join(" "))?;
        writeln!(f_counts, "{og}\t{}", per_sp.iter().map(|v| v.len().to_string()).collect::<Vec<_>>().join("\t"))?;
        writeln!(f_names, "{og}\t{}", per_sp.iter().map(|v| v.join(" ")).collect::<Vec<_>>().join("\t"))?;
        let set: rustc_hash::FxHashSet<u32> = com.iter().copied().collect();
        let links: usize = com.iter().map(|&u| adj[u as usize].iter().filter(|e| set.contains(&e.to)).count()).sum();
        let cc = links as f64 / (set.len() * (set.len() - 1)) as f64;
        writeln!(f_stats, "{og}\t{nb_sp}\t{}\t{}\t{:?}", com.len(), all.len(), (cc * 1e4).round() / 1e4)?;
        for &u in com.iter().filter(|&&u| is_chim[u as usize]) { chim_ogs.entry(u).or_default().push(og.clone()) }
    }

    let mut f = create("dir_step3/chimeric_proteins.txt")?;
    writeln!(f, "#species_file\tprotein_name\tnb_OG_fused\tlist_fused_OGs")?;
    let mut chim: Vec<_> = chim_ogs.into_iter().collect();
    chim.sort();
    for (u, ogs) in chim {
        writeln!(f, "{}\t{}\t{}\t{}", prots.files[sp[u as usize] as usize], prots.name[u as usize], ogs.len(), ogs.join(" "))?;
    }
    let mut f = create("dir_step3/unclassified_proteins.txt")?;
    for (m, _) in classified.iter().enumerate().filter(|x| !x.1) { writeln!(f, "{}", prots.name[m])? }

    let mut total = vec![0usize; nsp];
    for &s in sp { total[s as usize] += 1 }
    let mut f = create("dir_step3/statistics_per_species.txt")?;
    writeln!(f, "#species\tperc_prot_assigned\tnb_prot_assigned")?;
    for s in 0..nsp {
        let perc = 100.0 * assigned[s] as f64 / total[s].max(1) as f64;
        writeln!(f, "{}\t{:?}\t{}", prots.files[s], (perc * 10.0).round() / 10.0, assigned[s])?;
    }
    let mut f = create("dir_step3/statistics_nb_OGs_VS_nb_species.txt")?;
    writeln!(f, "#nb_species\tnb_OGs")?;
    for (i, v) in nb_ogs_per_nb_sp.iter().enumerate() { writeln!(f, "{i}\t{v}")? }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(edges: &[(u32, u32, f64)], n: usize) -> Vec<Vec<Edge>> {
        let mut adj: Vec<Vec<Edge>> = (0..n).map(|_| vec![]).collect();
        for &(a, b, w) in edges {
            adj[a as usize].push(Edge { to: b, w, start: 0, end: 0 });
            adj[b as usize].push(Edge { to: a, w, start: 0, end: 0 });
        }
        adj
    }

    #[test]
    fn lpa_splits_two_cliques() {
        let mut e = vec![];
        for g in [0u32, 4] { for i in 0..4 { for j in i + 1..4 { e.push((g + i, g + j, 1.0)) } } }
        e.push((3, 4, 0.2));
        let adj = graph(&e, 8);
        let mut coms = label_propagation(&(0..8).collect::<Vec<_>>(), &adj, 10);
        for c in coms.iter_mut() { c.sort() }
        coms.sort();
        assert_eq!(coms, vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]]);
    }
}
