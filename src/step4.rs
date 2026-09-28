//! Step 4: orthologous pairs. Within each orthologous group, trees are pruned to
//! the group and every node of the leaf-to-root walk votes ortho or para for the
//! pairs it joins; pairs with ortho / (ortho + para) > ratio_ortho are reported.

use crate::{create, lines, step3::key, Opts, Proteins, R};
use crate::tree::Tree;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::Write;

enum Group { Tree(String), Sim(Vec<u32>) }

pub fn run(o: &Opts) -> R<()> {
    let ratio: f64 = o.get("ratio_ortho")?;
    let not_same_sp = o.str("not_same_sp") == "true";
    println!(" --- STEP 4: orthologous pairs\n ratio ortho: {ratio}\n not same sp: {not_same_sp}");
    let prots = Proteins::load()?;
    let (sp, nsp) = (&prots.sp, prots.files.len());
    let members = prots.members();
    crate::fresh_dir("dir_step4")?;

    let mut groups: FxHashMap<u32, Group> = FxHashMap::default();
    for l in lines("dir_step2/phylomes.tsv")? {
        let l = l?;
        let f: Vec<&str> = l.splitn(3, '\t').collect();
        let g = if f[0] == "S" { Group::Sim(f[2].split(' ').map(|x| x.parse()).collect::<Result<_, _>>()?) } else { Group::Tree(f[2].to_string()) };
        groups.insert(f[1].parse()?, g);
    }
    let ogs: Vec<Vec<u32>> = lines("dir_step3/OGs_in_network.txt")?
        .map(|l| Ok(l?.split('\t').nth(1).unwrap_or("").split(' ').map(|x| x.parse()).collect::<Result<_, _>>()?))
        .collect::<R<_>>()?;
    println!(" analyse {} orthologous groups", ogs.len());

    let mut out = create("dir_step4/orthologous_pairs.txt")?;
    for chunk in ogs.chunks(10_000) { // bounded memory: results are written chunk by chunk
        let res: Vec<Vec<(u32, u32)>> = chunk.par_iter().map(|og| og_pairs(og, &groups, sp, nsp, ratio)).collect::<R<_>>()?;
        for (a, b) in res.into_iter().flatten() {
            if not_same_sp && sp[a as usize] == sp[b as usize] { continue }
            for &x in &members[a as usize] {
                for &y in &members[b as usize] { writeln!(out, "{}\t{}", prots.name[x as usize], prots.name[y as usize])? }
            }
        }
    }
    println!(" done\n");
    Ok(())
}

fn og_pairs(og: &[u32], groups: &FxHashMap<u32, Group>, sp: &[u32], nsp: usize, ratio: f64) -> R<Vec<(u32, u32)>> {
    let inside: FxHashSet<u32> = og.iter().copied().collect();
    let mut cnt: FxHashMap<u64, (u32, u32)> = FxHashMap::default(); // (ortho, para)
    for q in og {
        match groups.get(q) {
            Some(Group::Sim(ids)) => {
                let ids: Vec<u32> = ids.iter().copied().filter(|x| inside.contains(x)).collect();
                for (i, &a) in ids.iter().enumerate() { for &b in &ids[i + 1..] { cnt.entry(key(a, b)).or_default().0 += 1 } }
            }
            Some(Group::Tree(nwk)) => {
                Tree::parse(nwk)?.walk(*q, sp, nsp, |x| inside.contains(&x), |s| {
                    let ortho = s.nb + s.ns - s.ov == 1 || s.ov == 0 || (s.ov == 1 && s.ns - s.ov >= 2 && s.nb - s.ov >= 2);
                    for &a in s.browsed {
                        for &b in s.sister {
                            let c = cnt.entry(key(a, b)).or_default();
                            if ortho { c.0 += 1 } else { c.1 += 1 }
                        }
                    }
                });
            }
            None => {}
        }
    }
    let mut pairs: Vec<(u32, u32)> = cnt.into_iter()
        .filter(|(_, (o, p))| *o > 0 && *o as f64 / (*o + *p) as f64 > ratio)
        .map(|(k, _)| ((k >> 32) as u32, k as u32)).collect();
    pairs.sort_unstable();
    Ok(pairs)
}
