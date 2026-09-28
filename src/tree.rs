//! Minimal newick tree: parsing, midpoint rooting (same rules as ete3) and the
//! leaf-to-root walk used by the species-overlap analyses of steps 3 and 4.

const NONE: usize = usize::MAX;

pub struct Tree {
    parent: Vec<usize>,
    children: Vec<Vec<usize>>,
    name: Vec<Option<u32>>, // leaf protein id
    len: Vec<f64>,          // branch length to parent
}

/// One step of the walk from a leaf to the root: the leaves already browsed,
/// the new sister leaves, and species counts on both sides.
pub struct Step<'a> {
    pub browsed: &'a [u32],
    pub sister: &'a [u32],
    pub nb: usize, // nb species in browsed
    pub ns: usize, // nb species in sister
    pub ov: usize, // nb species in both
}

impl Tree {
    pub fn parse(s: &str) -> Result<Tree, String> {
        let mut t = Tree { parent: vec![NONE], children: vec![vec![]], name: vec![None], len: vec![0.0] };
        let b = s.trim().as_bytes();
        let end = |from: usize, delims: &[u8]| from + b[from..].iter().position(|c| delims.contains(c)).unwrap_or(b.len() - from);
        let (mut cur, mut i) = (0, 0);
        while i < b.len() {
            match b[i] {
                b'(' => cur = t.add(cur),
                b',' => { let p = t.parent[cur]; if p == NONE { return Err(format!("bad newick: {s}")) } cur = t.add(p) }
                b')' => { cur = t.parent[cur]; if cur == NONE { return Err(format!("bad newick: {s}")) } }
                b';' => break,
                b':' => {
                    let j = end(i + 1, b",();");
                    let v = std::str::from_utf8(&b[i + 1..j]).unwrap().trim();
                    t.len[cur] = v.parse().map_err(|_| format!("bad branch length '{v}'"))?;
                    i = j;
                    continue;
                }
                _ => {
                    let j = end(i, b":,();");
                    let tok = std::str::from_utf8(&b[i..j]).unwrap().trim();
                    // leaf names are protein ids; internal labels (supports) are ignored
                    if !tok.is_empty() && t.children[cur].is_empty() {
                        t.name[cur] = Some(tok.parse().map_err(|_| format!("bad leaf name '{tok}'"))?);
                    }
                    i = j;
                    continue;
                }
            }
            i += 1;
        }
        Ok(t)
    }

    fn add(&mut self, p: usize) -> usize {
        let n = self.parent.len();
        self.parent.push(p);
        self.children.push(vec![]);
        self.name.push(None);
        self.len.push(0.0);
        self.children[p].push(n);
        n
    }

    fn is_leaf(&self, n: usize) -> bool { self.children[n].is_empty() }

    /// Farthest leaf below `n` (first one in preorder on ties), distance excluding n's own branch.
    fn farthest_leaf(&self, n: usize) -> (usize, f64) {
        let (mut best, mut stack) = ((n, f64::NEG_INFINITY), vec![(n, 0.0)]);
        if self.is_leaf(n) { return (n, 0.0) }
        while let Some((x, d)) = stack.pop() {
            if x != n && self.is_leaf(x) && d > best.1 { best = (x, d) }
            for &c in self.children[x].iter().rev() { stack.push((c, d + self.len[c])) }
        }
        best
    }

    /// Topology-only newick of this tree rerooted at its midpoint, following ete3's
    /// get_midpoint_outgroup() + set_outgroup() as used by Broccoli v1.
    pub fn midpoint_newick(&self) -> String {
        let (a, _) = self.farthest_leaf(0);
        // ete3 get_farthest_node() starting from leaf a
        let (mut far, mut prev, mut cdist, mut cur) = (0.0, a, self.len[a], self.parent[a]);
        while cur != NONE {
            for &c in &self.children[cur] {
                if c != prev {
                    let d = if self.is_leaf(c) { 0.0 } else { self.farthest_leaf(c).1 } + self.len[c];
                    if cdist + d > far { far = cdist + d }
                }
            }
            prev = cur;
            cdist += self.len[prev];
            cur = self.parent[prev];
        }
        // climb from a until half of the longest path is exceeded
        let (mut d, mut og) = (0.0, a);
        loop {
            if og == NONE { og = self.children[0][0]; break }
            d += self.len[og];
            if d > far / 2.0 { break }
            og = self.parent[og];
        }
        if og == 0 { return self.newick(0, NONE) }
        format!("({},{});", self.newick(og, self.parent[og]), self.newick(self.parent[og], og))
    }

    /// Newick of the subtree hanging from node `n` when coming from its neighbour `from`.
    fn newick(&self, n: usize, from: usize) -> String {
        let nbrs: Vec<usize> = self.children[n].iter().copied()
            .chain(std::iter::once(self.parent[n]))
            .filter(|&x| x != NONE && x != from)
            .collect();
        if nbrs.is_empty() { return self.name[n].map(|x| x.to_string()).unwrap_or_default() }
        if nbrs.len() == 1 && from != NONE { return self.newick(nbrs[0], n) } // collapse unary nodes
        let inner: Vec<String> = nbrs.iter().map(|&x| self.newick(x, n)).collect();
        let s = format!("({})", inner.join(","));
        if from == NONE { s + ";" } else { s }
    }

    /// Kept leaves below node `n`, appended to `out`.
    fn leaves(&self, n: usize, keep: &impl Fn(u32) -> bool, out: &mut Vec<u32>) {
        let mut stack = vec![n];
        while let Some(x) = stack.pop() {
            match self.name[x] {
                Some(id) if self.is_leaf(x) => if keep(id) { out.push(id) },
                _ => stack.extend(&self.children[x]),
            }
        }
    }

    /// Walk from leaf `start` to the root, ignoring leaves not kept, and call `f`
    /// at every node bringing new (sister) leaves. Returns all browsed leaves in
    /// browsing order (`start` first).
    pub fn walk(&self, start: u32, sp: &[u32], nsp: usize, keep: impl Fn(u32) -> bool, mut f: impl FnMut(Step)) -> Vec<u32> {
        let mut browsed = vec![start];
        let Some(mut cur) = (0..self.name.len()).find(|&n| self.name[n] == Some(start)) else { return browsed };
        // flags per species: 1 = in browsed, 2 = in current sister set
        let mut flag = vec![0u8; nsp];
        flag[sp[start as usize] as usize] = 1;
        let (mut nb, mut sister) = (1, vec![]);
        while self.parent[cur] != NONE {
            let p = self.parent[cur];
            sister.clear();
            for &c in &self.children[p] { if c != cur { self.leaves(c, &keep, &mut sister) } }
            if !sister.is_empty() {
                let (mut ns, mut ov) = (0, 0);
                for &x in &sister {
                    let s = &mut flag[sp[x as usize] as usize];
                    if *s & 2 == 0 { *s |= 2; ns += 1; if *s & 1 == 1 { ov += 1 } }
                }
                f(Step { browsed: &browsed, sister: &sister, nb, ns, ov });
                nb += ns - ov;
                for &x in &sister { flag[sp[x as usize] as usize] = 1 }
                browsed.extend_from_slice(&sister);
            }
            cur = p;
        }
        browsed
    }
}

#[cfg(test)]
mod tests {
    use super::Tree;

    #[test]
    fn midpoint_like_ete3() {
        // ete3: (4:0.2,(3:0.3,(1:0.1,2:0.2)0.9:0.05)1:0.2);
        let t = Tree::parse("(1:0.1,2:0.2,(3:0.3,4:0.4)0.9:0.05);").unwrap();
        assert_eq!(t.midpoint_newick(), "(4,(3,(1,2)));");
        let t = Tree::parse(&t.midpoint_newick()).unwrap();
        assert_eq!(t.midpoint_newick(), "(4,(3,(1,2)));"); // no lengths: rooted above first leaf
    }

    #[test]
    fn walk_steps() {
        // species = id / 10
        let sp: Vec<u32> = (0..40).map(|i| i / 10).collect();
        let t = Tree::parse("((((1,11),21),2),31);").unwrap();
        let mut seen = vec![];
        let order = t.walk(1, &sp, 4, |x| x != 21, |s| seen.push((s.sister.to_vec(), s.nb, s.ns, s.ov)));
        assert_eq!(order, vec![1, 11, 2, 31]);
        assert_eq!(seen, vec![(vec![11], 1, 1, 0), (vec![2], 2, 1, 1), (vec![31], 2, 1, 0)]);
    }
}
