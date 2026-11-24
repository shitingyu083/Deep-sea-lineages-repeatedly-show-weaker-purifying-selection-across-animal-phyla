#!/usr/bin/env python3
"""Helpers for 19c_species_tree.sh.
  frag IN OUT MINFRAC MINLEN MINTAXA   drop sequences with < MINFRAC non-gap columns of the trimmed
                                       alignment, then all-gap columns; exit 2 if the locus ends up
                                       shorter than MINLEN or with < MINTAXA sequences (nothing written)
  shrink W                             apply TreeShrink removals: W/filt -> W/final_aln, write
                                       W/genetrees_shrunk.tre, W/partition.nex, W/locus_stats.tsv
  root IN OUT_WITH OUT_PRUNED          root on the outgroups (data/outgroups/outgroups.tsv); write the
                                       rooted tree with and without them (exit 3 if the ingroup is not
                                       monophyletic)
  compare W A B [TAG]                  RF distance + conflicting splits between trees A and B (split
                                       labels = whatever support each tree carries); W/compare_<TAG>.tsv
  topo IN OUT                          topology-only copy of a Newick tree (for IQ-TREE -te)"""
import csv, os, re, sys

R = os.path.expanduser("~/deepsea_conv")


def read_fa(f):
    seqs, name = {}, None
    for l in open(f):
        if l.startswith(">"):
            name = l[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append("".join(l.split()))
    return {k: "".join(v) for k, v in seqs.items()}


def write_fa(seqs, f):
    with open(f, "w") as o:
        for k, s in seqs.items():
            o.write(f">{k}\n{s}\n")


def drop_gap_cols(seqs):
    if not seqs:
        return seqs
    L = len(next(iter(seqs.values())))
    keep = [i for i in range(L) if any(s[i] not in "-?X" for s in seqs.values())]
    return {k: "".join(s[i] for i in keep) for k, s in seqs.items()}


def frag(fin, fout, minfrac, minlen, mintaxa):
    s = read_fa(fin)
    if not s:
        sys.exit(2)
    L = len(next(iter(s.values())))
    s = {k: v for k, v in s.items() if L and sum(c not in "-?X" for c in v) >= minfrac * L}
    s = drop_gap_cols(s)
    if len(s) < mintaxa or not s or len(next(iter(s.values()))) < minlen:
        sys.exit(2)
    write_fa(s, fout + ".tmp")
    os.replace(fout + ".tmp", fout)


def gt_label(lab):
    """Support of a gene-tree branch label: IQ-TREE writes 'SH-aLRT/aBayes' ('/0.997' with --abayes alone), so
    the part after the last '/'; plain numeric labels (FastTree) are kept; '' = no support (IQ-TREE leaves the
    branches it adds for identical sequences unlabelled)."""
    return lab.rsplit("/", 1)[-1]


def prune_unrooted(txt, drop):
    """Remove tips from an unrooted IQ-TREE gene tree, keeping the aBayes supports (gt_label; unlabelled
    branches stay unlabelled) and summed branch lengths, and leave a trifurcating root (one length/support for
    the edge across the root: the larger of the two supports, a lower bound for the merged split). Every tree
    is rewritten, also without removals, so the labels are plain numbers for wASTRAL -B and ete3 format 0."""
    from ete3 import Tree
    t = Tree(txt, format=1)
    for n in t.traverse():
        if not n.is_leaf():
            n.lab, n.name = gt_label(n.name), ""
    drop = set(drop)
    if drop:
        t.prune([x for x in t.get_leaf_names() if x not in drop], preserve_branch_length=True)
        while len(t.children) == 1:
            t = t.children[0].detach()
        if len(t.children) == 2:
            a, b = t.children
            if a.is_leaf():
                a, b = b, a
            if not a.is_leaf():               # merge the two root edges into one
                b.dist += a.dist
                if not b.is_leaf():
                    s = [float(x) for x in (a.lab, b.lab) if x != ""]
                    b.lab = f"{max(s):g}" if s else ""
                a.delete(prevent_nondicotomic=False)
    t.lab = ""
    return nwk(t)


def shrink(W):
    """TreeShrink was run on W/genetrees.tre (one tree per line, order in W/genetrees.list).
    Removals are applied to the original IQ-TREE trees here (so aBayes supports survive regardless of
    how TreeShrink writes its trees)."""
    loci = [l.strip() for l in open(f"{W}/genetrees.list") if l.strip()]
    rm = [l.split() for l in open(f"{W}/treeshrink/output.txt")]
    orig = [l.strip() for l in open(f"{W}/genetrees.tre") if l.strip()]
    assert len(loci) == len(rm) == len(orig), (len(loci), len(rm), len(orig))
    trees = [prune_unrooted(t, r) for t, r in zip(orig, rm)]
    os.makedirs(f"{W}/final_aln", exist_ok=True)
    with open(f"{W}/genetrees_shrunk.tre", "w") as gt, open(f"{W}/locus_stats.tsv", "w") as st, \
         open(f"{W}/treeshrink_removed.tsv", "w") as tr:
        st.write("locus\tn_taxa\tn_removed\tlength\tmodel\n")
        tr.write("locus\tcode\n")
        parts = []
        for b, r, t in zip(loci, rm, trees):
            s = read_fa(f"{W}/filt/{b}.faa")
            for c in r:
                s.pop(c, None); tr.write(f"{b}\t{c}\n")
            s = drop_gap_cols(s)
            if len(s) < 4:
                continue
            write_fa(s, f"{W}/final_aln/{b}.faa")
            gt.write(t + "\n")
            m = None
            for l in open(f"{W}/gt/{b}.iqtree"):
                if l.startswith("Best-fit model"):
                    m = l.split()[-1]; break
            L = len(next(iter(s.values())))
            st.write(f"{b}\t{len(s)}\t{len(r)}\t{L}\t{m}\n")
            parts.append((b, m))
    with open(f"{W}/partition.nex", "w") as o:
        o.write("#nexus\nbegin sets;\n")
        for b, m in parts:
            o.write(f"  charset {b} = {W}/final_aln/{b}.faa: *;\n")
        o.write("  charpartition mine = " + ", ".join(f"{m}:{b}" for b, m in parts) + ";\nend;\n")
    print(f"{len(parts)} loci in partition; {sum(len(r) for r in rm)} sequences removed by TreeShrink")


def outgroups():
    return [r["code"] for r in csv.DictReader(open(f"{R}/data/outgroups/outgroups.tsv"), delimiter="\t")]


def load(f):
    """Internal-node labels (UFBoot, aBayes, ASTRAL support, IQ-TREE cf labels like '100/97.2/88.1') go
    to the node feature .lab (ete3 3.1 only accepts floats in .support)."""
    from ete3 import Tree
    txt = open(f).read().strip()
    t = Tree(txt, format=1, quoted_node_names=True)
    for n in t.traverse():
        if not n.is_leaf():
            n.lab, n.name = n.name, ""
    return t


def edge_labels(t):
    """Labels are branch attributes: key each by its bipartition (the side without a reference tip) so
    they can be re-attached after re-rooting moves nodes around."""
    leaves = frozenset(t.get_leaf_names()); ref = min(leaves)
    key = lambda n: (lambda s: s if ref not in s else leaves - s)(frozenset(n.get_leaf_names()))
    return key, {key(n): getattr(n, "lab", "") for n in t.traverse() if not n.is_leaf() and not n.is_root()}


def nwk(n, top=True):
    """Newick writer for trees from load(): leaf names, internal labels from .lab, lengths."""
    if n.is_leaf():
        s = n.name
    else:
        s = "(" + ",".join(nwk(c, False) for c in n.children) + ")" + str(getattr(n, "lab", "") or "")
    if not top:
        s += f":{n.dist:.10g}"
    return s + (";" if top else "")


def root(fin, fout, fpruned):
    t = load(fin)
    key, lab = edge_labels(t)
    og = [x for x in outgroups() if t.search_nodes(name=x)]
    if not og:
        sys.exit("no outgroup in tree")
    ing = [x for x in t.get_leaf_names() if x not in og]
    t.set_outgroup(t & ing[0])
    if len(og) > 1 and t.check_monophyly(values=og, target_attr="name", unrooted=True)[0]:
        t.set_outgroup(t.get_common_ancestor(og))
    else:
        t.set_outgroup(t & og[0])
        if len(og) > 1:
            print(f"outgroups {og} not monophyletic; rooted on {og[0]}")
    for n in t.traverse():                 # re-attach every label to its own bipartition
        if not n.is_leaf():
            n.lab = "" if n.is_root() else lab.get(key(n), "")
    # ingroup must be monophyletic once rooted
    mono = t.check_monophyly(values=ing, target_attr="name")[0]
    print(f"rooted on {og}; ingroup monophyletic: {mono}; tips {len(t)}")
    open(fout, "w").write(nwk(t) + "\n")
    t.prune(ing, preserve_branch_length=True)
    while len(t.children) == 1:            # drop the unary root left behind by the outgroup clade
        t = t.children[0].detach()
    t.lab = ""
    open(fpruned, "w").write(nwk(t) + "\n")
    if not mono:
        sys.exit(3)                        # outgroup nested in the ingroup: do not let dating proceed


def topo(fin, fout):
    """Topology-only Newick (ASTRAL/wASTRAL quoted support blocks, node labels and lengths removed)."""
    s = open(fin).read().strip()
    s = re.sub(r"'[^']*'", "", s)
    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r":[-+0-9.eE]+|:nan|:NaN", "", s)
    s = re.sub(r"\)[^,():;]+", ")", s)
    open(fout, "w").write(s + ("\n" if s.endswith(";") else ";\n"))


def splits(t, leaves):
    out = {}
    for n in t.traverse():
        if n.is_leaf() or n.is_root():
            continue
        s = frozenset(n.get_leaf_names())
        if 1 < len(s) < len(leaves) - 1:
            key = min(s, frozenset(leaves) - s, key=lambda x: sorted(x))
            out[key] = getattr(n, "lab", "")
    return out


def compare(W, fa, fb, tag=None):
    a, b = load(fa), load(fb)
    common = sorted(set(a.get_leaf_names()) & set(b.get_leaf_names()))
    a.prune(common); b.prune(common)
    sa, sb = splits(a, common), splits(b, common)
    only_a = [k for k in sa if k not in sb]; only_b = [k for k in sb if k not in sa]
    nrf = (len(only_a) + len(only_b)) / (len(sa) + len(sb))
    tag = tag or f"{os.path.basename(fa).split('.')[0]}_vs_{os.path.basename(fb).split('.')[0]}"
    with open(f"{W}/compare_{tag}.tsv", "w") as o:
        o.write("tree\tsupport_label\tn_tips\tclade\n")
        for k in only_a:
            o.write(f"A_only\t{sa[k]}\t{len(k)}\t{','.join(sorted(k))}\n")
        for k in only_b:
            o.write(f"B_only\t{sb[k]}\t{len(k)}\t{','.join(sorted(k))}\n")
    print(f"{tag}: tips {len(common)}; splits A {len(sa)} B {len(sb)}; RF {len(only_a) + len(only_b)}; nRF {nrf:.4f}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "frag":
        frag(sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))
    elif cmd == "shrink":
        shrink(sys.argv[2])
    elif cmd == "root":
        root(*sys.argv[2:5])
    elif cmd == "compare":
        compare(*sys.argv[2:6])
    elif cmd == "topo":
        topo(*sys.argv[2:4])
