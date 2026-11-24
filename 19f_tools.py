#!/usr/bin/env python3
"""Helpers for the fast DRAFT species tree (19f_fast_draft.sh / 19f_fast_post.sh / 19f_fast_bl.sh); the 19c
helpers (read_fa, write_fa, drop_gap_cols, prune_unrooted) are reused unchanged.
  shrink W MODEL   as 19c 'shrink', for FastTree gene trees (W/gt/<locus>.tre): apply the TreeShrink
                   removals to the trees and to the filtered alignments -> W/final_aln, W/genetrees_shrunk.tre,
                   W/partition.nex (every locus MODEL), W/locus_stats.tsv (+ alignment method per locus)
  concat W         concatenate W/final_aln in partition order (missing taxa = gaps) -> W/concat.faa,
                   W/concat_partition.txt (RAxML-style per-locus ranges), W/concat1.nex (one LG+G4 charset,
                   for dating with fixed branch lengths), W/wastral_topo_unrooted.tre (for FastTree -intree)
  samesplits A B   exit 0 if the unrooted topologies of A and B are identical (RF = 0), else 1"""
import importlib.util, os, sys

R = os.path.expanduser("~/deepsea_conv")
_spec = importlib.util.spec_from_file_location("t19c", f"{R}/scripts/19c_tools.py")
t19c = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(t19c)


def shrink(W, model):
    loci = [l.strip() for l in open(f"{W}/genetrees.list") if l.strip()]
    rm = [l.split() for l in open(f"{W}/treeshrink/output.txt")]
    orig = [l.strip() for l in open(f"{W}/genetrees.tre") if l.strip()]
    assert len(loci) == len(rm) == len(orig), (len(loci), len(rm), len(orig))
    os.makedirs(f"{W}/final_aln", exist_ok=True)
    parts = []
    with open(f"{W}/genetrees_shrunk.tre", "w") as gt, open(f"{W}/locus_stats.tsv", "w") as st, \
         open(f"{W}/treeshrink_removed.tsv", "w") as tr:
        st.write("locus\tn_taxa\tn_removed\tlength\tmodel\taln_method\n")
        tr.write("locus\tcode\n")
        for b, r, t in zip(loci, rm, orig):
            s = t19c.read_fa(f"{W}/filt/{b}.faa")
            for c in r:
                s.pop(c, None); tr.write(f"{b}\t{c}\n")
            s = t19c.drop_gap_cols(s)
            if len(s) < 4:
                continue
            t19c.write_fa(s, f"{W}/final_aln/{b}.faa")
            gt.write(t19c.prune_unrooted(t, r) + "\n")
            L = len(next(iter(s.values())))
            mf = f"{W}/aln/{b}.method"
            meth = open(mf).read().strip() if os.path.exists(mf) else "NA"
            st.write(f"{b}\t{len(s)}\t{len(r)}\t{L}\t{model}\t{meth}\n")
            parts.append(b)
    with open(f"{W}/partition.nex", "w") as o:
        o.write("#nexus\nbegin sets;\n")
        for b in parts:
            o.write(f"  charset {b} = {W}/final_aln/{b}.faa: *;\n")
        o.write("  charpartition mine = " + ", ".join(f"{model}:{b}" for b in parts) + ";\nend;\n")
    print(f"{len(parts)} loci in partition; {sum(len(r) for r in rm)} sequences removed by TreeShrink")


def concat(W):
    from ete3 import Tree
    rows = [l.rstrip("\n").split("\t") for l in open(f"{W}/locus_stats.tsv")][1:]
    loci, models = [r[0] for r in rows], [r[4] for r in rows]
    seqs = {b: t19c.read_fa(f"{W}/final_aln/{b}.faa") for b in loci}
    taxa = sorted(set().union(*[set(s) for s in seqs.values()]))
    t = Tree(f"{W}/wastral_topo.tre")
    tips = set(t.get_leaf_names())
    assert tips == set(taxa), (sorted(tips - set(taxa))[:5], sorted(set(taxa) - tips)[:5])
    cat, pos = {x: [] for x in taxa}, 1
    with open(f"{W}/concat_partition.txt", "w") as pf:
        for b, m in zip(loci, models):
            s = seqs[b]
            L = len(next(iter(s.values())))
            for x in taxa:
                cat[x].append(s.get(x, "-" * L))
            pf.write(f"{m.split('+')[0]}+G4, {b} = {pos}-{pos + L - 1}\n")
            pos += L
    full = {x: "".join(v) for x, v in cat.items()}
    t19c.write_fa(full, f"{W}/concat.faa")
    with open(f"{W}/concat1.nex", "w") as o:
        o.write(f"#nexus\nbegin sets;\n  charset all = {W}/concat.faa: *;\n"
                "  charpartition mine = LG+G4:all;\nend;\n")
    t.unroot()
    t.write(outfile=f"{W}/wastral_topo_unrooted.tre", format=9)
    seen = {}
    for x, v in full.items():
        seen.setdefault(v, []).append(x)
    dup = [g for g in seen.values() if len(g) > 1]
    occ = {x: sum(1 for b in loci if x in seqs[b]) for x in taxa}
    print(f"{len(taxa)} taxa x {pos - 1} columns from {len(loci)} loci; "
          f"loci per taxon min {min(occ.values())} median {sorted(occ.values())[len(occ) // 2]}; "
          f"identical concatenated sequences: {dup if dup else 'none'}")


def samesplits(fa, fb):
    from ete3 import Tree
    r = Tree(fa).robinson_foulds(Tree(fb), unrooted_trees=True)
    print(f"RF {r[0]} / {r[1]}")
    sys.exit(0 if r[0] == 0 else 1)


if __name__ == "__main__":
    if sys.argv[1] == "shrink":
        shrink(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "concat":
        concat(sys.argv[2])
    elif sys.argv[1] == "samesplits":
        samesplits(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown command {sys.argv[1]}")
