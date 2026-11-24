#!/usr/bin/env python3
"""Inputs for RERconverge: per-locus alignments and the species topology pruned to each locus's taxa.
Usage: python 36_rer_prep.py SRC TOPO OUT [MIN_TAXA=10]
  SRC   species-tree work dir with final_aln/<locus>.faa and locus_stats.tsv (post-TreeShrink alignments)
  TOPO  rooted species topology, tips = species codes (e.g. work/draft_v2/topo.tre)
  OUT   writes aln/<locus>.faa, topo/<locus>.tre (unrooted, no branch lengths) and jobs.tsv (locus, n_taxa, model)
Sequences of species absent from TOPO (excluded taxa, e.g. Euncla in meta/topo_fix.tsv) are dropped; loci with
fewer than MIN_TAXA remaining species are skipped (RERconverge uses genes with >= 10 species)."""
import csv, os, sys
from ete3 import Tree

def fasta(fn):
    name, seq = None, []
    for line in open(fn):
        line = line.strip()
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name:
        yield name, "".join(seq)

def main():
    src, topo, out = sys.argv[1:4]
    mint = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    for d in ("aln", "topo"):
        os.makedirs(f"{out}/{d}", exist_ok=True)
    sp = Tree(open(topo).read().strip(), format=1, quoted_node_names=True)
    tips = set(sp.get_leaf_names())
    rows = list(csv.DictReader(open(f"{src}/locus_stats.tsv"), delimiter="\t"))
    jobs, skipped, dropped = [], [], {}
    for r in rows:
        loc = r["locus"]
        seqs = [(n, s) for n, s in fasta(f"{src}/final_aln/{loc}.faa")]
        for n, _ in seqs:
            if n not in tips:
                dropped[n] = dropped.get(n, 0) + 1
        seqs = [(n, s) for n, s in seqs if n in tips]
        if len(seqs) < mint:
            skipped.append(loc)
            continue
        with open(f"{out}/aln/{loc}.faa", "w") as f:
            f.writelines(f">{n}\n{s}\n" for n, s in seqs)
        t = sp.copy()
        t.prune([n for n, _ in seqs])
        while len(t.children) == 1:          # a root left with one child after pruning
            t = t.children[0].detach()
        if len(t.children) == 2:
            t.unroot()
        t.write(outfile=f"{out}/topo/{loc}.tre", format=9)
        jobs.append((loc, len(seqs), r["model"] or "LG+G4"))
    with open(f"{out}/jobs.tsv", "w") as f:
        f.write("locus\tn_taxa\tmodel\n")
        f.writelines(f"{a}\t{b}\t{c}\n" for a, b, c in jobs)
    if dropped:
        print("dropped sequences of species not in the topology:",
              ", ".join(f"{k} ({v} loci)" for k, v in sorted(dropped.items())))
    if skipped:
        print(f"skipped {len(skipped)} loci with < {mint} species: " + ",".join(skipped))
    print(f"{len(jobs)} loci prepared from {len(rows)} (topology {len(tips)} tips)")

if __name__ == "__main__":
    main()
