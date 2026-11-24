#!/usr/bin/env python3
"""Helpers for the corrected draft species tree (19g_fixed_bl.sh).
  dropfa IN.faa OUT.faa CODE[,CODE]        copy an alignment without the listed taxa
  prunegt IN.tre OUT.tre CODE[,CODE]       remove taxa from every gene tree (keeps lengths and supports)
  prunetopo IN.tre OUT.tre CODE[,CODE]     remove taxa from one topology (rooting kept, no lengths)
  dropped SPEC                             print the comma-separated drop list of a topo_fix spec"""
import csv, sys
from ete3 import Tree


def dropfa(fin, fout, codes):
    keep = True
    with open(fout, "w") as o:
        for l in open(fin):
            if l.startswith(">"):
                keep = l[1:].split()[0] not in codes
            if keep:
                o.write(l)


def prune_one(t, codes):
    names = set(t.get_leaf_names())
    if not names & codes:
        return t
    t.prune([n for n in names if n not in codes], preserve_branch_length=True)
    return t


def prunegt(fin, fout, codes):
    n = k = 0
    with open(fout, "w") as o:
        for l in open(fin):
            if not l.strip():
                continue
            t = Tree(l.strip(), format=0)
            if set(t.get_leaf_names()) & codes:
                t = prune_one(t, codes); k += 1
            o.write(t.write(format=0) + "\n"); n += 1
    print(f"{n} gene trees, {k} pruned of {sorted(codes)}")


def prunetopo(fin, fout, codes):
    t = prune_one(Tree(open(fin).read().strip(), format=1, quoted_node_names=True), codes)
    t.write(outfile=fout, format=9)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "dropped":
        print(",".join(r["subject"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t") if r["op"] == "drop"))
    else:
        codes = set(sys.argv[4].split(","))
        {"dropfa": dropfa, "prunegt": prunegt, "prunetopo": prunetopo}[cmd](sys.argv[2], sys.argv[3], codes)
