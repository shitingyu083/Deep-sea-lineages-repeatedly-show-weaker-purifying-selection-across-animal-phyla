#!/usr/bin/env python3
"""Probe CAFE5 node-numbering convention against prep_report.txt ground truth.

CAFE5 change tables name columns as  <species><id>  for tips and  <id>  for internals.
We hypothesise: tips numbered 1..n left-to-right in the newick; internals numbered
n+1..n+m in DFS postorder. If that is right, replaying the deep/shallow label tree
through the same numbering must reproduce prep_report's
    "deep maximal subtrees K (tips T)"
for every clade. Nothing is written; read-only.
"""
import os
import re
import sys

from ete3 import Tree

W = os.path.expanduser("~/deepsea_conv/work/cafe")


def numbering(tree_path):
    """Return (idmap, leaves_in_order). idmap: ete3 node -> CAFE id."""
    t = Tree(tree_path, format=1)
    leaves = t.get_leaves()
    idmap = {}
    for i, nd in enumerate(leaves, 1):
        idmap[nd] = i
    k = len(leaves)
    for nd in t.traverse("postorder"):
        if not nd.is_leaf():
            k += 1
            idmap[nd] = k
    return t, idmap


def label_map(tree_path):
    """node CAFE-id -> 1/2 label, plus the tree."""
    t, idmap = numbering(tree_path)
    lab = {}
    for nd, i in idmap.items():
        lab[i] = int(nd.dist)
    return t, idmap, lab


def deep_components(t, idmap):
    """Maximal connected sets of deep(label==2) branches.

    A 'branch' belongs to the node it leads to (nd.dist). So the branch set of a
    deep component = nodes with label 2. Two nodes are adjacent if parent-child.
    Returns list of node-id lists.
    """
    deep = {idmap[nd] for nd in idmap if int(nd.dist) == 2}
    # union-find over deep nodes linked by parent-child
    parent = {x: x for x in deep}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for nd in idmap:
        if idmap[nd] in deep and nd.up is not None and idmap[nd.up] in deep:
            union(idmap[nd], idmap[nd.up])
    comp = {}
    for x in deep:
        comp.setdefault(find(x), []).append(x)
    return list(comp.values())


def count_tips_in_components(t, idmap, comps):
    rev = {v: k for k, v in idmap.items()}
    out = []
    for c in comps:
        n = 0
        for x in c:
            nd = rev[x]
            # each deep node contributes itself if it is a leaf; plus leaves below
            for lf in nd.get_leaves():
                if idmap[lf] in set(c) or nd is lf:
                    n += 1
            if nd.is_leaf():
                n = 1
                break
        out.append(n)
    return out


def main():
    clades = sys.argv[1:] or ["Annelida"]
    for cl in clades:
        d = os.path.join(W, cl)
        pr = os.path.join(d, "prep_report.txt")
        truth = ""
        if os.path.exists(pr):
            txt = open(pr).read()
            m = re.search(r"deep maximal subtrees (\d+) \(tips (\d+)\)", txt)
            if m:
                truth = "K=%s tips=%s" % (m.group(1), m.group(2))
        t, idmap, lab = label_map(os.path.join(d, "lambda_deep.nwk"))
        comps = deep_components(t, idmap)
        nd_tips = sum(1 for nd in idmap if int(nd.dist) == 2 and nd.is_leaf())
        # total tips under each component (correct way)
        rev = {v: k for k, v in idmap.items()}
        tip_tot = []
        for c in comps:
            s = set(c)
            n = 0
            for x in c:
                nd = rev[x]
                if nd.is_leaf():
                    n += 1
                else:
                    for lf in nd.get_leaves():
                        if idmap[lf] in s:
                            n += 1
            tip_tot.append(n)
        print("%-22s %s | replay: K=%d tips=%d | component tip sizes %s | nodes %d ids %s"
              % (cl, truth, len(comps), nd_tips,
                 sorted(tip_tot, reverse=True), len(idmap),
                 "%d..%d" % (min(idmap.values()), max(idmap.values()))))


if __name__ == "__main__":
    main()
