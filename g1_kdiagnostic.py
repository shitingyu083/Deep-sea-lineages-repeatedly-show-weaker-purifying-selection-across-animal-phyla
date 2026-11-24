#!/usr/bin/env python3
"""G1 K-diagnostic — is the raw-R permutation p interpretable?

PREREG_G1 section 4 fixes R = sum_f max(0, r(f)-1) with r(f) = max(e,c) counted
over K independent deep origins, and tests it against the existing nullp1..10
Base_change.tab files.  R is mechanically increasing in K: two null trees with
the same per-family signals but different K cannot be compared on raw R.

This script measures, for every tree that exists on disk, the quantities that
decide whether the comparison is valid:

    K            number of maximal deep subtrees (one per independent origin)
    n_branch     number of deep-labelled branches
    n_deep_tip   number of deep-labelled tips
    R            PREREG primary statistic, unnormalised
    rate         R / (N_fam * (K-1))   -- the PREREG normalisation
    n_binom      K-invariant recurrence count: families whose e/c split is
                 significantly imbalanced under Binom(K, 0.5), BH q<0.05.
                 E[false positives] = 0 regardless of K, so this column CAN be
                 compared between trees of different K.

Nothing is re-run and no family is filtered.  Read-only over CAFE5 output.

Writes  work/addon/G1/G1_kdiagnostic.tsv and G1_kdiagnostic.log
"""
import collections
import os
import sys

import numpy as np
import pandas as pd
from ete3 import Tree
from scipy.stats import binom

WORK = os.path.expanduser("~/deepsea_conv/work")
CAFE = os.path.join(WORK, "cafe")
OUT = os.path.join(WORK, "addon", "G1")
os.makedirs(OUT, exist_ok=True)

CLADES = ["Annelida", "Chondrichthyes", "Mollusca_rescue",
          "Cnidaria", "Echinodermata", "Actinopteri_m_rescue"]

LOG = open(os.path.join(OUT, "G1_kdiagnostic.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


def build_map(tree_path):
    t = Tree(tree_path, format=1)
    idmap = {}
    for i, nd in enumerate(t.get_leaves(), 1):
        idmap[nd] = i
    k = len(idmap)
    for nd in t.traverse("postorder"):
        if not nd.is_leaf():
            k += 1
            idmap[nd] = k
    return t, idmap


def deep_components(idmap):
    deep = {idmap[nd] for nd in idmap if int(nd.dist) == 2}
    parent = {x: x for x in deep}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for nd in idmap:
        a = idmap[nd]
        if a in deep and nd.up is not None and idmap[nd.up] in deep:
            ra, rb = find(a), find(idmap[nd.up])
            if ra != rb:
                parent[ra] = rb
    comp = collections.defaultdict(set)
    for x in deep:
        comp[find(x)].add(x)
    return list(comp.values())


def col_node(c):
    return int(c[c.rindex("<") + 1:c.rindex(">")])


def expected_cols(idmap):
    return ["%s<%d>" % (nd.name, i) if nd.is_leaf() else "<%d>" % i
            for nd, i in idmap.items()]


def recurrence(df, comps):
    cols_by_node = collections.defaultdict(list)
    for c in df.columns:
        cols_by_node[col_node(c)].append(c)
    sub = pd.DataFrame(index=df.index)
    for j, comp in enumerate(comps):
        cs = [c for n in comp for c in cols_by_node.get(n, [])]
        sub[j] = df[cs].sum(axis=1)
    e = (sub > 0).sum(axis=1).to_numpy()
    c_ = (sub < 0).sum(axis=1).to_numpy()
    r = np.maximum(e, c_)
    return int(np.maximum(0, r - 1).sum()), e, c_, r


def bh_count(p, q=0.05):
    """Number of families rejected by Benjamini-Hochberg at level q."""
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return 0
    o = np.argsort(p)
    ps = p[o]
    thr = q * np.arange(1, n + 1) / n
    below = ps <= thr
    if not below.any():
        return 0
    k = np.max(np.nonzero(below)[0]) + 1
    return int(k)


def binomial_recurrence(e, c_, K):
    """K-invariant statistic: two-sided exact binomial test of e vs c out of K."""
    m = np.maximum(e, c_)
    p = 2.0 * binom.sf(m - 1, K, 0.5)
    np.clip(p, 0.0, 1.0, out=p)
    return bh_count(p), int((p < 0.05 / max(1, p.size)).sum())


def scan(cl, tag, nwk, tab):
    if not (os.path.exists(nwk) and os.path.exists(tab)):
        return None
    _, im = build_map(nwk)
    comps = deep_components(im)
    K = len(comps)
    if K < 2:
        return None
    df = pd.read_csv(tab, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if set(df.columns) != set(expected_cols(im)):
        say("  %s/%s column mismatch, skipped" % (cl, tag))
        return None
    R, e, c_, r = recurrence(df, comps)
    nb = sum(1 for nd in im if int(nd.dist) == 2)
    ntp = sum(1 for nd in im if int(nd.dist) == 2 and nd.is_leaf())
    nb_bh, nb_bonf = binomial_recurrence(e, c_, K)
    return dict(clade=cl, tree=tag, K=K, n_branch=nb, n_deep_tip=ntp,
                n_fam=len(df), R=R, rate=R / (len(df) * (K - 1)),
                exp_R_free=float((len(df) * (K - 1)) / 2.0 ** K),
                n_binom_bh=nb_bh, n_binom_bonf=nb_bonf,
                mean_r=float(r.mean()))


def main():
    rows = []
    for cl in CLADES:
        d = os.path.join(CAFE, cl)
        if not os.path.isdir(d):
            continue
        say("\n" + "=" * 78)
        say("CLADE", cl)
        say("=" * 78)
        trees = [("obs", os.path.join(d, "lambda_deep.nwk"), os.path.join(d, "deep", "Base_change.tab"))]
        for i in range(1, 6):
            trees.append(("struct%d" % i, os.path.join(d, "lambda_null%d.nwk" % i),
                          os.path.join(d, "null%d" % i, "Base_change.tab")))
        for i in range(1, 11):
            trees.append(("perm%d" % i, os.path.join(d, "lambda_nullp%d.nwk" % i),
                          os.path.join(d, "nullp%d" % i, "Base_change.tab")))
        got = []
        for tag, nwk, tab in trees:
            r = scan(cl, tag, nwk, tab)
            if r is not None:
                got.append(r)
        hdr = ("  %-9s %4s %7s %10s %8s %9s %8s %10s %11s"
               % ("tree", "K", "branch", "deep_tip", "R", "rate", "E_R_free", "n_binom", "n_bonf"))
        say(hdr)
        for r in got:
            say("  %-9s %4d %7d %10d %8d %9.4f %8.1f %10d %11d"
                % (r["tree"], r["K"], r["n_branch"], r["n_deep_tip"], r["R"],
                   r["rate"], r["exp_R_free"], r["n_binom_bh"], r["n_binom_bonf"]))
        obs = next(x for x in got if x["tree"] == "obs")
        for bat in ("struct", "perm"):
            ns = [x for x in got if x["tree"].startswith(bat)]
            if not ns:
                say("  %-6s battery: absent" % bat)
                continue
            Ks = sorted({x["K"] for x in ns})
            say("  %-6s battery: n=%d  K range %d..%d  (observed K=%d)  %s"
                % (bat, len(ns), Ks[0], Ks[-1], obs["K"],
                   "K-MATCHED" if Ks == [obs["K"]] else "** K NOT MATCHED **"))
            for stat in ("R", "rate", "n_binom_bh"):
                vs = np.array([x[stat] for x in ns], dtype=float)
                ov = float(obs[stat])
                p = (1 + int((vs >= ov).sum())) / (1 + len(vs))
                say("    %-12s obs %10.4g  null %-40s  p = %.3f"
                    % (stat, ov, np.array2string(np.sort(vs), precision=4, max_line_width=200), p))
        rows.extend(got)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "G1_kdiagnostic.tsv"), sep="\t", index=False)
    say("\nwrote " + os.path.join(OUT, "G1_kdiagnostic.tsv"))


if __name__ == "__main__":
    main()
