#!/usr/bin/env python3
"""G1 — confirmatory claim C2: cross-clade directional recurrence.

PREREG_G1 section 8, C2 (frozen, verbatim):
    "The cross-clade count of eggNOG-bridged families that recur in the same
     direction in >= 2 independent clades, against the same permutation null
     restricted to bridged families."

IMPLEMENTATION OF THE FROZEN DEFINITION
---------------------------------------
For clade c and family f with K_c maximal deep subtrees, the component sums
s_1..s_K give e_f = #{s_i > 0} and c_f = #{s_i < 0}.  The family's direction in
that clade is

    d_c(f) = +1 if e_f > c_f      (recurrently expanding)
             -1 if c_f > e_f      (recurrently contracting)
              0 if e_f == c_f     (no direction: neutral, counts for neither)

Families are bridged across clades only through eggNOG_OGs, exactly as section 7
prescribes and as already computed into bridged_families.tsv.  Statistic:

    S = #{ OG : max( #{c : d_c(f)=+1}, #{c : d_c(f)=-1} ) >= 2 }

i.e. the number of bridged OGs for which at least two independent clades agree
on the same non-neutral direction.

NULL
----
The PREREG fixes the null as "the same permutation null restricted to bridged
families", i.e. the existing nullp1..nullp10 CAFE5 label-permutation runs.  For
each replicate r, d_c(f) is recomputed from nullp{r}/Base_change.tab using that
tree's own component structure, the bridging is held fixed (bridging is an
annotation property, not a result, so it is not permuted), and S is recomputed.
    p = (1 + #{S_r >= S_obs}) / (1 + 10)

KNOWN LIMITATION, STATED UP FRONT
---------------------------------
g1_kdiagnostic.py measured that the permutation trees do not preserve K (Annelida
obs K=8 vs perm 11..20; Mollusca obs 18 vs perm 26..32).  A larger K moves e_f and
c_f toward ties, which changes how often a direction is called at all.  C2 is
therefore *also* computed against the K-matched structure battery (null1..null5)
for the four clades that have one, and both numbers are reported.  Nothing is
promoted to a claim on the strength of the K-unmatched battery alone.

Read-only over CAFE5 output.  No family filtered, no model refitted.

Writes work/addon/G1/G1_C2.tsv, G1_C2_perOG.tsv, G1_C2.log
"""
import collections
import os

import numpy as np
import pandas as pd
from ete3 import Tree

WORK = os.path.expanduser("~/deepsea_conv/work")
CAFE = os.path.join(WORK, "cafe")
OUT = os.path.join(WORK, "addon", "G1")

CLADES = ["Annelida", "Chondrichthyes", "Mollusca_rescue",
          "Cnidaria", "Echinodermata", "Actinopteri_m_rescue"]

LOG = open(os.path.join(OUT, "G1_C2.log"), "w")


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


def directions(cl, nwk, tab):
    """family -> direction (+1/-1/0), plus K.  Reads one CAFE5 change table."""
    _, im = build_map(nwk)
    comps = deep_components(im)
    K = len(comps)
    df = pd.read_csv(tab, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if set(df.columns) != set(expected_cols(im)):
        raise SystemExit("numbering mismatch %s %s" % (cl, tab))
    cols_by_node = collections.defaultdict(list)
    for c in df.columns:
        cols_by_node[col_node(c)].append(c)
    S = np.empty((len(df), K), dtype=np.float64)
    for j, comp in enumerate(comps):
        cs = [c for n in comp for c in cols_by_node.get(n, [])]
        S[:, j] = df[cs].sum(axis=1).to_numpy(dtype=np.float64)
    e = (S > 0).sum(axis=1)
    c_ = (S < 0).sum(axis=1)
    d = np.sign(e - c_).astype(np.int8)
    return pd.Series(d, index=df.index), K


def cross_clade_S(dirs, bridge):
    """bridge: list of (clade, family, OG).  -> S over the bridged OGs."""
    og_dir = collections.defaultdict(dict)   # og -> {clade: direction}
    for cl, fam, og in bridge:
        d = dirs.get(cl, pd.Series(dtype=np.int8)).get(fam, 0)
        prev = og_dir[og].get(cl)
        # a clade may contribute several families to one OG; keep the strongest
        if prev is None or abs(int(d)) > abs(int(prev)):
            og_dir[og][cl] = int(d)
    S = 0
    per = {}
    for og, cd in og_dir.items():
        npos = sum(1 for v in cd.values() if v > 0)
        nneg = sum(1 for v in cd.values() if v < 0)
        m = max(npos, nneg)
        per[og] = (npos, nneg, len(cd))
        if m >= 2:
            S += 1
    return S, per


def pval(obs, nulls):
    nulls = np.asarray(nulls, dtype=float)
    k = int((nulls >= obs).sum())
    return k, (1 + k) / (1 + len(nulls))


def main():
    say("G1 confirmatory claim C2 — cross-clade directional recurrence")
    say("PREREG_G1 md5 8a19c101207c7016ee1affacb1ee8943 (frozen 2026-11-03)")

    br = pd.read_csv(os.path.join(OUT, "bridged_families.tsv"), sep="\t")
    bridge = list(zip(br["clade"], br["family"].astype(str), br["eggNOG_OG"]))
    say("bridged (clade,family,OG) rows: %d ; distinct OG: %d"
        % (len(bridge), br["eggNOG_OG"].nunique()))

    # ---------------- observation
    dirs_obs, Kobs = {}, {}
    for cl in CLADES:
        d, K = directions(cl, os.path.join(CAFE, cl, "lambda_deep.nwk"),
                          os.path.join(CAFE, cl, "deep", "Base_change.tab"))
        dirs_obs[cl] = d
        Kobs[cl] = K
        say("  %-20s K=%2d  +%-6d -%-6d 0%-6d (of %d bridged)"
            % (cl, K, int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum()),
               int(d.index.isin(
                   br.loc[br.clade == cl, "family"].astype(str)).sum())))
    S_obs, per_obs = cross_clade_S(dirs_obs, bridge)
    say("\nOBSERVED  S = %d  (bridged OGs whose direction is shared by >=2 clades)"
        % S_obs)

    # ---------------- permutation null (PREREG-specified)
    perm_S = []
    for r in range(1, 11):
        dirs = {}
        for cl in CLADES:
            nwk = os.path.join(CAFE, cl, "lambda_nullp%d.nwk" % r)
            tab = os.path.join(CAFE, cl, "nullp%d" % r, "Base_change.tab")
            if not (os.path.exists(nwk) and os.path.exists(tab)):
                continue
            try:
                d, _ = directions(cl, nwk, tab)
                dirs[cl] = d
            except SystemExit as e:
                say("  perm%d %s skipped: %s" % (r, cl, e))
        S_r, _ = cross_clade_S(dirs, bridge)
        perm_S.append(S_r)
        say("  perm%2d  S = %d" % (r, S_r))
    k, p_perm = pval(S_obs, perm_S)
    say("\npermutation battery: obs %d  null %s  p = %d/%d = %.3f"
        % (S_obs, sorted(perm_S), 1 + k, 11, p_perm))

    # ---------------- K-matched structure null (secondary, 4 clades)
    struct_cl = [c for c in CLADES
                 if os.path.exists(os.path.join(CAFE, c, "lambda_null1.nwk"))]
    say("\nK-matched structure battery available for: %s" % ", ".join(struct_cl))
    struct_S = []
    for r in range(1, 6):
        dirs = {}
        for cl in struct_cl:
            nwk = os.path.join(CAFE, cl, "lambda_null%d.nwk" % r)
            tab = os.path.join(CAFE, cl, "null%d" % r, "Base_change.tab")
            if not (os.path.exists(nwk) and os.path.exists(tab)):
                continue
            try:
                d, _ = directions(cl, nwk, tab)
                dirs[cl] = d
            except SystemExit as e:
                say("  struct%d %s skipped: %s" % (r, cl, e))
        # restrict the bridge to the clades this battery covers, both sides
        sub = [b for b in bridge if b[0] in struct_cl]
        S_r, _ = cross_clade_S(dirs, sub)
        struct_S.append(S_r)
        say("  struct%d  S = %d" % (r, S_r))
    sub_dirs = {c: dirs_obs[c] for c in struct_cl}
    S_obs_sub, _ = cross_clade_S(sub_dirs, [b for b in bridge if b[0] in struct_cl])
    k2, p_struct = pval(S_obs_sub, struct_S)
    say("\nstructure battery (%d clades): obs %d  null %s  p = %d/%d = %.3f"
        % (len(struct_cl), S_obs_sub, sorted(struct_S), 1 + k2, 6, p_struct))

    # ---------------- write
    per_rows = [dict(eggNOG_OG=og, n_pos=v[0], n_neg=v[1], n_clades=v[2],
                     shared=int(max(v[0], v[1]) >= 2))
                for og, v in per_obs.items()]
    pd.DataFrame(per_rows).sort_values(
        ["shared", "n_pos", "n_neg"], ascending=False).to_csv(
        os.path.join(OUT, "G1_C2_perOG.tsv"), sep="\t", index=False)

    summary = pd.DataFrame([
        dict(analysis="perm_battery_all_clades", n_clades=len(CLADES),
             S_obs=S_obs, null_sd=float(np.std(perm_S, ddof=1)) if len(perm_S) > 1 else 0.0,
             n_null=len(perm_S), p=p_perm, **{"null_min": min(perm_S), "null_max": max(perm_S)}),
        dict(analysis="struct_battery_Kmatched", n_clades=len(struct_cl),
             S_obs=S_obs_sub,
             null_sd=float(np.std(struct_S, ddof=1)) if len(struct_S) > 1 else 0.0,
             n_null=len(struct_S), p=p_struct,
             **{"null_min": min(struct_S), "null_max": max(struct_S)}),
    ])
    summary.to_csv(os.path.join(OUT, "G1_C2.tsv"), sep="\t", index=False)
    say("\nwrote G1_C2.tsv, G1_C2_perOG.tsv")


if __name__ == "__main__":
    main()
