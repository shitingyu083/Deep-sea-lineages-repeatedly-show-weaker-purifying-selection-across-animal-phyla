#!/usr/bin/env python3
"""G1 — C2 diagnostic: is the cross-clade "shared direction" excess real?

WHY THIS EXISTS (recorded in DEVIATIONS_g1.md)
----------------------------------------------
PREREG_G1 C2 fixes the null as "the same permutation null restricted to bridged
families".  g1_c2.py ran exactly that: S_obs = 5726 vs permutation null
[4843..5266], p = 1/11 = 0.091 -- the resolution floor of a 10-replicate battery.

That p is NOT usable as evidence, because the permutation null is miscalibrated in
a direction that inflates the observation.  Measured call rates (share of bridged
families that get a non-neutral direction) are, obs vs 10-replicate mean:

    Annelida             91.6 % vs 81.6 %
    Chondrichthyes       97.2 % vs 86.4 %
    Mollusca_rescue      91.7 % vs 86.4 %
    Cnidaria             90.9 % vs 81.4 %
    Echinodermata        94.9 % vs 67.9 %
    Actinopteri_m_rescue 94.4 % vs 91.4 %

A family that is neutral (e = c) cannot contribute to "two clades agree", so a
null with more ties mechanically produces a smaller S.  The permutation trees also
do not preserve K (g1_kdiagnostic.py), which is the same mechanism by another
route: larger K pushes e and c toward ties.

THE TEST RUN HERE
-----------------
A marginal-matched independence null.  For each clade c, take the *observed*
marginal over its bridged families of (prob of +, prob of -, prob of 0).  For each
eggNOG OG, hold fixed the set of clades in which it has a bridged family, and draw
each of those clades' direction independently from that clade's own observed
marginal.  Recompute S.  B = 2000.

This conditions on (i) each clade's call rate, (ii) each clade's expansion /
contraction balance, and (iii) each OG's clade membership, and destroys only the
alignment of directions *between* clades.  That is precisely the C2 hypothesis:
do the same families repeat the same direction across independent origins?

Reporting both is the point: the frozen permutation p says 0.091, this says
whether 0.091 stands up once its own bias is removed.

Writes work/addon/G1/G1_C2_marginal.tsv, G1_C2_marginal.log
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

B = int(os.environ.get("C2_B", "2000"))
SEED = 20261103
LOG = open(os.path.join(OUT, "G1_C2_marginal.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


def build_map(p):
    t = Tree(p, format=1)
    im = {}
    for i, nd in enumerate(t.get_leaves(), 1):
        im[nd] = i
    k = len(im)
    for nd in t.traverse("postorder"):
        if not nd.is_leaf():
            k += 1
            im[nd] = k
    return t, im


def deep_components(im):
    deep = {im[nd] for nd in im if int(nd.dist) == 2}
    par = {x: x for x in deep}

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a

    for nd in im:
        a = im[nd]
        if a in deep and nd.up is not None and im[nd.up] in deep:
            ra, rb = find(a), find(im[nd.up])
            if ra != rb:
                par[ra] = rb
    comp = collections.defaultdict(set)
    for x in deep:
        comp[find(x)].add(x)
    return list(comp.values())


def col_node(c):
    return int(c[c.rindex("<") + 1:c.rindex(">")])


def expected_cols(im):
    return ["%s<%d>" % (nd.name, i) if nd.is_leaf() else "<%d>" % i
            for nd, i in im.items()]


def directions(cl, nwk, tab):
    _, im = build_map(nwk)
    comps = deep_components(im)
    K = len(comps)
    df = pd.read_csv(tab, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if set(df.columns) != set(expected_cols(im)):
        raise SystemExit("numbering mismatch %s %s" % (cl, tab))
    cb = collections.defaultdict(list)
    for c in df.columns:
        cb[col_node(c)].append(c)
    S = np.empty((len(df), K), dtype=np.float64)
    for j, comp in enumerate(comps):
        cs = [c for n in comp for c in cb.get(n, [])]
        S[:, j] = df[cs].sum(axis=1).to_numpy(dtype=np.float64)
    e = (S > 0).sum(axis=1)
    c_ = (S < 0).sum(axis=1)
    return pd.Series(np.sign(e - c_).astype(np.int8), index=df.index), K


def main():
    say("G1 C2 diagnostic — marginal-matched independence null  (B=%d)" % B)

    br = pd.read_csv(os.path.join(OUT, "bridged_families.tsv"), sep="\t")
    br["family"] = br["family"].astype(str)

    # one clade can hold several families for one OG; the bridge is many-to-one.
    # keep the highest-support family per (OG, clade) and use ONLY that one, so
    # each (OG, clade) contributes exactly one direction.
    best = (br.sort_values("support", ascending=False)
              .drop_duplicates(["eggNOG_OG", "clade"]))
    og_clade_fam = collections.defaultdict(dict)   # og -> {clade: family}
    for cl, fam, og in zip(best["clade"], best["family"], best["eggNOG_OG"]):
        og_clade_fam[og][cl] = fam
    cl_fams = {c: set() for c in CLADES}
    for cl, fam in zip(best["clade"], best["family"]):
        cl_fams[cl].add(fam)
    say("bridged OGs %d ; (OG,clade) slots %d ; up to %d families per slot kept"
        % (len(og_clade_fam), sum(len(v) for v in og_clade_fam.values()),
           len(br) - len(best)))

    obs_dirs, Ks, marg = {}, {}, {}
    for cl in CLADES:
        d, K = directions(cl, os.path.join(CAFE, cl, "lambda_deep.nwk"),
                          os.path.join(CAFE, cl, "deep", "Base_change.tab"))
        obs_dirs[cl], Ks[cl] = d, K
        s = d[d.index.isin(cl_fams[cl])]
        n = len(s)
        marg[cl] = np.array([(s > 0).sum(), (s < 0).sum(), (s == 0).sum()],
                            dtype=float) / max(1, n)
        say("  %-20s K=%2d  bridged n=%-6d  marginal  +%.3f  -%.3f  0%.3f"
            % (cl, K, n, marg[cl][0], marg[cl][1], marg[cl][2]))

    # observed direction per (OG, clade) slot
    og_clade_dir = {og: {cl: int(obs_dirs[cl].get(fam, 0))
                         for cl, fam in cd.items()}
                    for og, cd in og_clade_fam.items()}

    S_obs = sum(1 for cd in og_clade_dir.values()
                if max(sum(1 for v in cd.values() if v > 0),
                       sum(1 for v in cd.values() if v < 0)) >= 2)
    say("\nOBSERVED S = %d  over %d bridged OGs" % (S_obs, len(og_clade_dir)))

    # ---- null: per OG, draw each present clade's direction from its own marginal.
    # Vectorised: one (B x nslot) uniform array decides every slot at once via
    # cumulative-marginal cut points, so no Python loop over 19k OGs x B.
    rng = np.random.default_rng(SEED)
    # marg is [+, -, 0]: u < p(+) -> +1 ; u < p(+)+p(-) -> -1 ; else 0
    cuts = {cl: (marg[cl][0], marg[cl][0] + marg[cl][1]) for cl in CLADES}

    ogs = list(og_clade_dir.keys())
    slot_cl = [c for og in ogs for c in og_clade_dir[og]]
    slot_og = np.repeat(np.arange(len(ogs)), [len(og_clade_dir[og]) for og in ogs])
    # index for each slot's categorical draw: shape (nslot, B)
    U = rng.random((len(slot_cl), B))
    V = np.empty((len(slot_cl), B), dtype=np.int8)
    for cl in CLADES:
        idx = [i for i, c in enumerate(slot_cl) if c == cl]
        if not idx:
            continue
        c0, c1 = cuts[cl][0], cuts[cl][1]
        u = U[idx]
        V[idx] = np.where(u < c0, 1, np.where(u < c1, -1, 0)).astype(np.int8)

    npos = np.zeros((len(ogs), B), dtype=np.int16)
    nneg = np.zeros((len(ogs), B), dtype=np.int16)
    np.add.at(npos, slot_og, (V > 0).astype(np.int16))
    np.add.at(nneg, slot_og, (V < 0).astype(np.int16))
    S_null = (np.maximum(npos, nneg) >= 2).sum(axis=0).astype(np.int64)
    say("(null built over %d slots, %d OGs, B=%d)" % (len(slot_cl), len(ogs), B))
    del U, V

    p = (1 + int((S_null >= S_obs).sum())) / (1 + B)
    say("\nmarginal-matched null: mean %.1f  sd %.1f  [%.0f, %.0f]"
        % (S_null.mean(), S_null.std(ddof=1), S_null.min(), S_null.max()))
    say("  obs %d   p = %d/%d = %.4f   z = %.2f"
        % (S_obs, 1 + int((S_null >= S_obs).sum()), B + 1, p,
           (S_obs - S_null.mean()) / (S_null.std(ddof=1) + 1e-12)))

    # ---- pairwise concordance among doubly-called pairs (clean sub-statistic)
    pairs = []
    ognames = list(og_clade_dir.keys())
    for og in ognames:
        cd = og_clade_dir[og]
        cls = [c for c, v in cd.items() if v != 0]
        for i in range(len(cls)):
            for j in range(i + 1, len(cls)):
                pairs.append((cls[i], cd[cls[i]], cls[j], cd[cls[j]]))
    agree = sum(1 for a, va, b, vb in pairs if va == vb)
    # The pair statistic is conditioned on BOTH clades being called, so the
    # expectation must be conditioned the same way.  Under independent marginals:
    #   P(agree | both called) = (p+a*p+b + p-a*p-b) / ((p+a+p-a)*(p+b+p-b))
    exp = 0.0
    for a, va, b, vb in pairs:
        num = marg[a][0] * marg[b][0] + marg[a][1] * marg[b][1]
        den = (marg[a][0] + marg[a][1]) * (marg[b][0] + marg[b][1])
        exp += num / den
    say("\npairwise concordance among doubly-called (OG, clade) pairs:")
    say("  pairs %d  agree %d (%.3f)  expected under marginals %.3f"
        % (len(pairs), agree, agree / max(1, len(pairs)), exp / max(1, len(pairs))))

    pd.DataFrame([dict(statistic="S_shared_dir_ge2", observed=S_obs,
                       null_mean=float(S_null.mean()),
                       null_sd=float(S_null.std(ddof=1)),
                       null_min=int(S_null.min()), null_max=int(S_null.max()),
                       B=B, p=p),
                  dict(statistic="pairwise_concordance", observed=agree / max(1, len(pairs)),
                       null_mean=exp / max(1, len(pairs)), null_sd=float("nan"),
                       null_min=int(pd.Series([0]).iloc[0]), null_max=0,
                       B=B, p=float("nan"))]
                 ).to_csv(os.path.join(OUT, "G1_C2_marginal.tsv"), sep="\t", index=False)
    say("\nwrote G1_C2_marginal.tsv")


if __name__ == "__main__":
    main()
