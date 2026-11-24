#!/usr/bin/env python3
"""G1 exploratory E4 -- do gene-family changes concentrate on early or late
deep-sea branches?

PREREG_G1 section 8 already registers this as exploratory, verbatim:
    "E4 Timing: whether recurring changes concentrate on early versus late deep
     branches."

No new pre-registration is written and no frozen file is edited.  Every number
below carries the label "exploratory, not corrected for multiple testing,
hypothesis-generating" and none of it can become a headline (PREREG_G1 s8).

v2 -- WHAT WAS WRONG IN THE FIRST RUN (2026-11-03)
--------------------------------------------------
The first run computed branch age as the root-to-DISTAL-end distance, i.e. the
depth of the node the branch leads to:

    a_b = dist(node) / height                     # first run, discarded

On a strictly ultrametric tree every terminal branch then has a_b exactly 1.0,
whatever its real timing, so the statistic collapses into "terminal branch or
internal branch".  Echinodermata has 6 deep branches, all six terminal: its ages
all came out 1.000, A_obs = -0.00000, N1 null sd = 0.00000, and N2 (six
singleton components) is the identity permutation, which is why its p was
exactly 1.0000.  That is a property of the tree plus the formula, not biology.

The docstring of the first run already said PROXIMAL end; the code did not.  v2
uses the proximal end, which is what the statistic needs:

    a_b = dist(parent(node)) / height             # v2

a_b is then the age at which the branch began, in [0, 1): 0 = begins at the
root, -> 1 = begins just before a tip.  Terminal and internal branches are
treated identically.  The first run's outputs are kept only as
E4_timing_v0_DISCARDED_distal.* and are not reported.

WHY G1's OWN NULLS CANNOT BE REUSED HERE
-----------------------------------------
Both G1 null batteries change the deep-branch set.  Measured on Annelida: the
observed deep tree carries 34 deep branches, the nullp battery 22-31 over its 40
replicates, the structure-matched battery 22 over its 5.  Under either one, the
branches whose ages get compared are a different set of branches, so the age
statistic would not be a statistic about the observed clade at all.

E4 is also not asking G1's question again.  G1 C1 asked how many independent
origins expand/contract and found a clean null.  E4 asks a conditional question:
given the deep branches this clade has, is the change on them spread evenly in
time or piled at one end?  The null that answers it holds the deep-branch set
fixed and permutes AGE among those branches.  That null is fixed here, before
the run, and no alternative null is computed afterwards.

STATISTIC
---------
    w_b   = sum over the included families of |Delta| on branch b
    A_obs = sum_b w_b a_b / sum_b w_b  -  mean_b a_b

A_obs > 0 means the change mass sits on LATER branches, < 0 on EARLIER.  The
second term is the same quantity with weights ignored, so A_obs is the part of
the weighted mean age that the weights explain.  Two-sided permutation p:
|A_null - mean A_null| >= |A_obs - mean A_null|.

Readouts (fixed in advance):
  ALL   -- w_b over every family, i.e. turnover of any kind.
  RECUR -- w_b over families with r >= 2 only, i.e. the RECURRING changes,
           which is the form the question is phrased in.

Nulls for the deep arm (fixed in advance):
  N1 global      -- ages permuted across all deep branches of the clade.
  N2 within-comp -- ages permuted only within each origin component, so the
                    component-level age structure is preserved.  N2 is the more
                    conservative of the two when branches inside a component
                    share a rate; it is reported beside N1 rather than chosen
                    after the fact.  N2 is undefined (p = 1 by construction) for
                    a clade whose deep branches are all singleton components;
                    the script reports the number of non-singleton components so
                    this is visible.

N2 IS INVALID AS A TEST -- MEASURED, NOT ASSUMED (added after the v2 run)
-----------------------------------------------------------------------
N2 permutes ages within an origin component.  A component is a chain: the deep
branches on one descending lineage.  Their proximal ages are therefore monotonically
increasing along the chain by construction, and a within-component permutation
destroys a constraint the tree guarantees.  Measured on all clades with a
non-singleton component: 35 of 35 chains have strictly increasing proximal ages
(Annelida 6/6, Chondrichthyes 4/4, Mollusca_rescue 6/6, Cnidaria 5/5,
Echinodermata 0 chains, Actinopteri_m_rescue 8/8).

The consequence is visible in the v2 numbers: N2's null mean sits at about -0.08
rather than 0, because shuffled-within-chain ages are systematically younger than
the observed ones, so every clade's z turns spuriously positive.  N2 is retained
in the table as a diagnostic and is NOT used for inference; N1 is the test, and
every N2 z is reported with the label "null invalid -- see header".  This is
recorded so that a later reader cannot pick N2 because it gives a bigger number.

SHALLOW-BRANCH CONTROL (fixed in advance)
-----------------------------------------
The same statistic is computed over the clade's SHALLOW branches, with N1 age
permutation among them.  This is the control I said was required before write-up:
if shallow branches show the same late-concentration, a deep-specific timing
interpretation is not supported and what is being measured is a property of how
CAFE allocates turnover near tips, not a deep-sea timing result.  The control is
reported beside the deep arm in every table; it is not a denominator and the
deep number is not divided by it.

OUTCOME OF THE TOPOLOGY CONTROL (added after the v2 run, before write-up)
------------------------------------------------------------------------
81-98% of the deep change mass sits on deep branches that lead DIRECTLY to a tip
(Annelida 42982/52864, Chondrichthyes 18118/20281, Mollusca 142034/164243,
Cnidaria 82878/96542, Echinodermata 21403/21403, Actinopteri 153263/155939).
Restricting the same statistic to deep branches that lead to INTERNAL nodes --
where a late start cannot be produced by tip topology alone -- leaves no signal in
any clade, and both clades that were nominally significant under N1 in v2 change
SIGN: Mollusca A -0.06121 (p=0.0173) -> +0.07902 (p=0.0911), Cnidaria A -0.06630
(p=0.0223) -> +0.00949 (p=0.6066).  Echinodermata has no internal-leading deep
branch at all, so the control is undefined there.

So E4 does not find a timing signal.  What it measures is where the deep-node
labels fall relative to tips, not when deep-sea turnover happened.  The shallow
control is flat in every clade (|z| <= 1.26), which rules out a general
CAFE-tip-turnover explanation but does not rescue the deep number: the
topology-restricted deep arm already is the honest test, and it is null.

WHY A FIXED AGE PERMUTATION IS THE RIGHT TEST, AND ITS LIMIT
------------------------------------------------------------
Under "turnover is independent of when the branch began", a_b is exchangeable
over the clade's branches, so permuting it is exact.  The limit is stated now:
branches are not evolutionarily independent, and if a fast-evolving subclade
both starts late and turns over more, the permutation null reads that as a
timing effect.  N2 partially controls for this on the deep arm; the shallow arm
has no component structure and only N1.

Usage:  ./envs/dc/bin/python scripts/g1_e4_timing.py
Env:    E4_B (10000)
Writes  work/addon/G1/E4_timing.tsv, E4_timing_perclade.tsv, E4_timing.log
"""
import collections
import os

import numpy as np
import pandas as pd
from ete3 import Tree

WORK = os.path.expanduser("~/deepsea_conv/work")
CAFE = os.path.join(WORK, "cafe")
OUT = os.path.join(WORK, "addon", "G1")
os.makedirs(OUT, exist_ok=True)

CLADES = ["Annelida", "Chondrichthyes", "Mollusca_rescue",
          "Cnidaria", "Echinodermata", "Actinopteri_m_rescue"]
B = int(os.environ.get("E4_B", "10000"))
SEED = 20261103
LOG = open(os.path.join(OUT, "E4_timing.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


# ---------------------------------------------------------------- tree handling
def build_map(p):
    """Verbatim from g1_family_convergence.py -- CAFE5 numbering replay."""
    t = Tree(p, format=1)
    leaves = t.get_leaves()
    im = {}
    for i, nd in enumerate(leaves, 1):
        im[nd] = i
    k = len(leaves)
    for nd in t.traverse("postorder"):
        if not nd.is_leaf():
            k += 1
            im[nd] = k
    return t, im


def deep_components(im):
    """Verbatim from g1_family_convergence.py."""
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
    return [comp[k] for k in comp]


def col_node(c):
    return int(c[c.rindex("<") + 1:c.rindex(">")])


def recurrence(df, comps):
    """Verbatim statistic of g1_family_convergence.recurrence (PREREG_G1 s4)."""
    node_of = {c: col_node(c) for c in df.columns}
    cols_by_node = collections.defaultdict(list)
    for c, n in node_of.items():
        cols_by_node[n].append(c)
    sub = pd.DataFrame(index=df.index)
    for j, comp in enumerate(comps):
        cs = [c for n in comp for c in cols_by_node.get(n, [])]
        sub[j] = df[cs].sum(axis=1)
    e = (sub > 0).sum(axis=1)
    c_ = (sub < 0).sum(axis=1)
    r = np.maximum(e, c_)
    R = int(np.maximum(0, r - 1).sum())
    fam = pd.DataFrame({"e": e, "c": c_, "r": r})
    return R, fam


# ---------------------------------------------------------------- the statistic
def age_stat(w, a):
    """Weighted mean age minus the unweighted mean age.  w, a in a fixed order."""
    sw = w.sum()
    if sw <= 0:
        return np.nan
    return float((w * a).sum() / sw - a.mean())


def perm_null(w, a, comp_idx, mode, B, rng):
    """Array of B null A values.

    mode 'global'      : permute a across all branches.
    mode 'within_comp' : permute a within each component separately.
    """
    nb = len(a)
    out = np.empty(B)
    if mode == "global":
        for b in range(B):
            out[b] = age_stat(w, a[rng.permutation(nb)])
        return out
    groups = collections.defaultdict(list)
    for j, g in enumerate(comp_idx):
        groups[g].append(j)
    groups = [np.array(v) for v in groups.values()]
    for b in range(B):
        ap = np.empty(nb)
        for idx in groups:
            ap[idx] = a[rng.permutation(idx)]
        out[b] = age_stat(w, ap)
    return out


def report(w, a, comp_idx, tag, rng, B, nulls):
    A_obs = age_stat(w, a)
    if not np.isfinite(A_obs):
        say("    %-8s no change mass on these branches, undefined" % tag)
        return None
    row = dict(readout=tag, A_obs=A_obs, n_branch=int(len(a)),
               n_unique_age=int(len(np.unique(np.round(a, 9)))))
    for mode, mtag in nulls:
        null = perm_null(w, a, comp_idx, mode, B, rng)
        mu, sd = float(null.mean()), float(null.std(ddof=1))
        p = (1 + int((np.abs(null - mu) >= abs(A_obs - mu)).sum())) / (1 + B)
        z = (A_obs - mu) / (sd + 1e-12) if sd > 0 else np.nan
        row.update({mtag + "_null_mean": mu, mtag + "_null_sd": sd,
                    mtag + "_p": p, mtag + "_z": z})
        note = ""
        if mode == "within_comp":
            note = "  [null INVALID as a test: ages are monotone along a chain, " \
                   "see header; diagnostic only]"
        say("    %-8s A=%+.5f | %s null %+.5f+-%.5f  p=%.4f  z=%s%s"
            % (tag, A_obs, mtag, mu, sd, p,
               "%+.2f" % z if np.isfinite(z) else "n/a (degenerate null)", note))
    cut = float(np.median(a))
    early = float(w[a < cut].sum() / w.sum()) if w.sum() > 0 else np.nan
    row.update(median_age=cut, mass_early=early)
    say("    %-8s mass on early half (a < %.3f) = %.3f ; branches %d"
        % (tag, cut, early, len(a)))
    return row


def run_clade(cl, B):
    lp = os.path.join(CAFE, cl, "lambda_deep.nwk")
    tp = os.path.join(CAFE, cl, "tree.nwk")
    cp = os.path.join(CAFE, cl, "deep", "Base_change.tab")
    for p in (lp, tp, cp):
        if not os.path.exists(p):
            say("%s: missing %s, skipped" % (cl, p))
            return None
    t_lab, im_lab = build_map(lp)
    t_tim, im_tim = build_map(tp)

    names_lab = [nd.name for nd in t_lab.get_leaves()]
    names_tim = [nd.name for nd in t_tim.get_leaves()]
    if names_lab != names_tim:
        say("%s: TIP ORDER MISMATCH, skipped" % cl)
        return None
    if len(im_lab) != len(im_tim):
        say("%s: node count mismatch %d vs %d, skipped"
            % (cl, len(im_lab), len(im_tim)))
        return None

    # the two trees hold different node objects; join them on the shared CAFE id
    tim_by_id = {im_tim[nd]: nd for nd in im_tim}
    lab2tim = {im_lab[nd]: tim_by_id[im_lab[nd]] for nd in im_lab}
    height = max(nd.get_distance(t_tim) for nd in t_tim.get_leaves())
    rtd = [nd.get_distance(t_tim) for nd in t_tim.get_leaves()]
    if max(rtd) - min(rtd) > 0.02 * height:
        say("%s: time tree not ultrametric (root-to-tip %.3f-%.3f), skipped"
            % (cl, min(rtd), max(rtd)))
        return None

    comps = deep_components(im_lab)
    K = len(comps)
    if K < 2:
        say("%s: K<2, skipped" % cl)
        return None

    df = pd.read_csv(cp, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    exp = set("%s<%d>" % (nd.name, i) if nd.is_leaf() else "<%d>" % i
              for nd, i in im_lab.items())
    if set(df.columns) != exp:
        say("%s: COLUMN MISMATCH, skipped" % cl)
        return None

    R, fam = recurrence(df, comps)

    # proximal age: the depth of the parent, i.e. when the branch began
    def prox_age(nid):
        nd = lab2tim[nid]
        pu = nd.up
        return 0.0 if pu is None else pu.get_distance(t_tim) / height

    deep_nd = [nd for nd in im_lab if int(nd.dist) == 2]
    shal_nd = [nd for nd in im_lab if int(nd.dist) == 1]
    deep_ids = sorted({im_lab[nd] for nd in deep_nd})
    shal_ids = sorted({im_lab[nd] for nd in shal_nd})
    nid2comp = {}
    for j, cs in enumerate(comps):
        for n in cs:
            nid2comp[n] = j
    n_multi = sum(1 for cs in comps if len(cs) > 1)

    a_deep = np.array([prox_age(n) for n in deep_ids])
    a_shal = np.array([prox_age(n) for n in shal_ids])
    comp_idx = np.array([nid2comp[n] for n in deep_ids])

    col2id = {c: col_node(c) for c in df.columns}
    cols_by_id = collections.defaultdict(list)
    for c, nid in col2id.items():
        cols_by_id[nid].append(c)

    M = df.abs()
    rec = fam["r"].to_numpy() >= 2

    def wvec(ids, mask=None):
        out = []
        for n in ids:
            cs = cols_by_id.get(n, [])
            if not cs:
                out.append(0.0)
            elif mask is None:
                out.append(float(M[cs].to_numpy(float).sum()))
            else:
                out.append(float(M.loc[mask, cs].to_numpy(float).sum()))
        return np.array(out)

    wdeep_all = wvec(deep_ids)
    wdeep_rec = wvec(deep_ids, rec)
    wshal_all = wvec(shal_ids)
    wshal_rec = wvec(shal_ids, rec)

    say("\n%s  K=%d origins (%d non-singleton)  height=%.1f Myr  R=%d  "
        "recurring families=%d/%d"
        % (cl, K, n_multi, height, R, int(rec.sum()), len(df)))
    say("  DEEP    branches=%d  age min %.4f med %.4f max %.4f  unique=%d"
        % (len(a_deep), a_deep.min(), np.median(a_deep), a_deep.max(),
           len(np.unique(np.round(a_deep, 9)))))
    say("  SHALLOW branches=%d  age min %.4f med %.4f max %.4f  unique=%d"
        % (len(a_shal), a_shal.min(), np.median(a_shal), a_shal.max(),
           len(np.unique(np.round(a_shal, 9)))))
    say("  mass  deep all %.0f  deep recur %.0f  shallow all %.0f  shallow recur %.0f"
        % (wdeep_all.sum(), wdeep_rec.sum(), wshal_all.sum(), wshal_rec.sum()))

    # ---- terminal-vs-internal diagnostic (added after the v2 run)
    # On an ultrametric tree a DEEP branch that leads to a tip can start late for a
    # purely topological reason: the tip's parent is a young node.  The statistic
    # must not be carried by that class alone, so both the mass split and the
    # restricted arm are reported.
    leaf_ids = {im_lab[nd] for nd in t_lab.get_leaves()}
    dleaf = np.array([n in leaf_ids for n in deep_ids])
    sleaf = np.array([n in leaf_ids for n in shal_ids])
    say("  topology  deep: %d/%d lead to a tip, %d/%d to an internal node ; "
        "mass all %.0f vs %.0f"
        % (int(dleaf.sum()), len(dleaf), int((~dleaf).sum()), len(dleaf),
           wdeep_all[dleaf].sum(), wdeep_all[~dleaf].sum()))
    say("            shallow: %d/%d lead to a tip, %d/%d to an internal node ; "
        "mass all %.0f vs %.0f"
        % (int(sleaf.sum()), len(sleaf), int((~sleaf).sum()), len(sleaf),
           wshal_all[sleaf].sum(), wshal_all[~sleaf].sum()))
    if dleaf.any() and (~dleaf).any():
        say("            mean age deep->tip %.4f  deep->internal %.4f"
            % (a_deep[dleaf].mean(), a_deep[~dleaf].mean()))

    rng = np.random.default_rng(SEED)
    rows = []
    say("  -- deep arm (target of the question) --")
    for w, tag in ((wdeep_all, "ALL"), (wdeep_rec, "RECUR")):
        if w.sum() <= 0:
            say("    %-8s no change mass on deep branches" % tag)
            continue
        r = report(w, a_deep, comp_idx, tag, rng, B,
                   (("global", "N1"), ("within_comp", "N2")))
        if r:
            r.update(arm="deep", clade=cl, K=K, n_component_multi=n_multi, R=R)
            rows.append(r)
    if dleaf.any() and (~dleaf).any() and (~dleaf).sum() >= 5:
        say("  -- deep arm restricted to branches leading to INTERNAL nodes "
            "(topology-free; N1 only) --")
        for w, tag in ((wdeep_all, "ALL"), (wdeep_rec, "RECUR")):
            ww = w[~dleaf]
            if ww.sum() <= 0:
                continue
            r = report(ww, a_deep[~dleaf], comp_idx[~dleaf],
                       tag, rng, B, (("global", "N1"),))
            if r:
                r.update(arm="deep_int", clade=cl, K=K,
                         n_component_multi=n_multi, R=R)
                rows.append(r)
    say("  -- shallow arm (control; N1 only) --")
    for w, tag in ((wshal_all, "ALL"), (wshal_rec, "RECUR")):
        if w.sum() <= 0:
            say("    %-8s no change mass on shallow branches" % tag)
            continue
        r = report(w, a_shal, np.zeros(len(a_shal), dtype=int), tag, rng, B,
                   (("global", "N1"),))
        if r:
            r.update(arm="shallow", clade=cl, K=K, n_component_multi=n_multi, R=R)
            rows.append(r)
    return rows


def _phi(x):
    from math import erf, sqrt
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def main():
    say("G1 exploratory E4 v3 -- turnover timing across branches (proximal ages)")
    say("PREREG_G1 md5 8a19c101207c7016ee1affacb1ee8943 | s8 E4, no edit made")
    say("EXPLORATORY: not corrected for multiple testing; hypothesis-generating; "
        "cannot be a headline (PREREG_G1 s8)")
    say("B=%d seed=%d  deep nulls: N1 global, N2 within-component; shallow: N1"
        % (B, SEED))
    all_rows = []
    for cl in CLADES:
        r = run_clade(cl, B)
        if r:
            all_rows.extend(r)
    if not all_rows:
        say("\nnothing written")
        return
    df = pd.DataFrame(all_rows)
    cols = ["clade", "arm", "readout", "K", "n_component_multi", "R",
            "n_branch", "n_unique_age", "A_obs", "median_age", "mass_early",
            "N1_null_mean", "N1_null_sd", "N1_p", "N1_z",
            "N2_null_mean", "N2_null_sd", "N2_p", "N2_z"]
    df = df[[c for c in cols if c in df.columns]]
    df.to_csv(os.path.join(OUT, "E4_timing.tsv"), sep="\t", index=False)

    say("\n" + "=" * 78)
    say("CROSS-CLADE (exploratory) -- equal-weight Stouffer on signed z")
    say("positive z = change mass on LATER branches")
    cross = []
    for arm in ("deep", "deep_int", "shallow"):
        sub_arm = df[df["arm"] == arm]
        if sub_arm.empty:
            continue
        for readout in ("ALL", "RECUR"):
            sub = sub_arm[sub_arm["readout"] == readout]
            if sub.empty:
                continue
            for mtag in ("N1", "N2"):
                if mtag + "_z" not in sub.columns:
                    continue
                z = sub[mtag + "_z"].to_numpy(float)
                z = z[np.isfinite(z)]
                if len(z) == 0:
                    say("  %-7s %-5s %s : no finite z" % (arm, readout, mtag))
                    continue
                Z = z.sum() / np.sqrt(len(z))
                p = 2 * (1 - _phi(abs(Z)))
                say("  %-7s %-5s %s : n=%d  z=%s  Stouffer Z=%+.3f  p=%.4f"
                    % (arm, readout, mtag, len(z),
                       " ".join("%+.2f" % x for x in z), Z, p))
                cross.append(dict(arm=arm, readout=readout, null=mtag,
                                  n_clade=len(z), Z=float(Z), p_two_sided=float(p)))
            nlate = int((sub["A_obs"] > 0).sum())
            say("  %-7s %-5s sign split: %d/%d clades late-weighted"
                % (arm, readout, nlate, len(sub)))
    if cross:
        pd.DataFrame(cross).to_csv(os.path.join(OUT, "E4_timing_perclade.tsv"),
                                   sep="\t", index=False)
    say("\nwrote E4_timing.tsv" + (" + E4_timing_perclade.tsv" if cross else ""))


if __name__ == "__main__":
    main()
