#!/usr/bin/env python3
"""Summarise the backbone MCMCTree runs (2026-10-25; PREREG_addon A6, convergence as PREREG_fig4 C).
Node numbering (PAML): tips = alignment order, internal nodes in order of '(' in tree.tre (root = ns + 1), so
mcmc.txt column t_n<ns+1+i> = i-th internal node in preorder; checked against each run's FigTree.tre (posterior
mean heights, matched by clade). Each backbone node is defined in the full tree by two tips (first tip of its first and
of its second child) -> LSD2 age of that MRCA (Ma; rate-scaled height x root age / root height).
Convergence: per clock model, chains r1 vs r2: Pearson r of node means > 0.999, max |diff| / mean < 5 %, per-chain
ESS (initial positive sequence) >= 200 for every node; PASS = clock 2 passes (clock 3 reported).
Draws for 50c: 100 equally spaced samples of the pooled clock-2 chains (node ages, Ma).
Usage: python 50b_mcmctree_collect.py D [--suffix x]  (D = OUT/n<N>; suffix x = the extension chains)
Writes D/{node_ages,draws,convergence}<sfx>.tsv|txt and D/mcmc_summary<sfx>.txt."""
import os, re, sys
import numpy as np
from ete3 import Tree

D = os.path.abspath(sys.argv[1])
SFX = sys.argv[sys.argv.index("--suffix") + 1] if "--suffix" in sys.argv else ""
R = os.path.expanduser("~/deepsea_conv")
FULL = f"{R}/work/final_v3/date/lsd.timetree.nwk"
RUNS = [f"c2_r1{SFX}", f"c2_r2{SFX}", f"c3_r1{SFX}", f"c3_r2{SFX}", "prior"]
clean = lambda s: re.sub(r"\[[^\]]*\]", "", s.strip())


def ess(x):
    x = x - x.mean(); n = len(x)
    if n < 10 or x.var() == 0:
        return float("nan")
    f = np.fft.rfft(x, 2 * n); ac = np.fft.irfft(f * np.conj(f))[:n] / (x.var() * n)
    s = 0.0
    for k in range(1, n - 1, 2):
        p = ac[k] + ac[k + 1]
        if p < 0:
            break
        s += p
    return n / (1 + 2 * s)


def load(r, ns, nint):
    fn = f"{D}/{r}/mcmc.txt"
    if not os.path.exists(fn):
        return None
    h = open(fn).readline().split()
    X = np.loadtxt(fn, skiprows=1, ndmin=2)
    want = [f"t_n{ns + 1 + i}" for i in range(nint)]
    if any(w not in h for w in want):
        sys.exit(f"{r}: mcmc.txt lacks {[w for w in want if w not in h][:3]}")
    return X[:, [h.index(w) for w in want]] * 100                  # Ma


def figtree_heights(r):
    fn = f"{D}/{r}/FigTree.tre"
    if not os.path.exists(fn):
        return None
    txt = open(fn).read()
    m = re.search(r"tree\s+\S+\s*=\s*(?:\[[^\]]*\]\s*)?(\(.*?;)", txt, re.I | re.S)
    if not m:
        return None
    t = Tree(clean(m.group(1)), format=1, quoted_node_names=True)
    for l in t:
        l.name = re.sub(r"^\d+_", "", l.name) if l.name not in TIPS else l.name
    return t
L = open(f"{D}/tree.tre").read().split("\n", 1)[1]
bb = Tree(re.sub(r"'[^']*'", "", clean(L)).strip(), format=9)    # calibration labels removed; topology only
TIPS = set(bb.get_leaf_names()); ns = len(TIPS)
inner = [n for n in bb.traverse("preorder") if not n.is_leaf()]
defs = [(n.children[0].get_leaf_names()[0], n.children[1].get_leaf_names()[0], len(n)) for n in inner]
full = Tree(clean(open(FULL).read()), format=1, quoted_node_names=True)
lsdf = FULL.replace(".nwk", ".lsd")
tm = [float(m.group(1)) for line in open(lsdf) for m in [re.search(r"tMRCA\s+(-?[\d.]+)", line)] if m]
H = max(full.get_distance(l) for l in full); sc = abs(tm[-1]) / H
lsd = [(H - full.get_distance(full.get_common_ancestor(a, b))) * sc for a, b, _ in defs]
X = {r: load(r, ns, len(inner)) for r in RUNS}
summ = []
chk = {}
for r, x in X.items():                                          # numbering check against FigTree.tre
    if x is None:
        continue
    try:
        ft = figtree_heights(r)
        hh = max(ft.get_distance(l) for l in ft)
        fa = np.array([(hh - ft.get_distance(ft.get_common_ancestor(a, b))) * 100 for a, b, _ in defs])
    except Exception as e:                                       # reported; the numbering is then UNVERIFIED
        chk[r] = f"numbering UNVERIFIED (FigTree.tre: {type(e).__name__} {str(e)[:60]})"; continue
    d = np.abs(fa - x.mean(0)); chk[r] = f"max |FigTree - mean| {d.max():.3f} Ma"
    if d.max() > max(0.5, 0.01 * x.mean(0).max()):
        sys.exit(f"{r}: node numbering check FAILED ({chk[r]})")
with open(f"{D}/node_ages{SFX}.tsv", "w") as f:
    f.write("node\ttipA\ttipB\tn_tips\tlsd2_Ma\t" + "\t".join(f"{r}_mean\t{r}_lo\t{r}_hi\t{r}_ess" for r in RUNS) + "\n")
    for j, (a, b, nt) in enumerate(defs):
        v = []
        for r in RUNS:
            if X[r] is None:
                v += ["NA"] * 4; continue
            c = X[r][:, j]
            v += [f"{c.mean():.2f}", f"{np.percentile(c, 2.5):.2f}", f"{np.percentile(c, 97.5):.2f}", f"{ess(c):.0f}"]
        f.write(f"t_n{ns + 1 + j}\t{a}\t{b}\t{nt}\t{lsd[j]:.2f}\t" + "\t".join(v) + "\n")
crit, verdict = {}, "FAIL"
for ck in ("c2", "c3"):
    a, b = X[f"{ck}_r1{SFX}"], X[f"{ck}_r2{SFX}"]
    if a is None or b is None:
        crit[ck] = (False, "chain missing"); continue
    ma, mb = a.mean(0), b.mean(0)
    r_ = np.corrcoef(ma, mb)[0, 1]; rd = np.max(np.abs(ma - mb) / ((ma + mb) / 2))
    e = min(np.nanmin([ess(c) for c in a.T]), np.nanmin([ess(c) for c in b.T]))
    ok = r_ > 0.999 and rd < 0.05 and e >= 200
    crit[ck] = (ok, f"r {r_:.5f}, max rel diff {rd:.4f}, min per-chain ESS {e:.0f}, root {ma[0]:.1f} / {mb[0]:.1f} Ma")
if crit.get("c2", (False,))[0]:
    verdict = "PASS"
with open(f"{D}/convergence{SFX}.txt", "w") as f:
    f.write(f"{verdict} (PREREG_fig4 C on the clock-2 chains)\n")
    for ck, (ok, s) in crit.items():
        f.write(f"{ck}: {'pass' if ok else 'fail'}; {s}\n")
pool = [X[r] for r in (f"c2_r1{SFX}", f"c2_r2{SFX}") if X[r] is not None]
if pool:
    P = np.vstack(pool); ix = np.linspace(0, len(P) - 1, 100).round().astype(int)
    with open(f"{D}/draws{SFX}.tsv", "w") as f:
        f.write("draw\tnode\ttipA\ttipB\tage_Ma\n")
        for k, i in enumerate(ix, 1):
            for j, (a, b, _) in enumerate(defs):
                f.write(f"{k}\tt_n{ns + 1 + j}\t{a}\t{b}\t{P[i, j]:.4f}\n")
with open(f"{D}/mcmc_summary{SFX}.txt", "w") as f:
    for r in RUNS:
        x = X[r]
        f.write(f"{r}: {'missing' if x is None else f'{x.shape[0]} samples, min ESS {np.nanmin([ess(c) for c in x.T]):.0f}, root {x[:, 0].mean():.1f} Ma'}"
                f"{'; ' + chk[r] if r in chk else ''}\n")
    for ck in ("c2", "c3"):
        x = X[f"{ck}_r1{SFX}"]
        if x is not None:
            f.write(f"{ck}_r1 vs LSD2 (backbone nodes): Pearson r {np.corrcoef(x.mean(0), lsd)[0, 1]:.4f}; "
                    f"root {x[:, 0].mean():.1f} vs {lsd[0]:.1f} Ma\n")
    if X[f"c2_r1{SFX}"] is not None and X[f"c3_r1{SFX}"] is not None:
        f.write(f"clock 2 vs 3: Pearson r {np.corrcoef(X[f'c2_r1{SFX}'].mean(0), X[f'c3_r1{SFX}'].mean(0))[0, 1]:.4f}\n")
    if X["prior"] is not None and X[f"c2_r1{SFX}"] is not None:
        w0 = X["prior"]; w1 = X[f"c2_r1{SFX}"]
        sh = np.median((np.percentile(w1, 97.5, 0) - np.percentile(w1, 2.5, 0)) / (np.percentile(w0, 97.5, 0) - np.percentile(w0, 2.5, 0)))
        f.write(f"posterior / prior 95 % interval width, median over nodes {sh:.3f}\n")
print(open(f"{D}/convergence{SFX}.txt").read() + open(f"{D}/mcmc_summary{SFX}.txt").read(), end="")
