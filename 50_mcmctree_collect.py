#!/usr/bin/env python3
"""Summarise the MCMCTree runs (2026-10-25): per node posterior mean and 95 % interval of each run (mcmc.txt t_n*
columns, x 100 Ma), chain agreement (c2_r1 vs c2_r2: max |mean difference| / mean, Pearson r), ESS (initial positive
sequence) of every node age, and prior vs posterior. (LSD2 comparison: separate step on FigTree.tre.)
Usage: python 50_mcmctree_collect.py OUT  -> OUT/node_ages.tsv, OUT/mcmc_summary.txt"""
import os, sys
import numpy as np

OUT = sys.argv[1]
RUNS = ("c2_r1", "c2_r2", "c3_r1", "prior")


def load(r):
    fn = f"{OUT}/{r}/mcmc.txt"
    if not os.path.exists(fn):
        return None
    h = open(fn).readline().split()
    X = np.loadtxt(fn, skiprows=1)
    cols = [i for i, c in enumerate(h) if c.startswith("t_n")]
    return [h[i] for i in cols], X[:, cols] * 100


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


D = {r: load(r) for r in RUNS}
names = next((d[0] for d in D.values() if d), [])
with open(f"{OUT}/node_ages.tsv", "w") as f:
    f.write("node\t" + "\t".join(f"{r}_mean\t{r}_lo\t{r}_hi\t{r}_ess" for r in RUNS) + "\n")
    for j, n in enumerate(names):
        v = []
        for r in RUNS:
            if D[r] is None:
                v += ["NA"] * 4; continue
            x = D[r][1][:, j]
            v += [f"{x.mean():.2f}", f"{np.percentile(x, 2.5):.2f}", f"{np.percentile(x, 97.5):.2f}", f"{ess(x):.0f}"]
        f.write(n + "\t" + "\t".join(v) + "\n")
with open(f"{OUT}/mcmc_summary.txt", "w") as f:
    for r in RUNS:
        f.write(f"{r}: {'missing' if D[r] is None else f'{D[r][1].shape[0]} samples, min ESS {np.nanmin([ess(c) for c in D[r][1].T]):.0f}'}\n")
    if D["c2_r1"] and D["c2_r2"]:
        a, b = D["c2_r1"][1].mean(0), D["c2_r2"][1].mean(0)
        f.write(f"c2 chains: Pearson r {np.corrcoef(a, b)[0, 1]:.5f}; max |diff|/mean {np.max(np.abs(a - b) / ((a + b) / 2)):.4f}; "
                f"root {a[0]:.1f} vs {b[0]:.1f} Ma\n")
    if D["c2_r1"] and D["c3_r1"]:
        a, c = D["c2_r1"][1].mean(0), D["c3_r1"][1].mean(0)
        f.write(f"clock 2 vs 3: Pearson r {np.corrcoef(a, c)[0, 1]:.4f}; root {a[0]:.1f} vs {c[0]:.1f} Ma\n")
print(open(f"{OUT}/mcmc_summary.txt").read())
