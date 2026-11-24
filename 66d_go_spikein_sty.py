#!/usr/bin/env python3
"""R2 (PREREG_nature.md L56-57): power of the D3 competitive GO test (62_func_competitive_sty.py, unchanged) under
spike-in. For 200 random GO terms tested in >= 3 clades (seed 20261028), a fraction f of the term's HOGs in each
clade is made truly affected: I_real := max(I_real, Bernoulli(pi_c)), pi_c = simulated RELAX power of clade c
(R2 power_by_clade.tsv, rate_qequiv) for scenario K50 or K75. The term is re-tested (clade mean/variance recomputed)
and its q obtained by BH with the other terms' original p values (same namespace). Test = 'relax'.
Usage: python 66d_go_spikein_sty.py FUNC_DIR"""
import csv, os, sys
from collections import defaultdict
import numpy as np
from scipy.stats import norm
B = os.path.expanduser("~/deepsea_conv"); FUNC = sys.argv[1]; O = f"{B}/work/addon/R2/go_spikein"; os.makedirs(O, exist_ok=True)
CL = "Actinopteri_m,Chondrichthyes,Mollusca,Annelida,Arthropoda,Cnidaria,Porifera,Echinodermata".split(",")
KMIN, KFRAC, MINCL, T = 5, 0.5, 3, "relax"
def rd(fn): return list(csv.DictReader(open(fn), delimiter="\t"))
pw = defaultdict(dict)
for r in rd(f"{B}/work/addon/R2/main/power_by_clade.tsv"): pw[r["scen"]][r["clade"]] = float(r["rate_qequiv"])
names = {}
for c in CL:
    fn = f"{FUNC}/{c}/enrich_{T}.tsv"
    if os.path.exists(fn):
        for r in rd(fn): names[r["term"]] = r["namespace"]
CD = {}
for c in CL:
    go = {r["hog"]: [g for g in r["go"].split(",") if g] for r in rd(f"{FUNC}/{c}/hog_annot.tsv")}
    sel = [r for r in rd(f"{FUNC}/{c}/hog_sel.tsv") if r["test"] == T and go.get(r["hog"])]
    if len(sel) < 20: continue
    Ir = np.array([float(r["I_real"]) for r in sel]); In = np.array([float(r["I_null_mean"]) for r in sel])
    mem = defaultdict(list)
    for i, r in enumerate(sel):
        for g in set(go[r["hog"]]): mem[g].append(i)
    CD[c] = (Ir, In, {g: np.array(ix) for g, ix in mem.items() if KMIN <= len(ix) <= KFRAC * len(sel)})
def zc(d, ix):
    N, k = len(d), len(ix); mu, var = d.mean(), d.var()
    if var == 0: return None
    return (d[ix].mean() - mu) / np.sqrt(var / k * (N - k) / (N - 1))
def term_p(g, Ir_over=None):
    z, w = [], []
    for c, (Ir, In, mem) in CD.items():
        if g not in mem: continue
        Iu = Ir_over[c] if Ir_over and c in Ir_over else Ir
        v = zc(Iu - In, mem[g])
        if v is None: continue
        z.append(v); w.append(np.sqrt(len(mem[g])))
    if len(z) < MINCL: return None
    z, w = np.array(z), np.array(w); return float(norm.sf((w * z).sum() / np.sqrt((w ** 2).sum())))
P0 = {g: p for g in {g for v in CD.values() for g in v[2]} if g in names for p in [term_p(g)] if p is not None}
print("terms (orig) ", len(P0), "; q<0.05 check below")
def qof(g, p):
    ns = names[g]; ps = np.array([q for h, q in P0.items() if names[h] == ns and h != g] + [p]); m = len(ps)
    o = np.sort(ps); r = np.searchsorted(o, p, side="right")
    qs = o * m / np.arange(1, m + 1); qs = np.minimum.accumulate(qs[::-1])[::-1]
    return float(min(1.0, qs[r - 1]))
rng = np.random.default_rng(20261028)
terms = sorted(P0); pick = rng.choice(len(terms), size=min(200, len(terms)), replace=False)
rows = []
for j in pick:
    g = terms[j]; ktot = sum(len(v[2][g]) for v in CD.values() if g in v[2])
    for sc in ("K50", "K75"):
        for f in (1.0, 0.5):
            over = {}
            for c, (Ir, In, mem) in CD.items():
                if g not in mem: continue
                ix = mem[g]; ch = rng.choice(ix, size=max(1, int(round(f * len(ix)))), replace=False)
                I2 = Ir.copy(); I2[ch] = np.maximum(I2[ch], rng.random(len(ch)) < pw[sc].get(c, 0.0)); over[c] = I2
            p = term_p(g, over)
            if p is None: continue
            rows.append((g, names[g], ktot, sc, f, P0[g], p, qof(g, p)))
with open(f"{O}/spikein_terms.tsv", "w") as fo:
    fo.write("term\tnamespace\tk_total\tscen\tfrac\tp_orig\tp_spike\tq_spike\n")
    for r in rows: fo.write("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in r) + "\n")
with open(f"{O}/spikein_power.tsv", "w") as fo:
    fo.write("scen\tfrac\tsize_bin\tn_terms\tpower_q05\n")
    for sc in ("K50", "K75"):
        for f in (1.0, 0.5):
            for lab, lo, hi in (("all", 0, 1e9), ("k<30", 0, 30), ("30-99", 30, 100), ("k>=100", 100, 1e9)):
                x = [r for r in rows if r[3] == sc and r[4] == f and lo <= r[2] < hi]
                if x: fo.write(f"{sc}\t{f}\t{lab}\t{len(x)}\t{sum(r[7] < 0.05 for r in x) / len(x):.3f}\n")
print(open(f"{O}/spikein_power.tsv").read())
open(f"{O}/SPIKE.DONE", "w").close()
