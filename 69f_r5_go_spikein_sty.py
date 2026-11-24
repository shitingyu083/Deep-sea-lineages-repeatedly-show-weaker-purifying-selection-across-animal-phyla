#!/usr/bin/env python3
"""R5 D2b (PREREG_r5.md): the 66d GO spike-in (code copied unchanged: same terms, test, clade z, Stouffer, BH within
namespace) extended to scenarios K50, K75 (R2 power) and K85, K90, K75h (R5 D2a power), fractions 1.0/0.5/0.25 of a
term's HOGs, and an effect present in only m of the term's k tested clades (m = 1..k, m clades drawn at random per term).
The 200 terms are drawn exactly as in 66d (default_rng(20261028)); the spike-in draws use default_rng(20261030).
Usage: python 69f_r5_go_spikein_sty.py FUNC_DIR"""
import csv, os, sys
from collections import defaultdict
import numpy as np
from scipy.stats import norm
B = os.path.expanduser("~/deepsea_conv"); FUNC = sys.argv[1]; O = f"{B}/work/addon/R5/go_spikein"; os.makedirs(O, exist_ok=True)
CL = "Actinopteri_m,Chondrichthyes,Mollusca,Annelida,Arthropoda,Cnidaria,Porifera,Echinodermata".split(",")
KMIN, KFRAC, MINCL, T = 5, 0.5, 3, "relax"
SCN, FRS = ("K50", "K75", "K85", "K90", "K75h"), (1.0, 0.5, 0.25)
def rd(fn): return list(csv.DictReader(open(fn), delimiter="\t"))
pw = defaultdict(dict)
for r in rd(f"{B}/work/addon/R5/D2/power_by_clade.tsv"): pw[r["scen"]][r["clade"]] = float(r["rate_qequiv"])
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
def qof(g, p):
    ns = names[g]; ps = np.array([q for h, q in P0.items() if names[h] == ns and h != g] + [p]); m = len(ps)
    o = np.sort(ps); r = np.searchsorted(o, p, side="right")
    qs = o * m / np.arange(1, m + 1); qs = np.minimum.accumulate(qs[::-1])[::-1]
    return float(min(1.0, qs[r - 1]))
terms = sorted(P0); pick = np.random.default_rng(20261028).choice(len(terms), size=min(200, len(terms)), replace=False)
rng = np.random.default_rng(20261030); rows = []
for j in pick:
    g = terms[j]; ktot = sum(len(v[2][g]) for v in CD.values() if g in v[2])
    tc = sorted(c for c, v in CD.items() if g in v[2] and zc(v[0] - v[1], v[2][g]) is not None); k = len(tc)
    for sc in SCN:
        for f in FRS:
            for m in range(1, k + 1):
                eff = set(rng.choice(tc, size=m, replace=False)); over = {}
                for c in eff:
                    Ir, In, mem = CD[c]; ix = mem[g]; ch = rng.choice(ix, size=max(1, int(round(f * len(ix)))), replace=False)
                    I2 = Ir.copy(); I2[ch] = np.maximum(I2[ch], rng.random(len(ch)) < pw[sc].get(c, 0.0)); over[c] = I2
                p = term_p(g, over)
                if p is None: continue
                rows.append((g, names[g], ktot, k, m, sc, f, P0[g], p, qof(g, p)))
with open(f"{O}/spikein_terms.tsv", "w") as fo:
    fo.write("term\tnamespace\tk_total\tk_clades\tm_effect\tscen\tfrac\tp_orig\tp_spike\tq_spike\n")
    for r in rows: fo.write("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in r) + "\n")
with open(f"{O}/spikein_power.tsv", "w") as fo:
    fo.write("scen\tfrac\tsubset\tn_terms\tpower_q05\n")
    for sc in SCN:
        for f in FRS:
            sub = [("all_clades", lambda r: r[4] == r[3])] + [(f"m={m}", lambda r, m=m: r[4] == m) for m in range(1, 9)] + \
                  [(f"m/k>={x}", lambda r, x=x: r[4] < r[3] and r[4] / r[3] >= x) for x in (0.5,)]
            for lab, fn in sub:
                x = [r for r in rows if r[5] == sc and r[6] == f and fn(r)]
                if x: fo.write(f"{sc}\t{f}\t{lab}\t{len(x)}\t{sum(r[9] < 0.05 for r in x) / len(x):.3f}\n")
# integrity: K50/K75 at frac 1.0/0.5 with all clades are the 66d design (different spike draws, so not bit-identical)
print(open(f"{O}/spikein_power.tsv").read()); open(f"{O}/SPIKE.DONE", "w").close()
