#!/usr/bin/env python3
"""R3 (PREREG_nature.md L69-73): GO rank enrichment of the genome-scale fish RER (binary deep label).
Gene statistic s = sign(Rho) * -log10(P) from work/addon/R3/fishrer/lab/binary.tsv; per GO term (5-500 tested HOGs,
HOG->GO from eggNOG as in 47: FUNC/<clade>/hog_annot.tsv) two-sided Mann-Whitney term vs rest; BH within namespace.
Also reports the four Wang et al. 2019 terms. Read-only on inputs. Usage: python 66c_fishrer_go_sty.py FUNC_DIR"""
import csv, os, sys, math
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu
B = os.path.expanduser("~/deepsea_conv"); FUNC = sys.argv[1]; O = f"{B}/work/addon/R3/fishrer/go"; os.makedirs(O, exist_ok=True)
WANG = ["GO:0031214", "GO:0001503", "GO:0006643", "GO:0006457"]
def rd(fn): return list(csv.DictReader(open(fn), delimiter="\t"))
ann = None
for c in ("Actinopteri_m", "Actinopteri"):
    fn = f"{FUNC}/{c}/hog_annot.tsv"
    if os.path.exists(fn): ann = fn; cl = c; break
if ann is None: print("no hog_annot"); sys.exit(0)
go = {r["hog"]: set(g for g in r["go"].split(",") if g) for r in rd(ann)}
names = {}
for t in ("relax", "bph_strict", "hyphy_any", "cs_qual"):
    fn = f"{FUNC}/{cl}/enrich_{t}.tsv"
    if os.path.exists(fn):
        for r in rd(fn): names[r["term"]] = (r["namespace"], r["name"])
rer = [r for r in rd(f"{B}/work/addon/R3/fishrer/lab/binary.tsv") if r["P"] not in ("NA", "") and r["Rho"] not in ("NA", "")]
s = {r["gene"]: math.copysign(-math.log10(max(float(r["P"]), 1e-300)), float(r["Rho"])) for r in rer}
genes = [g for g in s if go.get(g)]
print(f"annot {ann}; RER genes {len(s)}; with GO {len(genes)}")
v = np.array([s[g] for g in genes]); mem = defaultdict(list)
for i, g in enumerate(genes):
    for t in go[g]: mem[t].append(i)
rows = []
for t, ix in mem.items():
    k = len(ix)
    if k < 5 or k > 500: continue
    m = np.zeros(len(v), bool); m[ix] = True
    u, p = mannwhitneyu(v[m], v[~m], alternative="two-sided")
    ns, nm = names.get(t, ("NA", "NA"))
    rows.append([t, ns, nm, k, float(np.median(v[m])), float(np.median(v[~m])), float(u / (k * (len(v) - k))), float(p), float("nan")])
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); r = p[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]; q = np.empty(n); q[o] = np.minimum(r, 1); return q
for ns in {r[1] for r in rows}:
    ix = [i for i, r in enumerate(rows) if r[1] == ns]
    for i, q in zip(ix, bh([rows[i][7] for i in ix])): rows[i][8] = q
rows.sort(key=lambda r: r[7])
with open(f"{O}/fishrer_go_rank.tsv", "w") as f:
    f.write("term\tnamespace\tname\tk\tmedian_s_term\tmedian_s_rest\tAUC\tp\tq\n")
    for r in rows: f.write("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in r) + "\n")
print(f"terms {len(rows)}; q<0.05 {sum(r[8] < 0.05 for r in rows)}; q<0.05 & AUC>0.5 {sum(r[8] < 0.05 and r[6] > 0.5 for r in rows)}")
for r in rows[:15]: print("  ", r[0], r[1], r[2][:50], "k", r[3], "AUC %.3f p %.2g q %.2g" % (r[6], r[7], r[8]))
print("Wang 2019 terms:")
for w in WANG:
    hit = [r for r in rows if r[0] == w]
    print("  ", w, ("k %d AUC %.3f p %.3g q %.3g" % (hit[0][3], hit[0][6], hit[0][7], hit[0][8])) if hit else f"not tested (k={len(mem.get(w, []))})")
open(f"{O}/GO.DONE", "w").close()
