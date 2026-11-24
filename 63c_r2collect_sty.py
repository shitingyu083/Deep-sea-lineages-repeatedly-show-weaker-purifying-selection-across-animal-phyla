#!/usr/bin/env python
"""R2 collect (PREREG_nature.md): per-clade FPR/power at nominal p<0.05 and at the q<0.05-equivalent nominal p of
the real data; informative flag (power >= 0.8); recurrence power P(>=3 clades | k clades covered)."""
import os, csv, itertools, collections
B = os.path.expanduser("~/deepsea_conv"); R = f"{B}/work/addon/R2/main"
res = [r for r in csv.DictReader(open(f"{R}/res_all.tsv"), delimiter="\t")]
sel = [r for r in csv.DictReader(open(f"{B}/work/addon/A4/sel_features.tsv"), delimiter="\t") if r["set"] == "real"]

def fl(x):
    try: return float(x)
    except Exception: return None

def qthr(C, pc, qc):  # largest nominal p with q<0.05; if none, BH first-rank threshold 0.05/n
    v = [(fl(r[pc]), fl(r[qc])) for r in sel if r["clade"] == C and fl(r[pc]) is not None]
    s = [p for p, q in v if q is not None and q < 0.05]
    return (max(s) if s else 0.05 / max(1, len(v))), len(v)

out = []; power = collections.defaultdict(dict)
for C in sorted(set(r["clade"] for r in res)):
    tr, nr = qthr(C, "relax_p", "relax_q"); tb, nb = qthr(C, "bph_p", "bph_q")
    for sc in ["H0r", "K75", "K50", "K50h", "H0b", "P1", "P2"]:
        rr = [r for r in res if r["clade"] == C and r["scen"] == sc and fl(r["p"]) is not None]
        if sc in ("H0r", "K75", "K50", "K50h"):
            th = tr; det = [fl(r["p"]) for r in rr if fl(r["K"]) is not None and fl(r["K"]) < 1]  # relaxation direction
        else:
            th = tb; det = [fl(r["p"]) for r in rr]
        n = len(rr)
        a05 = sum(p < 0.05 for p in det) / n if n else float("nan"); aq = sum(p < th for p in det) / n if n else float("nan")
        Kmed = sorted(fl(r["K"]) for r in rr if fl(r["K"]) is not None)
        out.append([C, sc, n, f"{a05:.3f}", f"{th:.3g}", f"{aq:.3f}", f"{Kmed[len(Kmed)//2]:.3f}" if Kmed else ""])
        power[sc][C] = aq
with open(f"{R}/power_by_clade.tsv", "w") as fo:
    fo.write("clade\tscen\tn\trate_p05\tq_equiv_p\trate_qequiv\tmedian_K\n")
    for r in out: fo.write("\t".join(map(str, r)) + "\n")

def pb_ge(ps, m):  # Poisson-binomial P(X >= m)
    d = [1.0]
    for p in ps:
        d = [a * (1 - p) + (d[i - 1] * p if i else 0) for i, a in enumerate(d)] + [d[-1] * p]
    return sum(d[m:])

with open(f"{R}/recurrence_power.tsv", "w") as fo:
    fo.write("scen\tk_clades\tmean_P_ge3\tmin\tmax\n")
    for sc in ["K75", "K50", "K50h", "P1", "P2"]:
        cl = sorted(power[sc])
        for k in range(3, len(cl) + 1):
            v = [pb_ge([power[sc][c] for c in s], 3) for s in itertools.combinations(cl, k)]
            fo.write(f"{sc}\t{k}\t{sum(v)/len(v):.3f}\t{min(v):.3f}\t{max(v):.3f}\n")
with open(f"{R}/informative.tsv", "w") as fo:
    fo.write("clade\tscen\tpower_qequiv\tinformative(>=0.8)\n")
    for C, sc, n, a05, th, aq, km in out:
        if sc not in ("H0r", "H0b"): fo.write(f"{C}\t{sc}\t{aq}\t{'yes' if float(aq) >= 0.8 else 'no'}\n")
print(open(f"{R}/power_by_clade.tsv").read()); print(open(f"{R}/recurrence_power.tsv").read())
