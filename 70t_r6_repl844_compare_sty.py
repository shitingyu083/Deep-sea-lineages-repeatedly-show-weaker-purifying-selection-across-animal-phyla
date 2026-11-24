#!/usr/bin/env python3
# DEVIATIONS_r6 (e) evaluation: original R6-1 fits (r6_1_fits.tsv) vs one replicate (refit/repl844/).
# Rules exactly as (e): random stratum decides; failed replicate = disagreement; relaxed = P < 0.05 & K < 1.
# Writes refit/repl844_compare.tsv (one row per job) for 70t_r6_repl844_logor_sty.R (Delta log OR).
import csv, json, math, os
O = os.path.expanduser("~/deepsea_conv/work/addon/R6"); F = f"{O}/refit/repl844"
fits = {(r["clade"], r["hog"], r["set"], r["arm"]): r for r in csv.DictReader(open(f"{O}/r6_1_fits.tsv"), delimiter="\t")}
strata = list(csv.DictReader(open(f"{O}/repl844_strata.tsv"), delimiter="\t"))
cl = lambda K: math.log10(min(max(K, 1e-3), 50))
fl = lambda x: float(x) if x not in ("", "NA", "nan", None) else float("nan")
def call(p, K, kind):
    if p != p or K != K: return 0
    return int(p < 0.05 and (K < 1 if kind == "relax" else K > 1))
def wilson(x, n, z=1.959964):
    if n == 0: return (float("nan"),) * 2
    ph = x / n; d = 1 + z * z / n; c = ph + z * z / (2 * n); h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)
rows = []
for s in strata:
    k = (s["clade"], s["hog"], s["set"], s["arm"]); r = fits[k]
    po, Ko = fl(r["relax_p"]), fl(r["relax_K"]); ok_o = r["fit"] == "ok"
    j = f"{F}/{k[0]}/{k[1]}.{k[2]}.{k[3]}.relax.json"
    if os.path.exists(j):
        t = json.load(open(j))["test results"]; pr, Kr, st = float(t["p-value"]), float(t["relaxation or intensification parameter"]), "ok"
    else:
        pr = Kr = float("nan"); st = "missing"
    rows.append(dict(clade=k[0], hog=k[1], set=k[2], arm=k[3], stratum=s["stratum"], orig_ok=int(ok_o), repl=st,
                     p_o=po, K_o=Ko, p_r=pr, K_r=Kr,
                     rel_o=call(po, Ko, "relax"), rel_r=call(pr, Kr, "relax"),
                     int_o=call(po, Ko, "int"), int_r=call(pr, Kr, "int")))
with open(f"{O}/refit/repl844_compare.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
fails = set()
rc = f"{F}/logs/rc.txt"
if os.path.exists(rc):
    for ln in open(rc):
        a = ln.split()
        if len(a) >= 5 and a[4] != "0": fails.add(tuple(a[:4]))
done = lambda r: r["repl"] == "ok" or (r["clade"], r["hog"], r["set"], r["arm"]) in fails
for st in ("random", "boundary"):
    R = [r for r in rows if r["stratum"] == st]; D = [r for r in R if done(r)]
    agr = sum(r["repl"] == "ok" and r["rel_o"] == r["rel_r"] for r in D)
    ia = sum(r["repl"] == "ok" and r["int_o"] == r["int_r"] for r in D)
    nf = sum(r["repl"] != "ok" for r in D)
    lo, hi = wilson(agr, len(D))
    dk = sorted(abs(cl(r["K_o"]) - cl(r["K_r"])) for r in D if r["repl"] == "ok" and r["K_o"] == r["K_o"])
    med = dk[len(dk) // 2] if dk else float("nan"); p90 = dk[int(0.9 * (len(dk) - 1))] if dk else float("nan")
    print(f"{st}: planned {len(R)}  finished {len(D)}  replicate failures {nf}")
    print(f"  relaxed-call agreement {agr}/{len(D)} = {agr / max(len(D), 1):.3f}  Wilson 95% {lo:.3f}-{hi:.3f}")
    print(f"  intensified-call agreement {ia}/{len(D)} = {ia / max(len(D), 1):.3f}")
    print(f"  |d clamped log10 K| median {med:.3f}  90th pct {p90:.3f}  max {dk[-1] if dk else float('nan'):.3f}")
    fl_ = [r for r in D if r["repl"] == "ok"]
    print(f"  relaxed calls orig {sum(r['rel_o'] for r in fl_)}  repl {sum(r['rel_r'] for r in fl_)}; "
          f"orig->not {sum(r['rel_o'] and not r['rel_r'] for r in fl_)}, not->relaxed {sum(r['rel_r'] and not r['rel_o'] for r in fl_)}")
print("END present:", os.path.exists(f"{F}/END"))
