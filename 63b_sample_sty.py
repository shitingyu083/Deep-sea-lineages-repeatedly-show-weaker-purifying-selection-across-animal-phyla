#!/usr/bin/env python
"""R2 main run set-up (PREREG_nature.md): check the pilot2 gate, then sample 20 HOGs per clade from the
null-HyPhy list (set real, in_nullhy == 1; seed 20261028) and write calib/job lists. Exit code 2 = gate fail."""
import os, sys, csv, random
B = os.path.expanduser("~/deepsea_conv"); P2 = f"{B}/work/addon/R2/pilot2"; R = f"{B}/work/addon/R2/main"
PY = f"{B}/envs/dc/bin/python"; S = f"{B}/scripts/63a2_powersim_sty.py"
SCN = ["H0r", "K75", "K50", "K50h", "H0b", "P1", "P2"]
ok = True; msg = []
for r in csv.DictReader(open(f"{P2}/pilot_res.tsv"), delimiter="\t"):
    if r["scen"] == "K50":
        k = float(r["K"]); g = 0.3 <= k <= 0.75; ok &= g; msg.append(f"{r['hog']} K50 K={k:.3f} {'ok' if g else 'FAIL'}")
for l in open(f"{P2}/treelen.tsv"):
    f = l.split("\t"); q = float(f[4]) / float(f[2]); g = 0.8 <= q <= 1.2; ok &= g
    msg.append(f"{f[0]} simL/fitL={q:.3f} scale={f[6].strip()} {'ok' if g else 'FAIL'}")
os.makedirs(R, exist_ok=True)
open(f"{R}/gate.txt", "w").write("\n".join(msg) + f"\nGATE\t{'PASS' if ok else 'FAIL'}\n")
if not ok: sys.exit(2)
rows = [r for r in csv.DictReader(open(f"{B}/work/addon/A4/sel_features.tsv"), delimiter="\t")
        if r["set"] == "real" and r["in_nullhy"] == "1"]
pick = []
for C in sorted(set(r["clade"] for r in rows)):
    hs = sorted(set(r["hog"] for r in rows if r["clade"] == C))
    s = random.Random(20261028).sample(hs, min(20, len(hs))); pick += [(C, h) for h in sorted(s)]
with open(f"{R}/hogs.tsv", "w") as fo:
    fo.write("clade\thog\n"); [fo.write(f"{C}\t{h}\n") for C, h in pick]
with open(f"{R}/calib_jobs.txt", "w") as fo:
    for C, h in pick: fo.write(f"{PY} {S} calib {C} {h} {R}/{C}/{h}/_fit\n")
with open(f"{R}/jobs.txt", "w") as fo:
    for C, h in pick:
        for sc in SCN: fo.write(f"{PY} {S} {C} {h} {sc} 0 {R}/{C}/{h}/{sc}\n")
print(open(f"{R}/gate.txt").read(), len(pick), "HOGs")
