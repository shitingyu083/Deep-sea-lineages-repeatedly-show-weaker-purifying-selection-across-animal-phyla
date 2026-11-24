#!/usr/bin/env python3
"""R6-3 extraction (PREREG_r6.md): per production fit (real = tiers/hylist.txt, null1 = tiers/nullhy.txt) the RELAX
p and K and the BUSTED-PH Unconstrained-model ω classes of Test (FG) and Background (BG), plus membership of the
pre-specified GO sets (func_v2 hog_annot.tsv; eggNOG GO already propagated). Read-only on production.
Per branch set: wpur = proportion-weighted mean ω of classes with ω < 1 (renormalised), ppos = proportion of classes
with ω > 1, wpos = proportion-weighted mean ω of those classes.
Writes OUT/r6_3_fits.tsv. Usage: python 70f_r6_omega_sty.py [OUT work/addon/R6/R6_3]"""
import json, os, sys

R = os.path.expanduser("~/deepsea_conv")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{R}/work/addon/R6/R6_3"
CL = "Actinopteri_m,Chondrichthyes,Mollusca,Annelida,Arthropoda,Cnidaria,Porifera,Echinodermata".split(",")
SETS = dict(LIGHT={"GO:0007601", "GO:0007602", "GO:0009583", "GO:0009881", "GO:0001750", "GO:0001654"},
            CIRCADIAN={"GO:0007623", "GO:0048511"}, PIGMENT={"GO:0043473"})
NA = "NA"


def load(p):
    try:
        with open(p) as f:
            return json.load(f)
    except Exception:
        return None


def classes(rd):
    c = [(float(v["omega"]), float(v["proportion"])) for v in rd.values()]
    pur = [(w, p) for w, p in c if w < 1]; pos = [(w, p) for w, p in c if w > 1]
    sp = sum(p for _, p in pur); pp = sum(p for _, p in pos)
    return (sum(w * p for w, p in pur) / sp if sp > 0 else NA, pp,
            sum(w * p for w, p in pos) / pp if pp > 0 else NA)


def fmt(x):
    return x if x == NA else f"{x:.6g}"


os.makedirs(OUT, exist_ok=True)
cols = ["clade", "hog", "label", "relax_p", "relax_K", "bph_ok", "wpur_T", "ppos_T", "wpos_T", "wpur_B", "ppos_B",
        "wpos_B", "has_go"] + list(SETS)
with open(f"{OUT}/r6_3_fits.tsv.tmp", "w") as out:
    out.write("\t".join(cols) + "\n")
    for c in CL:
        W = f"{R}/work/codon/{c}"
        go = {}
        for l in list(open(f"{R}/work/func_v2/{c}/hog_annot.tsv"))[1:]:
            p = l.rstrip("\n").split("\t")
            go[p[0]] = set(p[6].split(",")) if len(p) > 6 and p[6] else set()
        for lab, tier, tag in (("real", "hylist.txt", ""), ("null1", "nullhy.txt", ".null1")):
            for h in (l.strip() for l in open(f"{W}/tiers/{tier}")):
                if not h:
                    continue
                r = load(f"{W}/sel/{h}{tag}.relax.json") or {}
                tr = r.get("test results", {})
                rp, rk = tr.get("p-value"), tr.get("relaxation or intensification parameter")
                b = load(f"{W}/sel/{h}{tag}.bph.json") or {}
                rd = b.get("fits", {}).get("Unconstrained model", {}).get("Rate Distributions", {})
                try:
                    T, B = classes(rd["Test"]), classes(rd["Background"]); ok = 1
                except Exception:
                    T = B = (NA, NA, NA); ok = 0
                g = go.get(h, set())
                row = [c, h, lab, NA if rp is None else fmt(float(rp)), NA if rk is None else fmt(float(rk)), ok,
                       *map(fmt, T), *map(fmt, B), int(bool(g))] + [int(bool(g & s)) for s in SETS.values()]
                out.write("\t".join(map(str, row)) + "\n")
os.replace(f"{OUT}/r6_3_fits.tsv.tmp", f"{OUT}/r6_3_fits.tsv")
print("done")
