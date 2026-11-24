#!/usr/bin/env python3
"""R6-1 collect (PREREG_r6.md): one row per HOG x set (real/null1) x arm (O/U) with RELAX p, K, fit status and
alignment diagnostics (length, gap fraction, species), joined to build_report.tsv. Also writes the production
RELAX p/K of the same HOGs (set real/null1, arm P) for continuity. Writes OUT/r6_1_fits.tsv.
Usage: python 70e_r6_collect_sty.py"""
import csv, json, os

R = os.path.expanduser("~/deepsea_conv"); OUT = f"{R}/work/addon/R6"
NA = "NA"


def relax(p):
    try:
        with open(p) as f:
            tr = json.load(f)["test results"]
        return float(tr["p-value"]), float(tr["relaxation or intensification parameter"]), "ok"
    except FileNotFoundError:
        return NA, NA, "missing"
    except Exception:
        return NA, NA, "unreadable"


def aln(p):
    seqs, k = {}, None
    for l in open(p):
        if l.startswith(">"):
            k = l[1:].strip(); seqs[k] = []
        else:
            seqs[k].append(l.strip())
    s = ["".join(v) for v in seqs.values()]
    L = len(s[0]) if s else 0
    return L, (sum(x.count("-") for x in s) / (L * len(s)) if L else NA)


rep = list(csv.DictReader(open(f"{OUT}/build_report.tsv"), delimiter="\t"))
cols = ["clade", "hog", "set", "arm", "relax_p", "relax_K", "fit", "L", "gap_frac", "n", "nD", "nS", "n_prod",
        "drop_other", "drop_unmapped", "drop_invalid"]
with open(f"{OUT}/r6_1_fits.tsv.tmp", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rep:
        if r["status"] != "ok":
            continue
        c, h = r["clade"], r["hog"]
        base = [r[k] for k in ("n", "nD", "nS", "n_prod", "drop_other", "drop_unmapped", "drop_invalid")]
        for s in ("real", "null1"):
            for a in ("O", "U"):
                p, K, st = relax(f"{OUT}/sel/{c}/{h}.{s}.{a}.relax.json")
                L, g = aln(f"{OUT}/aln/{c}/{h}.{a}.fas")
                f.write("\t".join(map(str, [c, h, s, a, p, K, st, L, g] + base)) + "\n")
            p, K, st = relax(f"{R}/work/codon/{c}/sel/{h}{'' if s == 'real' else '.null1'}.relax.json")
            f.write("\t".join(map(str, [c, h, s, "P", p, K, st, NA, NA] + base)) + "\n")
os.replace(f"{OUT}/r6_1_fits.tsv.tmp", f"{OUT}/r6_1_fits.tsv")
print("done")
