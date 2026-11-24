#!/usr/bin/env python3
"""R6-4 lineage covariates (ADDENDUM_R6_4.md part 2). Reads species files only (no lineage effects).
Canonical deep lineages from R6/R6_2/canon.tsv (D lineages). Writes R6/R6_4/lineage_covariates.tsv and
R6/R6_4/cov_report.txt (stem-age cross-check against work/addon/A6_lsd/origins.tsv for single-tip origins).
Usage: envs/dc/bin/python 70p_r6_lineage_cov_sty.py"""
import csv, math, os, sys
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv"); R6 = f"{R}/work/addon/R6"; OUT = f"{R6}/R6_4"
tsv = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
num = lambda x: float(x) if x not in ("", "NA", "nan", "None", None) else float("nan")
os.makedirs(OUT, exist_ok=True)

lab = {r["code"]: r for r in tsv(f"{R}/meta/labels_v2.tsv")}
LAB = {c: r["lab"] for c, r in lab.items()}
pm = {r["code"]: r for r in tsv(f"{R}/work/physchem/species_metrics.tsv")}
ch = {r["code"]: r for r in tsv(f"{R}/meta/chemo_habitat.tsv")}
ht = {r["code"]: r for r in tsv(f"{R}/meta/habitat_temp.tsv")}
gs = {}
for r in tsv(f"{R}/meta/candidates_all.tsv"):
    v = num(r["total_len"])
    if v == v and v > 0:
        gs.setdefault(r["species"], []).append(v)

t = Tree(f"{R}/work/final_v3/date/lsd.timetree.nwk", format=1)
H = {}
for n in t.traverse("postorder"):
    H[n] = 0.0 if n.is_leaf() else max(H[c] + c.dist for c in n.children)
root_h = H[t]
scale = 1000.0 if root_h < 5 else 1.0          # LSD2 output in Gyr if root < 5
leaves = {l.name: l for l in t.get_leaves()}


def mean(v):
    v = [x for x in v if x == x]
    return sum(v) / len(v) if v else float("nan")


def enc_ratio(c):
    if c not in pm:
        return float("nan")
    s, e = num(pm[c]["busco_GC3"]), num(pm[c]["busco_ENC"])
    if s != s or e != e:
        return float("nan")
    ex = 2 + s + 29 / (s ** 2 + (1 - s) ** 2)
    return (ex - e) / ex


rows, miss = [], []
for r in tsv(f"{R6}/R6_2/canon.tsv"):
    if r["canon"][4] != "D":
        continue
    sp = r["species"].split(",")
    inn = [s for s in sp if s in leaves]
    if inn:
        m = leaves[inn[0]] if len(inn) == 1 else t.get_common_ancestor([leaves[s] for s in inn])
        crown = H[m] * scale; stem = H[m.up] * scale if m.up else float("nan")
        nS = sum(1 for l in m.get_leaves() if LAB.get(l.name) == "S")
    else:
        crown = stem = float("nan"); nS = -1
    dep = [num(lab[s]["depth_median"]) if s in lab else float("nan") for s in sp]
    rows.append({"clade": r["clade"], "canon": r["canon"], "n_species": len(sp), "n_in_timetree": len(inn),
                 "depth": mean([math.log10(x) for x in dep if x == x and x > 0]),
                 "age": math.log10(stem) if stem == stem and stem > 0 else float("nan"),
                 "stem_Ma": stem, "crown_Ma": crown, "nS_in_mrca": nS,
                 "enc_ratio": mean([enc_ratio(s) for s in sp]),
                 "chemo": int(any(ch.get(s, {}).get("chemo") == "1" for s in sp)),
                 "chemo_known": sum(1 for s in sp if s in ch),
                 "log_gsize": mean([math.log10(mean(gs[pm[s]["species"]])) for s in sp
                                    if s in pm and pm[s]["species"] in gs]),
                 "T_hab": mean([num(ht[s]["T_hab"]) for s in sp if s in ht]),
                 "species": ",".join(sp)})
    miss += [s for s in sp if s not in lab or s not in pm]

cols = list(rows[0].keys())
with open(f"{OUT}/lineage_covariates.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for x in rows:
        f.write("\t".join(f"{x[k]:.6g}" if isinstance(x[k], float) else str(x[k]) for k in cols) + "\n")

# cross-check: single-D-tip origins in A6_lsd/origins.tsv vs stem age here
o = {x["D_tips"]: num(x["stem_age"]) for x in tsv(f"{R}/work/addon/A6_lsd/origins.tsv") if x["n_D"] == "1"}
cmp = [(x["species"], x["stem_Ma"], o[x["species"]]) for x in rows if x["n_species"] == 1 and x["species"] in o]
with open(f"{OUT}/cov_report.txt", "w") as f:
    f.write(f"D lineages {len(rows)}; root height {root_h:.6g} scale {scale}\n")
    f.write(f"species missing labels/physchem: {len(set(miss))} {sorted(set(miss))[:30]}\n")
    for k in ("depth", "age", "enc_ratio", "log_gsize", "T_hab"):
        f.write(f"NA {k}: {sum(1 for x in rows if x[k] != x[k])}\n")
    f.write(f"chemo=1: {sum(x['chemo'] for x in rows)}; lineages with S in MRCA: "
            f"{sum(1 for x in rows if x['nS_in_mrca'] > 0)}; not in timetree: {sum(1 for x in rows if x['n_in_timetree'] == 0)}\n")
    f.write(f"stem-age cross-check vs A6_lsd single-tip origins: n {len(cmp)}; max abs diff "
            f"{max((abs(a - b) for _, a, b in cmp), default=float('nan')):.4g} Ma\n")
print(open(f"{OUT}/cov_report.txt").read())
