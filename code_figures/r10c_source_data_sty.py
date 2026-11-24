"""r10c (2026-11-07): Source_Data.xlsx of submission_ready_20261107 for the v5 Figs 2-4 (viz_v5_sty/out).
Fig2 / Fig3 / Fig4 sheets <- fig<k>_source_v5.tsv + fig<k>_facts_v5.tsv; new sheet Fig3ab_points with the per-
orthogroup omega-class changes plotted in Fig. 3a,b (R6/R6_3/r6_3_fits.tsv); Fig2a_points (all 11,200 fits, already
present) also covers Fig. 2b. Contents rows updated. Asserts the r10 workbook (md5 f3bde26f), so it runs once.
Usage (from paper/): python scripts/r10c_source_data_sty.py submission_ready_20261107"""
import csv, os, sys, hashlib
import numpy as np
import openpyxl
from openpyxl.styles import Font

PKG = sys.argv[1]
X = os.path.join(PKG, "source_data", "Source_Data.xlsx")
assert hashlib.md5(open(X, "rb").read()).hexdigest().startswith("f3bde26f")
V5 = os.path.join("..", "viz_v5_sty", "out")
wb = openpyxl.load_workbook(X)


def tsv(p):
    with open(p, encoding="utf-8") as f:
        return [r for r in csv.reader(f, delimiter="\t")]


for k in (2, 3, 4):
    name = "Fig%d" % k
    idx = wb.sheetnames.index(name)
    del wb[name]
    ws = wb.create_sheet(name, index=idx)
    ws.append(["Fig. %d plotted values (fig%d_source_v5.tsv)" % (k, k)]); ws["A1"].font = Font(bold=True)
    for r in tsv(os.path.join(V5, "fig%d_source_v5.tsv" % k)):
        ws.append(r)
    ws.append([]); ws.append(["Fig. %d facts (fig%d_facts_v5.tsv)" % (k, k)]); ws.cell(ws.max_row, 1).font = Font(bold=True)
    for r in tsv(os.path.join(V5, "fig%d_facts_v5.tsv" % k)):
        ws.append(r)
ws = wb.create_sheet("Fig3ab_points", index=wb.sheetnames.index("Fig3") + 1)
ws.append(["Fig. 3a,b: per-orthogroup change in BUSTED-PH omega classes, foreground minus background "
           "(R6/R6_3/r6_3_fits.tsv; real = deep foreground, null1 = matched null)"])
ws["A1"].font = Font(bold=True)
ws.append(["clade", "orthogroup", "foreground", "Delta_pur (log10)", "Delta_pos", "w_pur_T", "w_pur_B", "p_pos_T",
           "p_pos_B"])
f = lambda x: float(x) if x not in ("", "NA", "nan") else float("nan")
rows = [r for r in csv.DictReader(open(os.path.join("data_srv", "work", "addon", "R6", "R6_3", "r6_3_fits.tsv"),
                                       encoding="utf-8"), delimiter="\t") if r["bph_ok"] == "1"]
have = {}
for r in rows:
    have.setdefault((r["clade"], r["hog"]), set()).add(r["label"])
n = 0
for r in rows:
    if have[(r["clade"], r["hog"])] != {"real", "null1"}:
        continue
    dp = np.log10(f(r["wpur_T"]) + .001) - np.log10(f(r["wpur_B"]) + .001)
    dq = f(r["ppos_T"]) - f(r["ppos_B"])
    ws.append([r["clade"], r["hog"], r["label"], "%.6f" % dp, "%.6f" % dq, r["wpur_T"], r["wpur_B"], r["ppos_T"],
               r["ppos_B"]]); n += 1
assert n == 6400, n                                   # 400 paired orthogroups x 8 clades x 2 foregrounds
C = wb["Contents"]
for row in C.iter_rows(min_row=3):
    s = row[0].value
    if s in ("Fig2", "Fig3", "Fig4"):
        k = s[-1]
        row[1].value = "Fig. %s (all panels; summaries)" % k
        row[2].value = "viz_v5_sty/out/fig%s_source_v5.tsv, fig%s_facts_v5.tsv" % (k, k)
    if s == "Fig2a_points":
        row[1].value = "Fig. 2a,b (every plotted fit; pairs give 2b)"
pos = [r[0].value for r in C.iter_rows(min_row=1)].index("Fig3") + 1
C.insert_rows(pos + 1)
C.cell(pos + 1, 1, "Fig3ab_points"); C.cell(pos + 1, 2, "Fig. 3a,b (every plotted orthogroup)")
C.cell(pos + 1, 3, "addon/R6/R6_3/r6_3_fits.tsv")
wb.save(X)
print("sheets", wb.sheetnames[:8], "| Fig3ab_points rows", n)
