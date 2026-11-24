"""r10 (2026-11-07): update Source_Data.xlsx of submission_ready_20261107 in place for the v4 main figures.
Fig2/Fig3/Fig4 sheets <- viz_v4_sty/out fig<k>_source_v4.tsv + fig<k>_facts_v4.tsv; new sheet Fig2a_points with the
RELAX K and P of every fit plotted in Fig. 2a (R6/R6_3/r6_3_fits.tsv); old Fig3 sheet -> ED_Fig8; ED_Fig8 -> ED_Fig9;
ED_Fig9 -> ED_Fig10; Contents rewritten. Not idempotent: asserts the r9 sheet layout.
Usage (from paper/): python scripts/r10_source_data_sty.py submission_ready_20261107"""
import csv, os, sys
import openpyxl
from openpyxl.styles import Font

PKG = sys.argv[1]
X = os.path.join(PKG, "source_data", "Source_Data.xlsx")
V4 = os.path.join("..", "viz_v4_sty", "out")
wb = openpyxl.load_workbook(X)
assert "ED_Fig10" not in wb.sheetnames and wb["Fig3"]["A1"].value.startswith("Fig. 3 plotted values (fig3_source_v3")
TITLE = "Recurrent relaxation of purifying selection across independent deep-sea animal lineages"


def tsv(p):
    with open(p, encoding="utf-8") as f:
        return [r for r in csv.reader(f, delimiter="\t")]


# renames first (ED9 -> ED10, ED8 -> ED9, Fig3 -> ED8)
wb["ED_Fig9"].title = "ED_Fig10"
wb["ED_Fig8"].title = "ED_Fig9"
wb["Fig3"].title = "ED_Fig8"
ws = wb["ED_Fig8"]; ws["A1"] = ws["A1"].value.replace("Fig. 3 plotted values", "Extended Data Fig. 8 plotted values (former main Fig. 3 v3)")
for k in (2, 3, 4):
    name = "Fig%d" % k
    if name in wb.sheetnames:
        del wb[name]
    ws = wb.create_sheet(name, index=k)
    ws.append(["Fig. %d plotted values (fig%d_source_v4.tsv)" % (k, k)]); ws["A1"].font = Font(bold=True)
    for r in tsv(os.path.join(V4, "fig%d_source_v4.tsv" % k)):
        ws.append(r)
    ws.append([]); ws.append(["Fig. %d facts (fig%d_facts_v4.tsv)" % (k, k)]); ws.cell(ws.max_row, 1).font = Font(bold=True)
    for r in tsv(os.path.join(V4, "fig%d_facts_v4.tsv" % k)):
        ws.append(r)
# every plotted point of Fig. 2a
ws = wb.create_sheet("Fig2a_points", index=3)
ws.append(["Fig. 2a: RELAX K and P of every fit (R6/R6_3/r6_3_fits.tsv; real = deep foreground, null1 = matched null)"])
ws["A1"].font = Font(bold=True)
ws.append(["clade", "foreground", "orthogroup", "relax_K", "relax_P"])
n = 0
with open(os.path.join("data_srv", "work", "addon", "R6", "R6_3", "r6_3_fits.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        ws.append([r["clade"], r["label"], r["hog"], r["relax_K"], r["relax_p"]]); n += 1
assert n == 11200, n
# order: Contents, Fig1-4, Fig2a_points, ED_Fig1-10, rest
order = ["Contents", "Fig1", "Fig2", "Fig2a_points", "Fig3", "Fig4"] + ["ED_Fig%d" % i for i in range(1, 11)]
rest = [s for s in wb.sheetnames if s not in order]
wb._sheets = [wb[s] for s in order + rest]
C = wb["Contents"]
rows = [r for r in C.iter_rows(min_row=3, values_only=True)]
keep = {r[0]: r for r in rows}
C.delete_rows(1, C.max_row)
C.append(["Source Data for: " + TITLE]); C["A1"].font = Font(bold=True)
C.append(["Sheet", "Figure / panel", "Origin file (package build/)"])
new = {"Fig2": ("Fig. 2 (all panels; summaries)", "viz_v4_sty/out/fig2_source_v4.tsv, fig2_facts_v4.tsv"),
       "Fig2a_points": ("Fig. 2a (every plotted fit)", "addon/R6/R6_3/r6_3_fits.tsv"),
       "Fig3": ("Fig. 3 (all panels)", "viz_v4_sty/out/fig3_source_v4.tsv, fig3_facts_v4.tsv"),
       "Fig4": ("Fig. 4 (all panels)", "viz_v4_sty/out/fig4_source_v4.tsv, fig4_facts_v4.tsv"),
       "ED_Fig8": ("Extended Data Fig. 8 (former main Fig. 3)", "viz_v3_sty/out/fig3_source_v3.tsv, fig3_facts_v3.tsv"),
       "ED_Fig9": ("Extended Data Fig. 9", "ed_facts.tsv (former Extended Data Fig. 8)"),
       "ED_Fig10": ("Extended Data Fig. 10", "ed_facts.tsv (former Extended Data Fig. 9)")}
for s in wb.sheetnames[1:]:
    if s in new:
        C.append([s, *new[s]])
    else:
        C.append(list(keep[s]))
wb.save(X)
print("sheets", wb.sheetnames)
