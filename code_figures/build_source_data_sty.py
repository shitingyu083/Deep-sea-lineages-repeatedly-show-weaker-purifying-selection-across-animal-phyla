"""Collect the values plotted in every main and Extended Data figure into one Source Data workbook.

Inputs are the fact files written by fig1-4_sty.py / ed_sty.py and the src_*.tsv tables they read (build/).
Nothing is recomputed here; rows are copied verbatim.
Usage: python scripts/build_source_data_sty.py submission_ready_YYYYMMDD
"""
import csv, os, sys
from openpyxl import Workbook

pkg = sys.argv[1]
B = "build"
out = os.path.join(pkg, "source_data")
os.makedirs(out, exist_ok=True)


def rows(fn):
    with open(os.path.join(B, fn), encoding="utf-8") as f:
        return [r for r in csv.reader(f, delimiter="\t")]


wb = Workbook()
toc = wb.active
toc.title = "Contents"
toc.append(["Source Data for: Recurrent relaxed purifying selection in deep-sea animal lineages"])
toc.append(["Sheet", "Figure / panel", "Origin file (package build/)"])

for i in (1, 2, 3, 4):
    fn = "fig%d_facts.tsv" % i
    ws = wb.create_sheet("Fig%d" % i)
    for r in rows(fn):
        ws.append(r)
    toc.append(["Fig%d" % i, "Fig. %d (all panels)" % i, fn])

ed = rows("ed_facts.tsv")
for k in range(1, 9):
    sel = [r for r in ed if r and r[0].startswith("ED%d" % k)]
    ws = wb.create_sheet("ED_Fig%d" % k)
    ws.append(["panel", "item", "value", "detail/source"])
    for r in sel:
        ws.append(r)
    toc.append(["ED_Fig%d" % k, "Extended Data Fig. %d" % k, "ed_facts.tsv"])

for fn in sorted(os.listdir(B)):
    if fn.startswith("src_") and fn.endswith(".tsv"):
        name = fn[4:-4][:31]
        ws = wb.create_sheet(name)
        for r in rows(fn):
            ws.append(r)
        toc.append([name, "underlying table read by the figure scripts", fn])

p = os.path.join(out, "Source_Data.xlsx")
wb.save(p)
print("wrote", p, "sheets", len(wb.sheetnames))
