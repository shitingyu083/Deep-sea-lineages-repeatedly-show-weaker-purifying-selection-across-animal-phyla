"""Rebuild source_data/Source_Data.xlsx for the R6 figure plan (20261101 package).
Old sheets are renamed to the new numbering (Fig2 -> SFig1, Fig3 -> ED_Fig5, Fig4 -> SFig2, ED_Fig5-8 -> ED_Fig6-9);
new main Fig1-4 sheets are copied verbatim from the v2 figure outputs (viz_v2_sty/out) and the tables they read.
Nothing is recomputed. Usage (from deepsea_tmp/stage_sty): python <this> paper/submission_ready_20261101"""
import csv, os, sys
from openpyxl import load_workbook

pkg = sys.argv[1]
V = os.path.join("viz_v2_sty", "out")
D = os.path.join("paper", "data_srv")
fn = os.path.join(pkg, "source_data", "Source_Data.xlsx")
wb = load_workbook(fn)
assert "SFig1" not in wb.sheetnames, "already updated"


def tsv(p):
    with open(p, encoding="utf-8") as f:
        return [r for r in csv.reader(f, delimiter="\t")]


REN = [("ED_Fig8", "ED_Fig9"), ("ED_Fig7", "ED_Fig8"), ("ED_Fig6", "ED_Fig7"), ("ED_Fig5", "ED_Fig6"),
       ("Fig3", "ED_Fig5"), ("Fig2", "SFig1"), ("Fig4", "SFig2")]
for a, b in REN:
    wb[a].title = b
del wb["Fig1"]


def add(name, blocks):
    ws = wb.create_sheet(name)
    for title, rows in blocks:
        ws.append([title])
        for r in rows:
            ws.append(r)
        ws.append([])
    return ws


lab = tsv(os.path.join(D, "meta", "labels_v2.tsv"))
ci = {c: i for i, c in enumerate(lab[0])}
keep = ["code", "species", "phylum", "class", "order", "family", "clade", "lab"]
org = tsv(os.path.join(D, "work", "addon", "A2", "origins.tsv"))
oi = org[0].index("set")
add("Fig1", [("Fig. 1 summary values (fig1_facts_v2.tsv)", tsv(os.path.join(V, "fig1_facts_v2.tsv"))),
             ("Fig. 1a species habitat labels (D deep, S shallow, U unassigned; meta/labels_v2.tsv)",
              [keep] + [[r[ci[c]] for c in keep] for r in lab[1:]]),
             ("Fig. 1b disjoint deep origins, ages in Ma (addon/A2/origins.tsv)",
              [org[0]] + [r for r in org[1:] if r[oi] == "disjoint"])])
for k in (2, 3, 4):
    add("Fig%d" % k, [("Fig. %d plotted values (fig%d_source_v2.tsv)" % (k, k), tsv(os.path.join(V, "fig%d_source_v2.tsv" % k))),
                      ("Fig. %d summary values (fig%d_facts_v2.tsv)" % (k, k), tsv(os.path.join(V, "fig%d_facts_v2.tsv" % k)))])

order = ["Contents", "Fig1", "Fig2", "Fig3", "Fig4"] + ["ED_Fig%d" % i for i in range(1, 10)] + ["SFig1", "SFig2"]
rest = [s for s in wb.sheetnames if s not in order]
wb._sheets = [wb[s] for s in order + rest]

c = wb["Contents"]
c.delete_rows(3, c.max_row)
rows = [("Fig1", "Fig. 1 (all panels)", "viz_v2_sty/out/fig1_facts_v2.tsv; meta/labels_v2.tsv; addon/A2/origins.tsv")]
rows += [("Fig%d" % k, "Fig. %d (all panels)" % k, "viz_v2_sty/out/fig%d_source_v2.tsv, fig%d_facts_v2.tsv" % (k, k)) for k in (2, 3, 4)]
rows += [("ED_Fig%d" % i, "Extended Data Fig. %d" % i, "ed_facts.tsv" if i != 5 else "fig3_facts.tsv (former main Fig. 3)") for i in range(1, 10)]
rows += [("SFig1", "Supplementary Fig. 1", "fig2_facts.tsv (former main Fig. 2)"),
         ("SFig2", "Supplementary Fig. 2", "fig4_facts.tsv (former main Fig. 4)")]
rows += [(s, "underlying table read by the figure scripts", "src_%s.tsv" % s) for s in rest]
for r in rows:
    c.append(list(r))
wb.save(fn)
print("saved", fn, wb.sheetnames)
