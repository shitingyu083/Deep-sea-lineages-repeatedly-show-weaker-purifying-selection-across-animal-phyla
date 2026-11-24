"""Assemble supplementary/Supplementary_Information.md from the frozen extended manuscript
(build/main_extended_20261031.md): Supplementary Notes 1-6 (verbatim original Results
paragraphs), Supplementary Methods (verbatim original Methods), and legends for the
Supplementary Tables (from the xlsx Contents sheet) and Supplementary Data 1.
Usage: python build_si_text_sty.py <extended.md> <package_dir>"""
import os, re, sys
from openpyxl import load_workbook

ext, pkg = sys.argv[1], sys.argv[2]
L = open(ext, encoding="utf-8").read().split("\n")
line = lambda n: L[n - 1]          # 1-based
TITLE = "Recurrent relaxed purifying selection in deep-sea animal lineages"


def swap_ed(s):
    # Extended Data Figs 4 and 5 were renumbered to first-citation order in the main paper
    s = re.sub(r"Extended Data Fig\. 4\b", "@@ED5", s)
    s = re.sub(r"Extended Data Fig\. 5\b", "Extended Data Fig. 4", s)
    return s.replace("@@ED5", "Extended Data Fig. 5")


def note_text(s):
    s = re.sub(r"\(Methods\)", "(Supplementary Methods)", s)
    s = re.sub(r"\(Methods, ", "(Supplementary Methods, ", s)
    s = re.sub(r"; Methods\)", "; Supplementary Methods)", s)
    return swap_ed(s)


NOTES = [
    (1, "Depth habitat, phylogeny and repeated origins", [23, 25, 27]),
    (2, "Proteome composition and threonine enrichment", [31, 33, 35]),
    (3, "Selection regime on deep branches", [39, 41, 43]),
    (4, "Lineage-wide and annotation-route confounders", [45, 47]),
    (5, "Relative evolutionary rates and convergent substitution", [49, 51]),
    (6, "Shared genes, functions, power, published candidates and gene families", [55, 57, 59, 61, 63]),
]
assert line(73) == "## Methods" and line(134) == "## Data availability", (line(73), line(134))
for _, _, ns in NOTES:
    for n in ns:
        assert line(n) and not line(n).startswith("#"), n

out = ["# Supplementary Information", "",
       "**%s**" % TITLE, "",
       "Reference numbers in this document refer to the reference list of the main paper. "
       "Figure and Extended Data figure numbers refer to the main paper.", "",
       "## Contents",
       "- Supplementary Notes 1–6",
       "- Supplementary Methods",
       "- Legends for Supplementary Tables 1–10 (file Supplementary_Tables_1-10.xlsx)",
       "- Legend for Supplementary Data 1 (file Supplementary_Data_1_per_gene_selection.tsv.gz)", ""]
out += ["## Supplementary Notes", "",
        "The Supplementary Notes give the full results summarized in the main text, "
        "with all secondary statistics and sensitivity analyses.", ""]
for k, t, ns in NOTES:
    out += ["### Supplementary Note %d | %s" % (k, t), ""]
    for n in ns:
        out += [note_text(line(n)), ""]
out += ["## Supplementary Methods", ""]
for n in range(75, 133):
    s = line(n)
    out.append(swap_ed(s))
out += [""]

wb = load_workbook(os.path.join(pkg, "supplementary", "Supplementary_Tables_1-10.xlsx"), read_only=True)
rows = list(wb["Contents"].iter_rows(values_only=True))
out += ["## Supplementary Table legends", "",
        "Each table is a sheet of Supplementary_Tables_1-10.xlsx; the Contents sheet maps sheets to "
        "the analysis files they were built from.", ""]
n_leg = 0
for r in rows[3:]:
    if not r or not r[1]:
        continue
    cap = swap_ed(str(r[1]).strip())
    m = re.match(r"(Supplementary Table \w+) \| (.*)", cap)
    if m:
        out += ["**%s |** %s" % (m.group(1), m.group(2)), ""]
        n_leg += 1
out += ["## Supplementary Data legend", "",
        "**Supplementary Data 1 |** Per-gene RELAX, BUSTED-PH and CSUBST results for all 60,906 "
        "orthogroup-by-clade tests with real (deep) foreground labels in the eight test clades "
        "(tab-separated, gzip-compressed).", ""]
p = os.path.join(pkg, "supplementary", "Supplementary_Information.md")
open(p, "w", encoding="utf-8", newline="\n").write("\n".join(out))
print("written", p, "table legends", n_leg)
