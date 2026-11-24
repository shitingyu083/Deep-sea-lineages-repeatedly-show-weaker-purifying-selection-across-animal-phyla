"""Shared helpers for the v5 figures (Figs 2-4, 2026-11-07). Everything of v4 (and through it v3/v2: palette,
statistics, QA, berg/cell/terr helpers, clade colours) is imported unchanged; outputs go to viz_v5_sty/out only.

New benchmark forms (panel map: paper/figure_benchmark_map.md, section "Fig. 2-4 v5"; local copies _ref_sty/web/):
  lam6d  - Nature 2022 doi:10.1038/s41586-022-04418-5, Extended Data Fig. 6d: faceted volcano grid, boxed strip
           labels on top and (rotated) on the right, significant points in one highlight colour, others light grey,
           y ticks on the first column of each row, x ticks on the bottom row, no panel frames.
  stil2b - Stiller et al. 2024 Nature doi:10.1038/s41586-024-07323-1, Fig. 2b: time tree whose terminal branches run
           on as light clade-coloured lines with the clade name under them, a light violin per clade without outline
           and a darker clade-mean bar, silhouettes right of the violins, leader-line labels.
  stil5a - same paper, Fig. 5a: grouped box plots filled with one colour per series, black outlines, median and
           whiskers, black outlier dots, a row of stars (or exact P) over each group, framed legend with the box key.
  spid2b - Tong et al. 2022 Nat. Commun. doi:10.1038/s41467-022-34446-8, Fig. 2b: framed panel, horizontal filled
           violins with outline, black points, a black median dot on a line, row labels at the left, P text.
  cell3  - Xu et al. 2025 Cell doi:10.1016/j.cell.2025.01.002, Fig. 3: tree with clade labels and red node
           diamonds, depth-coloured squares, names, columns of numbers in value-sized circles, bar columns, grey row
           bands, pill headers, depth legend.
"""
import os, sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "viz_v4_sty"))
from figcommon_v4_sty import *  # noqa: F401,F403
import figcommon_v2_sty as _v2
V5 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(V5, exist_ok=True)
_v2.V2 = V5
