#!/usr/bin/env python
"""R1 depth-threshold sensitivity (PREREG_nature.md L42-44), step 1: per real A4 row, the FG tips of the row's
hyphy.nwk and their minimum / NA count of depth_median (meta/habitat_temp.tsv). Read-only on all inputs.
Writes work/addon/R1/depth_thresh/fg_depth_rows.tsv."""
import os, csv, math, importlib.util
B = os.path.expanduser("~/deepsea_conv"); O = f"{B}/work/addon/R1/depth_thresh"; os.makedirs(O, exist_ok=True)
sp = importlib.util.spec_from_file_location("ps", f"{B}/scripts/63a_powersim_sty.py"); ps = importlib.util.module_from_spec(sp)
sp.loader.exec_module(ps)
def fl(x):
    try: return float(x)
    except Exception: return float("nan")
dep = {r["code"]: fl(r.get("depth_median")) for r in csv.DictReader(open(f"{B}/meta/habitat_temp.tsv"), delimiter="\t")}
def fgtips(tr):
    fg = set()
    def walk(n, inf):
        f = inf or n.tag == "FG"
        if not n.ch and f: fg.add(n.name)
        for c in n.ch: walk(c, f)
    walk(tr, False); return fg
n = 0
with open(f"{O}/fg_depth_rows.tsv", "w") as fo:
    fo.write("clade\thog\tn_fg\tn_fg_na\tfg_min_depth\n")
    seen = set()
    for r in csv.DictReader(open(f"{B}/work/addon/A4/sel_features.tsv"), delimiter="\t"):
        if r["set"] != "real": continue
        k = (r["clade"], r["hog"])
        if k in seen: continue
        seen.add(k)
        p = f"{B}/work/codon/{r['clade']}/trees/{r['hog']}.hyphy.nwk"
        if not os.path.exists(p): continue
        fg = fgtips(ps.parse(open(p).read()))
        d = [dep.get(x, float("nan")) for x in fg]
        ok = [x for x in d if x == x]
        fo.write(f"{r['clade']}\t{r['hog']}\t{len(fg)}\t{len(d) - len(ok)}\t{min(ok) if ok else 'NA'}\n"); n += 1
print("rows", n)
