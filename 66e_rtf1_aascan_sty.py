#!/usr/bin/env python
"""R3 (PREREG_nature.md L66-68): RTF1 fish amino-acid scan. Fish RTF1 codon alignment (work/codon/Actinopteri/aln/
N0.HOG0011415.fas) translated; deep set = Actinopteri D species with depth_median > 3000 m present in the alignment.
Statistic = number of sites where >= 2 deep-set species share a residue absent from every shallow (S) fish.
Null = 1000 random shallow sets of the same size (seed 20261028), each scored against the remaining shallow fish.
Writes work/addon/R3/rtf1/aascan.tsv and aascan_sites.tsv."""
import os, csv, importlib.util, random
from collections import Counter
B = os.path.expanduser("~/deepsea_conv"); O = f"{B}/work/addon/R3/rtf1"
sp = importlib.util.spec_from_file_location("ps", f"{B}/scripts/63a_powersim_sty.py"); ps = importlib.util.module_from_spec(sp)
sp.loader.exec_module(ps)
HOG = "N0.HOG0011415"
def fl(x):
    try: return float(x)
    except Exception: return float("nan")
lab = {r["code"]: r["lab"] for r in csv.DictReader(open(f"{B}/meta/labels_v2.tsv"), delimiter="\t") if r["clade"] == "Actinopteri"}
dep = {r["code"]: fl(r.get("depth_median")) for r in csv.DictReader(open(f"{B}/meta/habitat_temp.tsv"), delimiter="\t")}
a = ps.read_fas(f"{B}/work/codon/Actinopteri/aln/{HOG}.fas")
aa = {k: "".join(ps.CODE.get(s[i:i + 3], "-") if s[i:i + 3] in ps.CODE else "-" for i in range(0, len(s) - 2, 3))
      for k, s in a.items() if k in lab}
deep = sorted(k for k in aa if lab[k] == "D" and dep.get(k, 0) == dep.get(k, 0) and dep.get(k, 0) > 3000)
shal = sorted(k for k in aa if lab[k] == "S")
L = len(next(iter(aa.values())))
print("fish in aln", len(aa), "deep>3000", deep, "shallow", len(shal), "sites", L)
def score(grp, ref):
    hits = []
    for i in range(L):
        refs = {aa[k][i] for k in ref} - {"-", "*", "X"}
        c = Counter(aa[k][i] for k in grp if aa[k][i] not in "-*X")
        sh = [x for x, n in c.items() if n >= 2 and x not in refs]
        if sh and refs: hits.append((i + 1, "".join(sorted(sh))))
    return hits
obs = score(deep, shal)
random.seed(20261028); nul = []
m = len(deep)
for _ in range(1000):
    g = random.sample(shal, m); rest = [k for k in shal if k not in g]
    nul.append(len(score(g, rest)))
p = (1 + sum(x >= len(obs) for x in nul)) / 1001
nul.sort()
with open(f"{O}/aascan.tsv", "w") as fo:
    fo.write("hog\tn_fish\tdeep_set\tn_shallow\tn_sites\tobs\tnull_mean\tnull_q95\tp_upper\n")
    fo.write(f"{HOG}\t{len(aa)}\t{','.join(deep)}\t{len(shal)}\t{L}\t{len(obs)}\t{sum(nul) / 1000:.3f}\t{nul[949]}\t{p:.4g}\n")
with open(f"{O}/aascan_sites.tsv", "w") as fo:
    fo.write("aa_site\tshared_deep_states\n"); fo.writelines(f"{i}\t{s}\n" for i, s in obs)
print(open(f"{O}/aascan.tsv").read())
