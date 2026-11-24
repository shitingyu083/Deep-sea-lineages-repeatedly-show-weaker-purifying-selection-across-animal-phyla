#!/usr/bin/env python
"""R3 (PREREG_nature.md), Chikina-style genome-scale fish RER: SRC dir for 36_rer_trees.sh from the Actinopteri codon
alignments (all real HOGs). Sequences restricted to Actinopteri species of meta/labels_v2.tsv; codons translated
(standard code; gapped/ambiguous codon -> '-', stop -> 'X'); HOGs with >= 20 fish species kept.
Writes work/addon/R3/fishrer/src/{final_aln/<hog>.faa, locus_stats.tsv}, topo_fish.tre (LSD2 time tree pruned)."""
import os, csv, glob, importlib.util
from ete3 import Tree
B = os.path.expanduser("~/deepsea_conv"); O = f"{B}/work/addon/R3/fishrer"; S = f"{O}/src"
os.makedirs(f"{S}/final_aln", exist_ok=True)
sp = importlib.util.spec_from_file_location("ps", f"{B}/scripts/63a_powersim_sty.py"); ps = importlib.util.module_from_spec(sp)
sp.loader.exec_module(ps)
fish = {r["code"] for r in csv.DictReader(open(f"{B}/meta/labels_v2.tsv"), delimiter="\t") if r["clade"] == "Actinopteri"}
t = Tree(open(f"{B}/work/final_v3/date/lsd.timetree.nwk").read().strip(), format=1)
keep = [x for x in t.get_leaf_names() if x in fish]
t.prune(keep, preserve_branch_length=True); t.write(outfile=f"{O}/topo_fish.tre", format=5)
rows = []; nsk = 0
for f in sorted(glob.glob(f"{B}/work/codon/Actinopteri/aln/*.fas")):
    h = os.path.basename(f)[:-4]
    if ".null" in h: continue
    a = {k: v for k, v in ps.read_fas(f).items() if k in keep}
    out = {}
    for k, s in a.items():
        aa = "".join(("X" if ps.CODE[s[i:i+3]] == "*" else ps.CODE[s[i:i+3]]) if s[i:i+3] in ps.CODE else "-"
                     for i in range(0, len(s) - 2, 3))
        if sum(c not in "-X" for c in aa) >= 30: out[k] = aa
    if len(out) < 20: nsk += 1; continue
    with open(f"{S}/final_aln/{h}.faa", "w") as fo:
        for k, v in out.items(): fo.write(f">{k}\n{v}\n")
    rows.append((h, len(out)))
with open(f"{S}/locus_stats.tsv", "w") as fo:
    fo.write("locus\tn_taxa\tmodel\n"); fo.writelines(f"{h}\t{n}\tLG+G4\n" for h, n in rows)
print(f"fish species in tree {len(keep)}; HOGs kept {len(rows)}; skipped (<20 fish) {nsk}")
