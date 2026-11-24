"""Tier lists for 24_selection_run.sh v3 (2026-10-20).
eligible: HOGs with >= 2 independent deep lineages (maximal all-deep subtrees; convergence needs >= 2)
          and >= 1 shallow-shallow null set written by 23_label_trees.py (needs at least as many
          maximal all-shallow lineages as deep ones), so every listed HOG has real and null labels
          (paired calibration of the genome-wide rates).
hylist:   N_HY eligible HOGs drawn at random (seed 20261020; all of them if fewer) -> HyPhy on the real
          labels and CSUBST on all null sets.
nullhy:   the first N_NULLHY HOGs of the shuffled hylist -> HyPhy on null set 1.
Writes work/codon/<clade>/tiers/{eligible,hylist,nullhy}.txt (sorted) and prints the counts.
Usage: python 43c_sel_lists.py CLADE [N_HY=1000] [N_NULLHY=400]
"""
import csv
import os
import random
import sys

R = os.path.expanduser("~/deepsea_conv")
clade = sys.argv[1]
NHY = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
NNH = int(sys.argv[3]) if len(sys.argv) > 3 else 400
base = f"{R}/work/codon/{clade}"
od = f"{base}/tiers"
os.makedirs(od, exist_ok=True)
rows = list(csv.DictReader(open(f"{base}/lineages.tsv"), delimiter="\t"))
fg2 = [r["hog"] for r in rows if int(r["n_fg_lineages"]) >= 2]


def has_null(h):
    p = f"{base}/trees/{h}.null1.fg.txt"
    return os.path.exists(p) and os.path.getsize(p) > 0


elig = sorted(h for h in fg2 if has_null(h))
shuf = elig[:]
random.Random(20261020).shuffle(shuf)
hy, nh = shuf[:NHY], shuf[:NHY][:NNH]
for name, lst in (("eligible", elig), ("hylist", hy), ("nullhy", nh)):
    with open(f"{od}/{name}.txt", "w") as f:
        f.write("".join(h + "\n" for h in sorted(lst)))
nfg = [int(r["n_fg_lineages"]) for r in rows]
print(f"{clade}: hogs {len(rows)}; n_fg 0/1/>=2: {nfg.count(0)}/{nfg.count(1)}/{len(fg2)}; "
      f"eligible (>=2 deep + null set) {len(elig)}; hylist {len(hy)}; nullhy {len(nh)}")
