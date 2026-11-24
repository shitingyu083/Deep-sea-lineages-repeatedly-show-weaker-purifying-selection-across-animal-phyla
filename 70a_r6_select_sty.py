#!/usr/bin/env python3
"""R6-1 step a (PREREG_r6.md): fixed HOG sample + list of annotation-route species to re-annotate.
Per main clade: 200 HOGs drawn (seed 20261031, sorted list, random.sample) from tiers/nullhy.txt, i.e. HOGs whose
production run has both a real-label and a null1 RELAX fit. For every sampled HOG the species of the production
alignment are classed by meta/species_codes.tsv source: miniprot (already uniform route), annot (ncbi_annot / local /
extra_local -> re-annotated by 59_route_reannot.py --one), other (tsa, ...: excluded from both arms).
Writes OUT/sample.tsv (clade, hog, n_aln, n_mp, n_annot, n_other) and OUT/reannot_codes.txt (annot species not yet
done in work/addon/A1_reannot). Production files are only read.
Usage: python r6a_select_sty.py [OUT=work/addon/R6]"""
import csv, os, random, sys
R = os.path.expanduser("~/deepsea_conv")
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{R}/work/addon/R6"
CLADES = ["Actinopteri_m", "Chondrichthyes", "Mollusca", "Annelida", "Arthropoda", "Cnidaria", "Porifera",
          "Echinodermata"]
ANNOT = ("ncbi_annot", "local", "extra_local")
N, SEED = 200, 20261031
os.makedirs(OUT, exist_ok=True)
src = {r["code"]: r["source"] for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}
a1 = f"{R}/work/addon/A1_reannot"
done_a1 = {c for c in os.listdir(a1) if os.path.exists(f"{a1}/{c}/REANNOT.DONE")}
need, rows = set(), []
for c in CLADES:
    W = f"{R}/work/codon/{c}"
    pool = sorted(l.strip() for l in open(f"{W}/tiers/nullhy.txt") if l.strip())
    pool = [h for h in pool if os.path.exists(f"{W}/aln/{h}.fas") and os.path.exists(f"{W}/sel/{h}.relax.json")
            and os.path.exists(f"{W}/sel/{h}.null1.relax.json")]
    rng = random.Random(f"{SEED}_{c}")
    for h in sorted(rng.sample(pool, min(N, len(pool)))):
        sp = [l[1:].split()[0] for l in open(f"{W}/aln/{h}.fas") if l.startswith(">")]
        k = {"mp": 0, "annot": 0, "other": 0}
        for s in sp:
            t = src.get(s, "")
            if t == "miniprot":
                k["mp"] += 1
            elif t in ANNOT:
                k["annot"] += 1; need.add(s)
            else:
                k["other"] += 1
        rows.append((c, h, len(sp), k["mp"], k["annot"], k["other"]))
with open(f"{OUT}/sample.tsv", "w") as f:
    f.write("clade\thog\tn_aln\tn_mp\tn_annot\tn_other\n")
    f.writelines("\t".join(map(str, r)) + "\n" for r in rows)
todo = sorted(need - done_a1)
with open(f"{OUT}/reannot_codes.txt", "w") as f:
    f.writelines(s + "\n" for s in todo)
print(f"sampled {len(rows)} HOGs; annot species {len(need)} (A1 done {len(need & done_a1)}, to do {len(todo)})")
for c in CLADES:
    rr = [r for r in rows if r[0] == c]
    print(c, len(rr), "mean n/mp/annot/other", *[round(sum(r[i] for r in rr) / max(1, len(rr)), 1) for i in (2, 3, 4, 5)])
