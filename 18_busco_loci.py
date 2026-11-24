#!/usr/bin/env python3
"""Animalia-wide species-tree loci from proteome-mode compleasm (metazoa_odb12).
Takes 'Single' BUSCO hits of every QC-pass species (meta/qc_table.tsv), keeps loci
present in >= MINOCC of species, writes work/tree/loci/<busco>.faa (header = code).
Usage: python 18_busco_loci.py [MINOCC=0.5]
Env: TREE_DIR (default work/tree), CODES=comma list to restrict to a clade's species."""
import csv, os, sys
from collections import defaultdict

R = os.path.expanduser("~/deepsea_conv")
MINOCC = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
TD = os.environ.get("TREE_DIR", f"{R}/work/tree")
OUT = f"{TD}/loci"; os.makedirs(OUT, exist_ok=True)
CODES = set(os.environ["CODES"].split(",")) if os.environ.get("CODES") else None


def read_fa(f):
    seqs, name = {}, None
    for l in open(f):
        if l.startswith(">"):
            name = l[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(l.strip())
    return {k: "".join(v) for k, v in seqs.items()}


qc = [r for r in csv.DictReader(open(f"{R}/meta/qc_table.tsv"), delimiter="\t") if r["status"] == "pass"
      and (CODES is None or r["code"] in CODES)]
loci = defaultdict(dict)
for r in qc:
    c = r["code"]
    hits = {}
    for l in open(f"{R}/work/protqc/{c}/full_table.tsv"):
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 3:  # Missing BUSCOs have no gene id
            continue
        b, st, sid = f[:3]
        if st == "Single":
            hits[b] = sid
    if not hits:
        continue
    fa = read_fa(f"{R}/data/proteomes/{c}.faa")
    for b, sid in hits.items():
        if sid in fa:
            loci[b][c] = fa[sid].rstrip("*")
n_sp = len(qc)
kept = {b: d for b, d in loci.items() if len(d) >= MINOCC * n_sp}
for f in os.listdir(OUT):
    os.remove(f"{OUT}/{f}")
for b, d in kept.items():
    with open(f"{OUT}/{b}.faa", "w") as o:
        for c, s in d.items():
            o.write(f">{c}\n{s}\n")
occ = defaultdict(int)
for d in kept.values():
    for c in d:
        occ[c] += 1
with open(f"{TD}/locus_occupancy.tsv", "w") as o:
    o.write("code\tn_loci\n")
    for r in qc:
        o.write(f"{r['code']}\t{occ.get(r['code'], 0)}\n")
print(f"species {n_sp}; loci total {len(loci)}; kept {len(kept)} (occ >= {MINOCC})")
print("lowest-occupancy species:", sorted(occ.items(), key=lambda x: x[1])[:8])
