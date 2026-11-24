#!/usr/bin/env python3
"""Species-tree loci v2: metazoa_odb12 single-copy BUSCOs (proteome-mode compleasm) of every
QC-pass species plus the choanoflagellate outgroups (data/outgroups/outgroups.tsv, compleasm in
work/protqc_outgroup/<code>/). A locus is kept if it is single-copy in >= MINOCC of the INGROUP
species (outgroup presence not required; rooting uses whichever loci carry them).
Usage: python 18c_busco_loci_v2.py TREE_DIR [MINOCC=0.5]
Env:   STATUSES=pass (comma list of meta/qc_table.tsv statuses to include; final tree = pass)
       EXCLUDE=code[,code] species left out although their status is included (final tree: Euncla, the
       identity-conflicted assembly of meta/topo_fix.tsv); default none
Writes TREE_DIR/loci/<busco>.faa (header = code), TREE_DIR/species.tsv, TREE_DIR/loci_table.tsv,
TREE_DIR/locus_occupancy.tsv"""
import csv, os, sys
from collections import defaultdict

R = os.path.expanduser("~/deepsea_conv")
TD = sys.argv[1]
MINOCC = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
STATUSES = set(os.environ.get("STATUSES", "pass").split(","))
EXCLUDE = {c for c in os.environ.get("EXCLUDE", "").split(",") if c}
OUT = f"{TD}/loci"; os.makedirs(OUT, exist_ok=True)


def read_fa(f):
    seqs, name = {}, None
    for l in open(f):
        if l.startswith(">"):
            name = l[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append("".join(l.split()))  # some lab proteomes (Hha, Rve) have spaces inside lines
    return {k: "".join(v) for k, v in seqs.items()}


AA = set("ACDEFGHIKLMNPQRSTVWYBZJX")


def clean(s):
    """Terminal stop dropped; internal stop ('*' or '.'), U/O and anything non-standard -> X."""
    s = s.upper().rstrip("*.")
    return "".join(c if c in AA else "X" for c in s)


def singles(table):
    hits = {}
    for l in open(table):
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) >= 3 and f[1] == "Single":  # Missing BUSCOs have no gene id
            hits[f[0]] = f[2]
    return hits


qc = [r for r in csv.DictReader(open(f"{R}/meta/qc_table.tsv"), delimiter="\t") if r["status"] in STATUSES]
if EXCLUDE:
    absent = EXCLUDE - {r["code"] for r in qc}
    qc = [r for r in qc if r["code"] not in EXCLUDE]
    print(f"excluded {sorted(EXCLUDE - absent)}" + (f"; not in the selected statuses: {sorted(absent)}" if absent else ""))
og = list(csv.DictReader(open(f"{R}/data/outgroups/outgroups.tsv"), delimiter="\t"))
taxa = [(r["code"], "ingroup", r["role"], r["phylum"], r["class"],
         f"{R}/work/protqc/{r['code']}/full_table.tsv", f"{R}/data/proteomes/{r['code']}.faa") for r in qc]
taxa += [(r["code"], "outgroup", "outgroup", r["group"], r["group"],
          f"{R}/work/protqc_outgroup/{r['code']}/full_table.tsv", f"{R}/data/outgroups/proteomes/{r['code']}.faa")
         for r in og]
loci = defaultdict(dict)
nsingle = {}
for code, grp, role, phy, cls, tab, faa in taxa:
    if not os.path.exists(tab):
        sys.exit(f"missing compleasm table for {code}: {tab}")
    hits = singles(tab); nsingle[code] = len(hits)
    fa = read_fa(faa) if hits else {}
    for b, sid in hits.items():
        if sid in fa:
            loci[b][code] = clean(fa[sid])
ing = {t[0] for t in taxa if t[1] == "ingroup"}
n_in = len(ing)
kept = {b: d for b, d in loci.items() if sum(c in ing for c in d) >= MINOCC * n_in}
for f in os.listdir(OUT):
    os.remove(f"{OUT}/{f}")
for b, d in sorted(kept.items()):
    with open(f"{OUT}/{b}.faa", "w") as o:
        for c, s in d.items():
            o.write(f">{c}\n{s}\n")
occ = defaultdict(int)
for d in kept.values():
    for c in d:
        occ[c] += 1
with open(f"{TD}/species.tsv", "w") as o:
    o.write("code\tset\trole\tphylum\tclass\tn_single\tn_loci\n")
    for code, grp, role, phy, cls, _, _ in taxa:
        o.write(f"{code}\t{grp}\t{role}\t{phy}\t{cls}\t{nsingle[code]}\t{occ.get(code, 0)}\n")
with open(f"{TD}/locus_occupancy.tsv", "w") as o:  # kept for 19_busco_tree compatibility
    o.write("code\tn_loci\n")
    for code, *_ in taxa:
        o.write(f"{code}\t{occ.get(code, 0)}\n")
with open(f"{TD}/loci_table.tsv", "w") as o:
    o.write("busco\tn_ingroup\tn_outgroup\tmean_len\tkept\n")
    for b, d in sorted(loci.items()):
        ni = sum(c in ing for c in d)
        ml = sum(len(s) for s in d.values()) / len(d)
        o.write(f"{b}\t{ni}\t{len(d) - ni}\t{ml:.0f}\t{int(b in kept)}\n")
print(f"ingroup {n_in} + outgroup {len(og)}; loci with any hit {len(loci)}; kept {len(kept)} "
      f"(single-copy in >= {MINOCC:.0%} of ingroup)")
low = sorted(((c, occ.get(c, 0)) for c, *_ in taxa), key=lambda x: x[1])[:10]
print("lowest-occupancy taxa:", low)
print("outgroup loci:", {r["code"]: occ.get(r["code"], 0) for r in og})
