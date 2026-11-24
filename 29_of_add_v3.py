#!/usr/bin/env python3
"""Stage the panel-v3 species of a clade for an OrthoFinder -b pass on top of its finished run (2026-10-20).
Panel-v3 species (meta/v3_hold.txt) were held out of the v1/v2 runs (16_orthofinder_setup.py, 29_of_add_v2.py)
but the QC-passing ones are in the final species tree (tree_v3). Members = the clade selector of
16_orthofinder_setup.py (focal taxon only, as in 29_of_add_v2.py); staged = QC-pass HOLD members missing from
the base run (BASE_WORKINGDIR/SpeciesIDs.txt) -> work/orthofinder/<clade>/add_v3/*.faa symlinks.
QC-pass members missing from the base run that are not HOLD species are reported, never staged.
Exit 0 = staged, 3 = nothing to add, 1 = error.  --dry: print the list only.
Usage: python 29_of_add_v3.py CLADE BASE_WORKINGDIR [--dry]"""
import csv, os, sys

R = os.path.expanduser("~/deepsea_conv")
clade, wd = sys.argv[1], sys.argv[2]
dry = "--dry" in sys.argv
# 16 runs at import time; read its CLADES table by exec-ing just that block (as 29_of_add_v2.py does)
src = open(f"{R}/scripts/16_orthofinder_setup.py").read()
ns = {}
exec(src[src.index("CLADES = {"):src.index("allow_pending")], ns)
sel, outg = ns["CLADES"][clade]
qc = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/qc_table.tsv"), delimiter="\t")}
rows = list(csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"))
HOLD = set(open(f"{R}/meta/v3_hold.txt").read().split())
# SpeciesIDs.txt lines "0: Code.faa" ('#' = species excluded from that run)
have = set()
for line in open(f"{wd}/SpeciesIDs.txt"):
    if ": " in line and not line.startswith("#"):
        have.add(line.split(": ", 1)[1].strip().rsplit(".", 1)[0])
miss = [r for r in rows if sel(r) and qc.get(r["code"], {}).get("status") == "pass" and r["code"] not in have]
odd = [r["code"] for r in miss if r["code"] not in HOLD]
new = [r for r in miss if r["code"] in HOLD]
print(f"{clade}\tbase species={len(have)}\tbase={wd}")
for r in new:
    q = qc[r["code"]]
    print(f"  add\t{r['code']}\t{r['species']}\t{r['role']}\t{r['class']}\t{r['order']}\t{r['source']}"
          f"\tprot_C={q.get('prot_C')}\tprot_M={q.get('prot_M')}")
if odd:
    print(f"  WARNING: QC-pass members missing from the base run but not HOLD (not staged): {','.join(odd)}")
bad = [r["code"] for r in new if not os.path.exists(f"{R}/data/proteomes/{r['code']}.faa")]
if bad:
    sys.exit(f"missing proteomes: {','.join(bad)}")
if not new:
    print(f"{clade}\tNOTHING_TO_ADD")
    sys.exit(3)
roles = [r["role"] for r in new]
print(f"{clade}\tADD\t{len(new)}\tdeep={roles.count('deep')}\t{','.join(r['code'] for r in new)}")
if dry:
    sys.exit(0)
d = f"{R}/work/orthofinder/{clade}/add_v3"
os.makedirs(d, exist_ok=True)
for f in os.listdir(d):
    os.remove(f"{d}/{f}")
for r in new:
    os.symlink(f"{R}/data/proteomes/{r['code']}.faa", f"{d}/{r['code']}.faa")
