#!/usr/bin/env python3
"""For a clade already run on v1 (work/orthofinder/<clade>/.done), stage the v2 additions.
Exit 0 = ready: work/orthofinder/<clade>/add_v2/*.faa symlinks written (QC-pass species not in input/).
Exit 2 = not ready (members still pending QC).  Exit 3 = nothing to add.
Readiness mirrors 16_orthofinder_setup.py (pending = no QC verdict, not a failed download / miniprot)."""
import csv, os, sys

R = os.path.expanduser("~/deepsea_conv")
clade = sys.argv[1]
# 16 runs at import time; read its CLADES table by exec-ing just that block
src = open(f"{R}/scripts/16_orthofinder_setup.py").read()
ns = {}
exec(src[src.index("CLADES = {"):src.index("allow_pending")], ns)
sel, outg = ns["CLADES"][clade]
failed_acc = {l.split()[1] for l in open(f"{R}/data/ncbi/failed.txt") if l.startswith("FAIL")} \
    if os.path.exists(f"{R}/data/ncbi/failed.txt") else set()
qc = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/qc_table.tsv"), delimiter="\t")}
rows = list(csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"))
# panel-v3 additions are held out of the v2 runs (added later in a separate -b pass)
HOLD = set(open(f"{R}/meta/v3_hold.txt").read().split()) if os.path.exists(f"{R}/meta/v3_hold.txt") else set()
rows = [r for r in rows if r["code"] not in HOLD]
members = [r for r in rows if sel(r)]
pend = [r["code"] for r in members if qc.get(r["code"], {}).get("status", "pending") == "pending"
        and r["accession"] not in failed_acc and not os.path.exists(f"{R}/work/miniprot/{r['code']}/.fail")]
if pend:
    print(f"{clade}\tNOT_READY\t{len(pend)} pending: {','.join(pend[:10])}")
    sys.exit(2)
have = {f[:-4] for f in os.listdir(f"{R}/work/orthofinder/{clade}/input")}
new = [r["code"] for r in members if qc.get(r["code"], {}).get("status") == "pass" and r["code"] not in have]
if not new:
    print(f"{clade}\tNOTHING_TO_ADD"); sys.exit(3)
d = f"{R}/work/orthofinder/{clade}/add_v2"
os.makedirs(d, exist_ok=True)
for f in os.listdir(d):
    os.remove(f"{d}/{f}")
for c in new:
    os.symlink(f"{R}/data/proteomes/{c}.faa", f"{d}/{c}.faa")
roles = [next(r["role"] for r in rows if r["code"] == c) for c in new]
print(f"{clade}\tADD\t{len(new)}\tdeep={roles.count('deep')}\t{','.join(new)}")
