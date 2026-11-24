#!/usr/bin/env python3
"""Build per-clade OrthoFinder input dirs (symlinks to QC-passing proteomes).
Clades = focal taxon (all QC-pass species) + fixed outgroups (backbone species).
Writes work/orthofinder/<clade>/input/*.faa and work/orthofinder/clades.tsv.
Usage: python 16_orthofinder_setup.py [--allow-pending]"""
import csv, os, sys

R = os.path.expanduser("~/deepsea_conv")
# clade: (selector on species_codes row, outgroup codes)
CLADES = {
    "Actinopteri":   (lambda r: r["class"] == "Actinopteri", ["Mobhyp", "Lamflu"]),
    "Chondrichthyes": (lambda r: r["class"] in ("Chondrichthyes", "Myxini", "Hyperoartia"), ["Hypnip", "Bralan"]),
    "Mollusca":      (lambda r: r["phylum"] == "Mollusca", ["Linana", "Galcae"]),
    "Annelida":      (lambda r: r["phylum"] in ("Annelida", "Nemertea"), ["Linana", "Mme"]),
    "Arthropoda":    (lambda r: r["phylum"] == "Arthropoda", ["Pricau", "Anisim"]),
    "Cnidaria":      (lambda r: r["phylum"] == "Cnidaria", ["Syccil", "Ptyfla", "Halasi"]),
    "Porifera":      (lambda r: r["phylum"] == "Porifera", ["Mnelei", "Pocdam", "Clyhem"]),
    "Echinodermata": (lambda r: r["phylum"] == "Echinodermata", ["Ptyfla", "Bralan"]),
}
allow_pending = "--allow-pending" in sys.argv
only = [a.split("=", 1)[1] for a in sys.argv if a.startswith("--clade=")]
failed_acc = {l.split()[1] for l in open(f"{R}/data/ncbi/failed.txt") if l.startswith("FAIL")} \
    if os.path.exists(f"{R}/data/ncbi/failed.txt") else set()
qc = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/qc_table.tsv"), delimiter="\t")}
rows = list(csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"))
# panel-v3 additions are held out of the v2 runs (added later in a separate -b pass)
HOLD = set(open(f"{R}/meta/v3_hold.txt").read().split()) if os.path.exists(f"{R}/meta/v3_hold.txt") else set()
rows = [r for r in rows if r["code"] not in HOLD]
ok = lambda c: c in qc and (qc[c]["status"] == "pass" or (allow_pending and qc[c]["status"] == "pending"
                                                          and int(qc[c]["n_genes"]) > 0))
summary = []
os.makedirs(f"{R}/work/orthofinder", exist_ok=True)
notready = []
for clade, (sel, outg) in CLADES.items():
    if only and clade not in only:
        continue
    members = [r for r in rows if sel(r) or r["code"] in outg]
    pend = [r["code"] for r in members if qc.get(r["code"], {}).get("status", "pending") == "pending"
            and r["accession"] not in failed_acc and not os.path.exists(f"{R}/work/miniprot/{r['code']}/.fail")]
    if pend and not allow_pending:
        notready.append(clade); print(f"{clade}\tNOT_READY\t{len(pend)} pending: {','.join(pend[:10])}")
        continue
    focal = [r["code"] for r in rows if sel(r) and ok(r["code"])]
    og = [c for c in outg if ok(c) and c not in focal]
    d = f"{R}/work/orthofinder/{clade}/input"
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        os.remove(f"{d}/{f}")
    for c in focal + og:
        os.symlink(f"{R}/data/proteomes/{c}.faa", f"{d}/{c}.faa")
    roles = [next(r["role"] for r in rows if r["code"] == c) for c in focal]
    summary.append((clade, len(focal), roles.count("deep"), len(og), ",".join(og),
                    ",".join(c for c in (r["code"] for r in rows if sel(r)) if not ok(c))))
with open(f"{R}/work/orthofinder/clades.tsv", "a") as f:
    for s in summary:
        f.write("\t".join(map(str, s)) + "\n")
        print(*s, sep="\t")
sys.exit(2 if notready else 0)
