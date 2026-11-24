#!/usr/bin/env python3
"""Collect genome + proteome compleasm QC (metazoa_odb12) and gene counts into
meta/qc_table.tsv, with a provisional pass flag for downstream orthology.
Pass rule (provisional, report in Methods): proteome C = S+D >= 60 % and M <= 25 %.
D > 25 % is flagged (possible haplotig / isoform redundancy) but not excluded.
Species listed in meta/excluded_assemblies.tsv (organelle-only records < 1 Mb)
get status "excluded" with the reason in flags. The table is replaced atomically
(several loops read it while others rewrite it).
Run from ~/deepsea_conv."""
import csv, os, re

R = os.path.expanduser("~/deepsea_conv")
MIN_C, MAX_M, FLAG_D = 60.0, 25.0, 25.0
EXCL = f"{R}/meta/excluded_assemblies.tsv"


def summ(path):
    if not os.path.exists(path):
        return {}
    d = {}
    for l in open(path):
        m = re.match(r"([SDFIMN]):([\d.]+)%?,?\s*(\d+)?", l.strip())
        if m and m.group(1) != "N":
            d[m.group(1)] = float(m.group(2))
    return d


def nseq(f):
    return sum(1 for l in open(f) if l.startswith(">")) if os.path.exists(f) else 0


excluded = {r["code"]: r["reason"] for r in csv.DictReader(open(EXCL), delimiter="\t")} \
    if os.path.exists(EXCL) else {}
rows = list(csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"))
out = []
for r in rows:
    c, acc = r["code"], r["accession"]
    gq = summ(f"{R}/work/qc/{acc}/summary.txt") or summ(f"{R}/work/qc/local_{c}/summary.txt")
    pq = summ(f"{R}/work/protqc/{c}/summary.txt")
    n = nseq(f"{R}/data/proteomes/{c}.faa")
    pc = pq.get("S", 0) + pq.get("D", 0) if pq else None
    status = "pending" if not pq else (
        "pass" if pc >= MIN_C and pq.get("M", 100) <= MAX_M else "fail")
    flags = []
    if pq and pq.get("D", 0) > FLAG_D:
        flags.append("highD")
    if os.path.exists(f"{R}/work/annotfail/{c}"):
        flags.append("annot->miniprot")
    if c in excluded:
        status = "excluded"; flags.append(f"excluded:{excluded[c]}")
    out.append(dict(code=c, species=r["species"], role=r["role"], phylum=r["phylum"],
                    **{"class": r["class"]}, source=r["source"], n_genes=n,
                    genome_C=round(gq.get("S", 0) + gq.get("D", 0), 2) if gq else "",
                    prot_S=pq.get("S", ""), prot_D=pq.get("D", ""), prot_F=pq.get("F", ""),
                    prot_M=pq.get("M", ""), prot_C=round(pc, 2) if pc is not None else "",
                    status=status, flags=",".join(flags)))
tmp = f"{R}/meta/.qc_table.tsv.{os.getpid()}"
with open(tmp, "w") as f:
    w = csv.DictWriter(f, list(out[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(out)
os.replace(tmp, f"{R}/meta/qc_table.tsv")
from collections import Counter
print(Counter(o["status"] for o in out))
for o in out:
    if o["status"] == "fail" or o["flags"]:
        print(o["code"], o["role"], o["n_genes"], o["prot_C"], o["prot_D"], o["prot_M"], o["status"], o["flags"], sep="\t")
