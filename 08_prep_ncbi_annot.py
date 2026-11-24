#!/usr/bin/env python3
"""Standardise NCBI-annotated assemblies (RefSeq / GenBank) into data/proteomes/.

Source: cds_from_genomic (gene + protein_id per CDS) and protein.faa.
Per gene keep the longest CDS whose translation matches its protein; drop
pseudo/partial-without-protein records and proteins < 30 aa.
Usage: python 08_prep_ncbi_annot.py CODE ACCESSION [CODE ACCESSION ...]
       python 08_prep_ncbi_annot.py --table meta/species_codes.tsv
Writes data/proteomes/CODE.faa, data/cds/CODE.fna, meta/idmap/CODE.tsv and
appends to meta/ncbi_prep_report.tsv
"""
import glob, gzip, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
loc = import_module("06_prep_local")  # reuse read_fa / translate / agree

R = os.path.expanduser("~/deepsea_conv")
MINLEN = 30


def find(acc, kind):
    d = f"{R}/data/ncbi/{acc}"
    flat = {"cds": "cds.fna.gz", "prot": "protein.faa.gz"}[kind]
    if os.path.exists(f"{d}/{flat}"):
        return f"{d}/{flat}"
    pat = {"cds": "cds_from_genomic.fna.gz", "prot": "protein.faa.gz"}[kind]
    hits = glob.glob(f"{d}/ncbi_dataset/data/{acc}/{pat}")
    return hits[0] if hits else None


def attr(h, k):
    m = re.search(r"\[" + k + r"=([^\]]+)\]", h)
    return m.group(1) if m else None


def prep(code, acc):
    cf, pf = find(acc, "cds"), find(acc, "prot")
    if not cf or not pf:
        return dict(code=code, acc=acc, status="missing_files")
    prot = {h.split()[0]: s for h, s in loc.read_fa(pf)}
    best, n_cds, n_bad = {}, 0, 0
    for h, s in loc.read_fa(cf):
        n_cds += 1
        if attr(h, "pseudo") == "true" or attr(h, "exception"):
            continue
        pid = attr(h, "protein_id")
        if not pid or pid not in prot:
            continue
        gx = attr(h, "db_xref") or ""
        m = re.search(r"GeneID:(\d+)", gx)
        g = m.group(1) if m else (attr(h, "locus_tag") or attr(h, "gene") or pid)
        p = prot[pid].replace("*", "")
        if len(p) < MINLEN:
            continue
        c = s.upper()
        if c[-3:] in ("TAA", "TAG", "TGA"):
            c = c[:-3]
        if len(c) % 3 or not loc.agree(p, c):
            n_bad += 1
            continue
        if g not in best or len(p) > best[g][0]:
            best[g] = (len(p), pid, p, c)
    for sub in ("proteomes", "cds"):
        os.makedirs(f"{R}/data/{sub}", exist_ok=True)
    os.makedirs(f"{R}/meta/idmap", exist_ok=True)
    with open(f"{R}/data/proteomes/{code}.faa", "w") as fa, \
         open(f"{R}/data/cds/{code}.fna", "w") as fn, \
         open(f"{R}/meta/idmap/{code}.tsv", "w") as mp:
        for i, g in enumerate(sorted(best), 1):
            L, pid, p, c = best[g]
            nid = f"{code}_{i:06d}"
            fa.write(f">{nid}\n{p}\n"); fn.write(f">{nid}\n{c}\n")
            mp.write(f"{nid}\t{pid}\t{g}\n")
    return dict(code=code, acc=acc, status="ok", n_cds=n_cds, n_bad=n_bad, n_genes=len(best))


def main():
    a = sys.argv[1:]
    if a and a[0] == "--table":
        pairs = []
        for line in open(a[1]):
            f = line.rstrip("\n").split("\t")
            if f[0] != "code" and len(f) > 3 and f[3] == "ncbi_annot":
                pairs.append((f[0], f[2]))
    else:
        pairs = list(zip(a[0::2], a[1::2]))
    rep = f"{R}/meta/ncbi_prep_report.tsv"
    new = not os.path.exists(rep)
    with open(rep, "a") as out:
        if new:
            out.write("code\tacc\tstatus\tn_cds\tn_bad\tn_genes\n")
        for code, acc in pairs:
            r = prep(code, acc)
            line = "\t".join(str(r.get(k, "")) for k in ("code", "acc", "status", "n_cds", "n_bad", "n_genes"))
            out.write(line + "\n"); out.flush(); print(line, flush=True)


if __name__ == "__main__":
    main()
