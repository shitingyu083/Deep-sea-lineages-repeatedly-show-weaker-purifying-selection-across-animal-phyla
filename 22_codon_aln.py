#!/usr/bin/env python3
"""Codon alignments for HOGs selected by 21_hog_select.py.
Per HOG: CDS must translate to its protein (terminal stop stripped; sequences with
internal stops or translation mismatch are dropped) -> MAFFT (auto) on protein ->
trimAl -automated1 -backtrans (codon-aware trimming) -> drop sequences with >50%
gaps -> keep if >=4 species, >=2 deep, >=2 shallow and >=300 nt.
Headers become species codes. Writes work/codon/<clade>/aln/<hog>.fas and
aln_report.tsv (hog, status, n, nD, nS, len_nt).
Usage: python 22_codon_aln.py CLADE NPROC"""
import os, subprocess, sys, tempfile
from multiprocessing import Pool
from Bio.Seq import Seq

R = os.path.expanduser("~/deepsea_conv"); E = f"{R}/envs/dc/bin"
sys.path.insert(0, f"{R}/scripts")
from importlib import import_module
lab = import_module("21_hog_select").labels()
PROT, CDS = {}, {}

def read_fa(path, keep):
    d, k = {}, None
    for l in open(path):
        if l.startswith(">"):
            k = l[1:].split()[0]; k = k if k in keep else None
            if k: d[k] = []
        elif k:
            d[k].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}

def read_aln(path):
    d, k = {}, None
    for l in open(path):
        if l.startswith(">"):
            k = l[1:].split()[0]; d[k] = []
        elif k:
            d[k].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}

def one(args):
    clade, hog, mem = args
    out = f"{R}/work/codon/{clade}/aln/{hog}.fas"
    if os.path.exists(out):
        return hog, "exists", 0, 0, 0, 0
    by_sp = {}
    for m in mem.split(","):
        sp, g = m.split(":", 1); by_sp.setdefault(sp, []).append(g)
    prot, cds = {}, {}
    for sp, gs in by_sp.items():
        g = gs[0]
        if g not in PROT or g not in CDS:
            continue
        nt = CDS[g].upper()
        if len(nt) % 3:
            continue
        aa = str(Seq(nt).translate())
        if aa.endswith("*"):
            aa, nt = aa[:-1], nt[:-3]
        q = PROT[g].rstrip("*")
        if "*" in aa or len(aa) != len(q) or sum(a != b for a, b in zip(aa[1:], q[1:])) > 0.02 * len(q):
            continue
        prot[sp], cds[sp] = aa, nt
    def counts(keys):
        return sum(lab.get(s) == "D" for s in keys), sum(lab.get(s) == "S" for s in keys)
    nD, nS = counts(prot)
    if len(prot) < 4 or nD < 2 or nS < 2:
        return hog, "too_few_valid_cds", len(prot), nD, nS, 0
    with tempfile.TemporaryDirectory(dir=f"{R}/work/codon/{clade}/tmp") as td:
        with open(f"{td}/p.fa", "w") as f:
            f.writelines(f">{s}\n{x}\n" for s, x in prot.items())
        with open(f"{td}/c.fa", "w") as f:
            f.writelines(f">{s}\n{x}\n" for s, x in cds.items())
        with open(f"{td}/p.aln", "w") as f:
            r = subprocess.run([f"{E}/mafft", "--auto", "--thread", "1", "--quiet", f"{td}/p.fa"], stdout=f)
        if r.returncode:
            return hog, "mafft_fail", len(prot), nD, nS, 0
        r = subprocess.run([f"{E}/trimal", "-in", f"{td}/p.aln", "-backtrans", f"{td}/c.fa",
                            "-automated1", "-out", f"{td}/c.aln"], capture_output=True)
        if r.returncode or not os.path.exists(f"{td}/c.aln"):
            return hog, "trimal_fail", len(prot), nD, nS, 0
        aln = read_aln(f"{td}/c.aln")
    L = len(next(iter(aln.values()))) if aln else 0
    aln = {s: x for s, x in aln.items() if L and x.count("-") / L <= 0.5}
    nD, nS = counts(aln)
    if len(aln) < 4 or nD < 2 or nS < 2 or L < 300:
        return hog, "post_trim_filter", len(aln), nD, nS, L
    with open(out + ".tmp", "w") as f:
        f.writelines(f">{s}\n{x}\n" for s, x in aln.items())
    os.replace(out + ".tmp", out)
    return hog, "ok", len(aln), nD, nS, L

def main():
    clade, npr = sys.argv[1], int(sys.argv[2])
    base = f"{R}/work/codon/{clade}"
    os.makedirs(f"{base}/aln", exist_ok=True); os.makedirs(f"{base}/tmp", exist_ok=True)
    jobs = [(clade, l.split("\t")[0], l.rstrip("\n").split("\t")[4])
            for l in list(open(f"{base}/hogs.tsv"))[1:]]
    # preload only the genes used, once; forked workers share these dicts
    need = {}
    for _, _, mem in jobs:
        for m in mem.split(","):
            sp, g = m.split(":", 1); need.setdefault(sp, set()).add(g)
    for sp, gs in need.items():
        PROT.update(read_fa(f"{R}/data/proteomes/{sp}.faa", gs))
        CDS.update(read_fa(f"{R}/data/cds/{sp}.fna", gs))
    with Pool(npr) as pool, open(f"{base}/aln_report.tsv", "a") as rep:
        for res in pool.imap_unordered(one, jobs, chunksize=4):
            rep.write("\t".join(map(str, res)) + "\n")
    print(clade, "done", len(jobs))

if __name__ == "__main__":
    main()
