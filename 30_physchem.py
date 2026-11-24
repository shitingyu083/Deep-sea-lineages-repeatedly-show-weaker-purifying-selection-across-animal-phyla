#!/usr/bin/env python3
"""Proteome / coding-sequence physicochemistry per species, on a SHARED gene set so that
annotation route (RefSeq / miniprot / TSA / author) does not drive composition differences.
Gene set = proteome-mode compleasm 'Single' hits (metazoa_odb12) of each QC-pass species.
Outputs
  work/physchem/gene_metrics.tsv.gz   code, busco, protein metrics, CDS metrics (long table)
  work/physchem/species_metrics.tsv   per-species medians over BUSCO genes + whole-proteome medians
Usage: python 30_physchem.py NPROC"""
import csv, gzip, math, os, sys
from collections import Counter
from multiprocessing import Pool
from statistics import median

R = os.path.expanduser("~/deepsea_conv")
OUT = f"{R}/work/physchem"; os.makedirs(OUT, exist_ok=True)
AA = "ACDEFGHIKLMNPQRSTVWY"
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
PK = dict(Cterm=3.55, Nterm=7.5, D=4.05, E=4.45, C=9.0, Y=10.0, H=5.98, K=10.0, R=12.0)  # EMBOSS-like
CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
SYN = {}
for c, a in CODON.items():
    SYN.setdefault(a, []).append(c)


def read_fa(f):
    op = gzip.open if f.endswith(".gz") else open
    seqs, name = {}, None
    for l in op(f, "rt"):
        if l.startswith(">"):
            name = l[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(l.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def charge(cnt, n_term, ph):
    pos = sum(cnt.get(a, 0) / (1 + 10 ** (ph - PK[a])) for a in "HKR") + 1 / (1 + 10 ** (ph - PK["Nterm"]))
    neg = sum(cnt.get(a, 0) / (1 + 10 ** (PK[a] - ph)) for a in "DECY") + 1 / (1 + 10 ** (PK["Cterm"] - ph))
    return pos - neg


def pI(cnt):
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if charge(cnt, None, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def prot_metrics(s):
    s = s.rstrip("*").upper()
    s = "".join(c for c in s if c in KD)
    n = len(s)
    if n < 30:
        return None
    cnt = Counter(s)
    f = {a: cnt.get(a, 0) / n for a in AA}
    m = {f"f_{a}": f[a] for a in AA}
    m["len"] = n
    m["gravy"] = sum(KD[c] for c in s) / n
    m["aliphatic"] = 100 * (f["A"] + 2.9 * f["V"] + 3.9 * (f["I"] + f["L"]))
    m["charged_DEKR"] = f["D"] + f["E"] + f["K"] + f["R"]
    m["acid_base"] = (f["D"] + f["E"]) / max(f["K"] + f["R"], 1e-9)
    m["R_RK"] = f["R"] / max(f["R"] + f["K"], 1e-9)
    m["IVYWREL"] = sum(f[a] for a in "IVYWREL")
    m["PG"] = f["P"] + f["G"]
    m["aromatic"] = f["F"] + f["W"] + f["Y"]
    m["ST"] = f["S"] + f["T"]
    m["net_charge7_per_res"] = charge(cnt, None, 7.0) / n
    m["pI"] = pI(cnt)
    return m


def cds_metrics(s):
    s = s.upper()
    if len(s) % 3 or len(s) < 90:
        return None
    cod = [s[i:i + 3] for i in range(0, len(s), 3)]
    if cod[-1] in ("TAA", "TAG", "TGA"):
        cod = cod[:-1]
    cod = [c for c in cod if c in CODON and CODON[c] != "*"]
    if len(cod) < 30:
        return None
    gc = sum(x in "GC" for x in s) / len(s)
    # GC3 at 4-fold degenerate sites only (codon families with 4+ synonyms sharing first two bases)
    four = [c for c in cod if CODON[c] in "AGPTV" or (c[:2] in ("CT", "CG", "TC") and CODON[c] in "LRS")]
    gc4 = sum(c[2] in "GC" for c in four) / len(four) if four else float("nan")
    gc3 = sum(c[2] in "GC" for c in cod) / len(cod)
    return dict(cds_GC=gc, GC3=gc3, GC4=gc4), Counter(cod)


def enc_wright(cc):
    """Effective number of codons (Wright 1990) from pooled codon counts."""
    Fk = {2: [], 3: [], 4: [], 6: []}
    for a, cs in SYN.items():
        k = len(cs)
        if a == "*" or k == 1:
            continue
        n = sum(cc.get(c, 0) for c in cs)
        if n < 2:
            continue
        F = (n * sum((cc.get(c, 0) / n) ** 2 for c in cs) - 1) / (n - 1)
        Fk[k].append(F)
    if not all(Fk[k] for k in (2, 4)):
        return float("nan")
    mean = lambda v: sum(v) / len(v)
    F3 = mean(Fk[3]) if Fk[3] else (mean(Fk[2]) + mean(Fk[4])) / 2
    F6 = mean(Fk[6]) if Fk[6] else mean(Fk[4])
    return min(61.0, 2 + 9 / mean(Fk[2]) + 1 / F3 + 5 / mean(Fk[4]) + 3 / F6)


def one(code):
    pf, cf = f"{R}/data/proteomes/{code}.faa", f"{R}/data/cds/{code}.fna"
    ft = f"{R}/work/protqc/{code}/full_table.tsv"
    if not (os.path.exists(pf) and os.path.exists(ft)):
        return code, None, None
    prot = read_fa(pf)
    cds = read_fa(cf) if os.path.exists(cf) else {}
    rows, busco_cod = [], Counter()
    for l in open(ft):
        if l.startswith("#"):
            continue
        x = l.rstrip("\n").split("\t")
        if len(x) < 3 or x[1] != "Single" or x[2] not in prot:
            continue
        pm = prot_metrics(prot[x[2]])
        if pm is None:
            continue
        r = dict(code=code, busco=x[0], gene=x[2], **pm)
        cm = cds_metrics(cds[x[2]]) if x[2] in cds else None
        if cm:
            r.update(cm[0]); busco_cod += cm[1]
        rows.append(r)
    whole = [m for m in (prot_metrics(s) for s in prot.values()) if m]
    sp = dict(code=code, n_busco_genes=len(rows), n_proteins=len(whole))
    keys = [k for k in rows[0] if k not in ("code", "busco", "gene")] if rows else []
    for k in keys:
        v = [r[k] for r in rows if k in r and not (isinstance(r[k], float) and math.isnan(r[k]))]
        sp[f"busco_{k}"] = median(v) if v else float("nan")
    sp["busco_ENC"] = enc_wright(busco_cod) if busco_cod else float("nan")
    for k in ("gravy", "pI", "net_charge7_per_res", "R_RK", "acid_base", "IVYWREL", "len"):
        sp[f"whole_{k}"] = median(m[k] for m in whole) if whole else float("nan")
    return code, rows, sp


if __name__ == "__main__":
    NP = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    qc = [r["code"] for r in csv.DictReader(open(f"{R}/meta/qc_table.tsv"), delimiter="\t") if r["status"] == "pass"]
    meta = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}
    with Pool(NP) as p:
        res = p.map(one, qc, chunksize=1)
    gene_cols, sps = None, []
    with gzip.open(f"{OUT}/gene_metrics.tsv.gz", "wt") as g:
        for code, rows, sp in res:
            if not rows:
                continue
            for r in rows:
                if gene_cols is None:
                    gene_cols = list(dict.fromkeys(list(r) + ["cds_GC", "GC3", "GC4"]))
                    g.write("\t".join(gene_cols) + "\n")
                g.write("\t".join(str(r.get(c, "NA")) for c in gene_cols) + "\n")
            m = meta.get(code, {})
            sp.update({k: m.get(k, "") for k in ("species", "role", "source", "phylum", "class", "order", "family",
                                                 "depth_median", "frac_gt200")})
            sps.append(sp)
    cols = list(dict.fromkeys(k for sp in sps for k in sp))
    with open(f"{OUT}/species_metrics.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for sp in sps:
            f.write("\t".join(str(sp.get(c, "NA")) for c in cols) + "\n")
    print(f"species {len(sps)} / qc-pass {len(qc)}; genes {sum(len(r[1] or []) for r in res)}")
