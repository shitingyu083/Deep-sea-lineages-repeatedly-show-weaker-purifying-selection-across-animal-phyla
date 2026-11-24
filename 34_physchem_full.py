#!/usr/bin/env python3
"""Full-proteome physicochemical metrics for every gene in the six G2 clades.

The v2 table (30_physchem.py) covers only compleasm BUSCO-Single genes because its
purpose was species-level comparison free of annotation route.  G2 compares genes
WITHIN a homologous family, so that constraint does not apply and BUSCO-only coverage
would leave ~700 usable families per clade instead of the ~10^4 PREREG_G2 section 1
assumed.

Output: one TSV per species in OUTDIR, gzipped, columns identical in name to
work/physchem/gene_metrics.tsv.gz for the metrics G2 uses.

Usage: python 34_physchem_full.py NPROC [OUTDIR]
"""
import glob
import gzip
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

R = os.path.expanduser("~/deepsea_conv")
AA = list("ACDEFGHIKLMNPQRSTVWY")
AIDX = {a: i for i, a in enumerate(AA)}
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
KDV = np.array([KD[a] for a in AA])

PK = dict(Cterm=3.55, Nterm=7.5, D=4.05, E=4.45, C=9.0, Y=10.0, H=5.98, K=10.0, R=12.0)
POS_AA = ["H", "K", "R"]
NEG_AA = ["D", "E", "C", "Y"]
POS_IDX = [AIDX[a] for a in POS_AA]
NEG_IDX = [AIDX[a] for a in NEG_AA]
POS_PK = np.array([PK[a] for a in POS_AA])
NEG_PK = np.array([PK[a] for a in NEG_AA])
POS_SIGN = np.array([1.0] * 3)
NEG_SIGN = np.array([1.0] * 4)

COLS = (["code", "gene", "len", "gravy", "aliphatic", "charged_DEKR", "acid_base",
         "R_RK", "IVYWREL", "PG", "aromatic", "ST", "net_charge7_per_res", "pI"]
        + ["f_" + a for a in AA])

# sequence -> integer matrix, built once per process
_TRANS = str.maketrans("ACDEFGHIKLMNPQRSTVWY",
                      "".join(chr(i) for i in range(len(AA))))


def read_fa(path):
    """name -> sequence (first token of header, as in 30_physchem.py)"""
    op = gzip.open if path.endswith(".gz") else open
    out, name, buf = {}, None, []
    with op(path, "rt") as f:
        for line in f:
            if line[0] == ">":
                if name is not None:
                    out[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            elif name is not None:
                buf.append(line.rstrip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def charge_vec(cnt, ph):
    """cnt (N,20) int array; ph (N,) float array -> net charge (N,)"""
    pos = (cnt[:, POS_IDX] / (1.0 + np.power(10.0, ph[:, None] - POS_PK[None, :]))).sum(axis=1)
    pos += 1.0 / (1.0 + np.power(10.0, ph - PK["Nterm"]))
    neg = (cnt[:, NEG_IDX] / (1.0 + np.power(10.0, NEG_PK[None, :] - ph[:, None]))).sum(axis=1)
    neg += 1.0 / (1.0 + np.power(10.0, PK["Cterm"] - ph))
    return pos - neg


def pI_vec(cnt):
    lo = np.zeros(cnt.shape[0])
    hi = np.full(cnt.shape[0], 14.0)
    for _ in range(40):
        mid = (lo + hi) / 2.0
        q = charge_vec(cnt, mid)
        up = q > 0
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    return (lo + hi) / 2.0


def metrics_for(seqs):
    """seqs: list of str -> DataFrame of the G2 columns (rows with len<30 dropped)"""
    genes, arrs = [], []
    for name, s in seqs:
        s = s.rstrip("*").upper()
        s = "".join(c for c in s if c in AIDX)
        if len(s) < 30:
            continue
        genes.append(name)
        arrs.append(np.frombuffer(s.translate(_TRANS).encode("latin-1"), dtype=np.uint8))
    if not genes:
        return None
    cat = np.concatenate(arrs)
    n = np.array([len(a) for a in arrs], dtype=np.int64)
    starts = np.zeros(len(arrs), dtype=np.int64)
    starts[1:] = np.cumsum(n)[:-1]
    cnt = np.zeros((len(arrs), 20), dtype=np.float64)
    np.add.at(cnt, (np.repeat(np.arange(len(arrs)), n), cat), 1.0)
    f = cnt / n[:, None]
    d = {"gene": genes, "len": n}
    for i, a in enumerate(AA):
        d["f_" + a] = f[:, i]
    gravy = (f @ KDV)
    d["gravy"] = gravy
    d["aliphatic"] = 100.0 * (f[:, AIDX["A"]] + 2.9 * f[:, AIDX["V"]]
                              + 3.9 * (f[:, AIDX["I"]] + f[:, AIDX["L"]]))
    d["charged_DEKR"] = f[:, [AIDX[x] for x in "DEKR"]].sum(axis=1)
    d["acid_base"] = (f[:, AIDX["D"]] + f[:, AIDX["E"]]) / np.maximum(
        f[:, AIDX["K"]] + f[:, AIDX["R"]], 1e-9)
    d["R_RK"] = f[:, AIDX["R"]] / np.maximum(f[:, AIDX["R"]] + f[:, AIDX["K"]], 1e-9)
    d["IVYWREL"] = f[:, [AIDX[x] for x in "IVYWREL"]].sum(axis=1)
    d["PG"] = f[:, AIDX["P"]] + f[:, AIDX["G"]]
    d["aromatic"] = f[:, [AIDX[x] for x in "FWY"]].sum(axis=1)
    d["ST"] = f[:, AIDX["S"]] + f[:, AIDX["T"]]
    d["net_charge7_per_res"] = charge_vec(cnt, np.full(len(genes), 7.0)) / n
    d["pI"] = pI_vec(cnt)
    return pd.DataFrame(d)


def one(code):
    pf = os.path.join(R, "data", "proteomes", code + ".faa")
    if not os.path.exists(pf):
        return code, 0
    df = metrics_for(list(read_fa(pf).items()))
    if df is None or not len(df):
        return code, 0
    df.insert(0, "code", code)
    return code, df


def dump(res, outdir):
    code, df = res
    if isinstance(df, int) or df is None or not len(df):
        return 0
    path = os.path.join(outdir, "%s.tsv.gz" % code)
    df[COLS].to_csv(path, sep="\t", index=False, compression="gzip",
                    float_format="%.10g")
    return len(df)


if __name__ == "__main__":
    NP = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(R, "work", "physchem_full")
    os.makedirs(OUT, exist_ok=True)

    codes = sorted(pd.read_csv(os.path.join(R, "work", "addon", "G2_gene_union.tsv"),
                               sep="\t")["code"])
    todo = [c for c in codes
            if not os.path.exists(os.path.join(OUT, c + ".tsv.gz"))]
    print("species %d, to do %d, procs %d" % (len(codes), len(todo), NP), flush=True)
    ntot = 0
    with Pool(NP) as p:
        for i, res in enumerate(p.imap_unordered(one, todo, chunksize=1), 1):
            ntot += dump(res, OUT)
            if i % 10 == 0 or i == len(todo):
                print("  %d/%d species, %d genes" % (i, len(todo), ntot), flush=True)
    open(os.path.join(OUT, "FULL.DONE"), "w").write("%d genes\n" % ntot)
    print("done: %d genes over %d species -> %s" % (ntot, len(codes), OUT))
