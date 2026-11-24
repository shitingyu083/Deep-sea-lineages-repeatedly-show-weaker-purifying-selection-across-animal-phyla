#!/usr/bin/env python3
"""54a (A4, PREREG_addon.md): per (HOG, set) features of the selection tests, joined to 25b's sel_results_v3.tsv.
Per row of work/codon/<C>/sel_results_v3.tsv (8 primary clades of 52's CL8): the FG / BG tips of
trees/<hog><tag>.hyphy.nwk (name{FG} / name{BG}; unlabelled tips, e.g. the unpaired deep lineages of Actinopteri_m,
are in neither), their annotation route (work/physchem/species_metrics.tsv source: ncbi_annot / local / extra_local =
annot, else miniprot / tsa; as 44 and addon_lib.R), the mean per-sequence gap fraction ('-', N, ?, X) of those tips in
the codon alignment aln/<hog>.fas, the tree tip count and the alignment length. fg_match = 1 when the tips of
trees/<hog><tag>.fg.txt (the CSUBST foreground) equal the HyPhy FG tips.
Usage: python 54a_sel_features.py [--clades "C1 C2 ..."] [--sub norescue] [--out work/addon/A4] [--np 4]
Writes OUT/sel_features.tsv (sel_features_<sub>.tsv with --sub; the sel_results of W/<sub>/, trees and alignments of
W) and prints per clade x set type (real, null1, null2+) the mean route and gap features."""
import csv
import os
import re
import sys
from multiprocessing import Pool

R = os.environ.get("DSC_ROOT", os.path.expanduser("~/deepsea_conv"))
A = sys.argv
opt = lambda k, d: A[A.index(k) + 1] if k in A else d
CL8 = "Chondrichthyes Porifera Arthropoda Echinodermata Cnidaria Annelida Mollusca Actinopteri_m"
CLADES = opt("--clades", CL8).split()
SUB = opt("--sub", "")
OUT = opt("--out", f"{R}/work/addon/A4")
NP = int(opt("--np", "4"))
ANNOT = {"ncbi_annot", "local", "extra_local"}
FEAT = ["n_tree", "n_aln", "aln_len", "fg_n", "bg_n", "un_n", "fg_nonannot", "bg_nonannot", "fg_mp", "bg_mp",
        "fg_tsa", "bg_tsa", "fg_gap", "bg_gap", "gapdiff", "fg_match", "n_noroute", "n_noaln"]


def routes():
    with open(f"{R}/work/physchem/species_metrics.tsv") as f:
        return {r["code"]: ("annot" if r["source"] in ANNOT else r["source"]) for r in csv.DictReader(f, delimiter="\t")}


RT = routes()


def tagged(nwk, lab):
    return {t for t in re.findall(r"([^(),:;{}\s]+)\{" + lab + r"\}", nwk) if not re.fullmatch(r"Node\d+", t)}


def tips(nwk):  # tip labels follow '(' or ','; internal node labels follow ')'
    return {re.sub(r"\{[^}]*\}$", "", t) for t in re.findall(r"[(,]\s*([^(),:;\s]+)", nwk)}


def read_fa(fn):
    d, k = {}, None
    for line in open(fn):
        if line.startswith(">"):
            k = line[1:].split()[0]; d[k] = []
        elif k:
            d[k].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}


def gapfrac(s):
    return (s.count("-") + s.count("N") + s.count("n") + s.count("?") + s.count("X") + s.count("x")) / len(s)


def mean(v):
    v = [x for x in v if x is not None]
    return sum(v) / len(v) if v else None


def one(job):
    C, hog, sets = job
    W = f"{R}/work/codon/{C}"
    try:
        aln = read_fa(f"{W}/aln/{hog}.fas")
    except OSError:
        aln = {}
    gap = {k: gapfrac(s) for k, s in aln.items() if s}
    L = max((len(s) for s in aln.values()), default=0)
    res = []
    for s in sets:
        tag = "" if s == "real" else "." + s
        try:
            nwk = open(f"{W}/trees/{hog}{tag}.hyphy.nwk").read()
        except OSError:
            res.append({}); continue
        fg, bg, al = tagged(nwk, "FG"), tagged(nwk, "BG"), tips(nwk)
        fgf = f"{W}/trees/{hog}{tag}.fg.txt"
        cs = {x.split("\t")[1].strip()[1:-1] for x in open(fgf) if "\t" in x} if os.path.exists(fgf) else None
        rt = lambda t, k: None if t not in RT else float(RT[t] in k)
        o = dict(n_tree=len(al), n_aln=len(aln), aln_len=L, fg_n=len(fg), bg_n=len(bg), un_n=len(al - fg - bg),
                 fg_nonannot=mean([rt(t, ("miniprot", "tsa")) for t in fg]),
                 bg_nonannot=mean([rt(t, ("miniprot", "tsa")) for t in bg]),
                 fg_mp=mean([rt(t, ("miniprot",)) for t in fg]), bg_mp=mean([rt(t, ("miniprot",)) for t in bg]),
                 fg_tsa=mean([rt(t, ("tsa",)) for t in fg]), bg_tsa=mean([rt(t, ("tsa",)) for t in bg]),
                 fg_gap=mean([gap.get(t) for t in fg]), bg_gap=mean([gap.get(t) for t in bg]),
                 fg_match=None if cs is None else int(cs == fg),
                 n_noroute=sum(t not in RT for t in al), n_noaln=sum(t not in gap for t in al))
        o["gapdiff"] = None if o["fg_gap"] is None or o["bg_gap"] is None else o["fg_gap"] - o["bg_gap"]
        res.append(o)
    return C, hog, res


def fmt(v):
    if v is None:
        return "NA"
    return f"{v:.6g}" if isinstance(v, float) else str(v)


def main():
    os.makedirs(OUT, exist_ok=True)
    fo = f"{OUT}/sel_features{'_' + SUB if SUB else ''}.tsv"
    hdr, rows, jobs = None, {}, []
    for C in CLADES:
        W = f"{R}/work/codon/{C}"
        fn = f"{W}/{SUB}/sel_results_v3.tsv" if SUB else f"{W}/sel_results_v3.tsv"
        with open(fn) as f:
            rd = csv.reader(f, delimiter="\t")
            h = next(rd)
            if hdr is None:
                hdr = h
            elif h != hdr:
                sys.exit(f"54a: {fn}: header differs from {CLADES[0]}")
            for r in rd:
                rows.setdefault((C, r[0]), []).append(r)
        jobs += [(C, hog, [r[1] for r in rs]) for (c, hog), rs in rows.items() if c == C]
    pool = Pool(NP) if NP > 1 else None
    it = pool.imap(one, jobs, chunksize=16) if pool else map(one, jobs)
    agg = {}
    with open(fo, "w", newline="\n") as f:
        f.write("\t".join(["clade"] + hdr + FEAT) + "\n")
        for C, hog, res in it:
            for r, o in zip(rows[(C, hog)], res):
                f.write("\t".join([C] + r + [fmt(o.get(k)) for k in FEAT]) + "\n")
                st = "real" if r[1] == "real" else ("null1" if r[1] == "null1" else "null2+")
                a = agg.setdefault((C, st), [0, [], [], [], []])
                a[0] += 1
                for i, k in enumerate(("fg_nonannot", "bg_nonannot", "gapdiff", "fg_match")):
                    a[i + 1].append(o.get(k))
    if pool:
        pool.close(); pool.join()
    print("clade\tset\tn\tfg_nonannot\tbg_nonannot\tgapdiff\tfg_match")
    for (C, st), a in sorted(agg.items()):
        print("\t".join([C, st, str(a[0])] + [fmt(mean(v)) for v in a[1:]]))
    print(f"54a: wrote {fo}")


if __name__ == "__main__":
    main()
