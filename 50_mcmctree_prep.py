#!/usr/bin/env python3
"""MCMCTree (approximate likelihood, dos Reis & Yang 2011) inputs: the robustness check of the LSD2 time tree with
soft-bounded fossil calibrations and relaxed clocks (2026-10-25).
Topology: the final rooted species tree as dated by LSD2 (work/final_v3/date/lsd.timetree.nwk; outgroups dropped).
Calibrations: the same Benton et al. 2015 node definitions as 33_date_lsd2.py (CALS, intruder check, intruded nodes
skipped), written as B(min, soft_max, 1e-300, 0.025) in units of 100 Ma (hard minimum, soft maximum, as Benton's
minima are hard and maxima soft). The root carries its own calibration if one maps to it, else RootAge < the largest
soft maximum of the calibrations (logged).
Data: the NLOC loci of the final tree (W/final_aln, W/locus_stats.tsv) with the highest taxon occupancy (ties: length),
split into NPART partitions by rate (quantiles of the mean pairwise p-distance on 300 random taxon pairs), missing
taxa as gaps (cleandata = 0).
rgene_gamma: shape 2, mean = substitutions per site per 100 Ma from the ML branch lengths of the same topology
(work/final_v3/bl.treefile: mean root-to-tip distance after midpoint-free rooting on the dated tree's root split)
divided by the LSD2 root age.
Usage: python 50_mcmctree_prep.py W [--out DIR] [--nloc 100] [--npart 4]
Writes OUT/{tree.tre, topo.nwk, aln.phy, part<k>.phy, loci.tsv, calibrations.tsv, prep_report.txt, params.sh}."""
import csv, glob, os, random, re, sys
from importlib import import_module
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv")
sys.path.insert(0, f"{R}/scripts")
A = sys.argv
W = os.path.abspath(A[1])
opt = lambda k, d: A[A.index(k) + 1] if k in A else d
OUT = opt("--out", f"{R}/work/mcmctree")
NLOC, NPART = int(opt("--nloc", "100")), int(opt("--npart", "4"))
TT = opt("--timetree", f"{R}/work/final_v3/date/lsd.timetree.nwk")
BL = opt("--bltree", f"{R}/work/final_v3/bl.treefile")
d33 = import_module("33_date_lsd2")


def clean(s):
    return re.sub(r"\[[^\]]*\]", "", s.strip())


def read_fa(fn):
    seq, cur = {}, None
    for line in open(fn):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]; seq[cur] = []
        elif cur:
            seq[cur].append(line)
    return {k: "".join(v) for k, v in seq.items()}


def pdist(a, b):
    n = d = 0
    for x, y in zip(a, b):
        if x in "-?X" or y in "-?X":
            continue
        n += 1; d += x != y
    return d / n if n >= 30 else None


def main():
    os.makedirs(OUT, exist_ok=True)
    tax = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}
    cal = {int(r["node_no"]): r for r in csv.DictReader(open(f"{R}/refs/calibrations_benton2015.tsv"), delimiter="\t")}
    t = Tree(clean(open(TT).read()), format=1, quoted_node_names=True)
    tips = t.get_leaf_names()
    # IQ-TREE/LSD2 writes the time tree in rate-scaled units (root height = |tMRCA| x rate); the root age in Ma is
    # the tMRCA reported in lsd.timetree.lsd (last 'tMRCA' line = with confidence intervals)
    lsdf = TT.replace(".nwk", ".lsd")
    tm = [float(m.group(1)) for line in open(lsdf) for m in [re.search(r"tMRCA\s+(-?[\d.]+)", line)] if m]
    root_age = abs(tm[-1])                                          # Ma
    root_h = max(t.get_distance(l) for l in t)
    for n in t.traverse():
        n.name = "" if not n.is_leaf() else n.name
    tp = [x for x in tips if x in tax]
    rows, labs = [], {}
    for k, (fa, fb) in d33.CALS.items():
        a = [x for x in tp if fa(tax[x])]; b = [x for x in tp if fb(tax[x])]
        if not a or not b:
            continue
        lo, hi = float(cal[k]["min_Ma"]), float(cal[k]["soft_max_Ma"])
        intr, desc, kids = d33.check(t, a, b, tax)
        m = t.get_common_ancestor(a + b)
        use = not intr
        if use and id(m) in labs:
            use = False; why = "same node as node %d" % labs[id(m)][0]
        else:
            why = "yes" if use else "no:intruders"
        if use:
            labs[id(m)] = (k, lo, hi)
        rows.append((k, cal[k]["node"], len(a), len(b), lo, hi, len(intr), m is t, why))
    for n in t.traverse():
        if id(n) in labs:
            n.name = f"CAL{labs[id(n)][0]}"
    s = t.write(format=8)
    bk = {k: f"'B({lo / 100:.4f},{hi / 100:.4f},1e-300,0.025)'" for k, lo, hi in labs.values()}
    s = s.replace("NoName", "")                                     # ete3 placeholder for unnamed internal nodes
    s = re.sub(r"CAL(\d+)\b", lambda m: bk[int(m.group(1))], s)
    if root_lab := labs.get(id(t)):                                 # ete3 may omit the root name
        if not s.rstrip(";").rstrip().endswith(")" + bk[root_lab[0]]):
            s = s.rstrip(";").rstrip() + bk[root_lab[0]] + ";"
    s = s.rstrip(";").rstrip() + ";"
    root_cal = id(t) in labs
    rootage = "" if root_cal else f"<{max(r[5] for r in rows if r[-1] == 'yes') / 100:.4f}"
    open(f"{OUT}/tree.tre", "w").write(f"{len(tips)} 1\n{s}\n")
    t2 = t.copy()
    for n in t2.traverse():
        if not n.is_leaf():
            n.name = ""
    open(f"{OUT}/topo.nwk", "w").write(t2.write(format=9) + "\n")
    with open(f"{OUT}/calibrations.tsv", "w") as f:
        f.write("node_no\tnode\tn_sideA\tn_sideB\tmin_Ma\tsoft_max_Ma\tn_intruders\tis_root\tapplied\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    # loci
    st = list(csv.DictReader(open(f"{W}/locus_stats.tsv"), delimiter="\t"))
    key = list(st[0].keys())
    loc_col, nt_col = key[0], next(c for c in key if c.startswith("n_taxa") or c == "n")
    len_col = next((c for c in key if c.startswith("len")), None)
    st.sort(key=lambda r: (-int(r[nt_col]), -int(r[len_col]) if len_col else 0, r[loc_col]))
    rng = random.Random(20261025)
    sel = []
    for r in st:
        fs = [x for x in glob.glob(f"{W}/final_aln/{r[loc_col]}.*") if not x.endswith((".log", ".time"))]
        if not fs:
            continue
        sq = read_fa(fs[0])
        sq = {k: v.upper().replace("*", "-").replace(".", "-") for k, v in sq.items() if k in set(tips)}
        if len(sq) < 4:
            continue
        ks = sorted(sq); ds = []
        for _ in range(300):
            a, b = rng.sample(ks, 2); d = pdist(sq[a], sq[b])
            if d is not None:
                ds.append(d)
        sel.append((r[loc_col], len(sq), len(next(iter(sq.values()))), sum(ds) / len(ds) if ds else 0.0, sq, fs[0]))
        if len(sel) == NLOC:
            break
    sel.sort(key=lambda x: x[3])
    q = [sel[i * len(sel) // NPART:(i + 1) * len(sel) // NPART] for i in range(NPART)]
    w = max(len(x) for x in tips) + 2
    with open(f"{OUT}/aln.phy", "w") as fa, open(f"{OUT}/loci.tsv", "w") as fl:
        fl.write("locus\tpartition\tn_taxa\tlength\tmean_pdist\tfile\n")
        for i, part in enumerate(q, 1):
            L = sum(x[2] for x in part)
            blk = [f"{len(tips)} {L}"]
            for tip in tips:
                blk.append(tip.ljust(w) + "".join(x[4].get(tip, "-" * x[2]) for x in part))
            txt = "\n".join(blk) + "\n\n"
            fa.write(txt); open(f"{OUT}/part{i}.phy", "w").write(txt)
            for x in part:
                fl.write(f"{x[0]}\t{i}\t{x[1]}\t{x[2]}\t{x[3]:.4f}\t{x[5]}\n")
    # rate prior
    b = Tree(clean(open(BL).read()), format=1, quoted_node_names=True)
    b.prune([x for x in b.get_leaf_names() if x in set(tips)], preserve_branch_length=True)
    side = t.children[0].get_leaf_names()
    if len(side) == 1:
        b.set_outgroup(b & side[0])
    else:
        b.set_outgroup(next(l for l in b if l.name not in side)); b.set_outgroup(b.get_common_ancestor(side))
    rtt = sum(b.get_distance(l) for l in b) / len(b)
    rate = rtt / (root_age / 100)
    with open(f"{OUT}/params.sh", "w") as f:
        f.write(f"NDATA={NPART}\nROOTAGE='{rootage}'\nRGENE='2 {2 / rate:.4f} 1'\n")
    with open(f"{OUT}/prep_report.txt", "w") as f:
        f.write(f"taxa {len(tips)}; LSD2 root age {root_age:.1f} Ma (tree height {root_h:.4f} rate-scaled); calibrations applied "
                f"{sum(r[-1] == 'yes' for r in rows)}/{len(rows)} ({', '.join(str(r[0]) for r in rows if r[-1] == 'yes')}); "
                f"root calibrated {root_cal}; RootAge {rootage or '(in tree)'}\n"
                f"loci {len(sel)} (occupancy {min(x[1] for x in sel)}-{max(x[1] for x in sel)}), partitions "
                f"{[sum(x[2] for x in p) for p in q]} sites, mean p-dist {[round(sum(x[3] for x in p) / len(p), 3) for p in q]}\n"
                f"ML mean root-to-tip {rtt:.4f} subs/site -> rate {rate:.4f} per 100 Ma; rgene_gamma 2 {2 / rate:.4f}\n")
    print(open(f"{OUT}/prep_report.txt").read(), end="", flush=True)


if __name__ == "__main__":
    main()
