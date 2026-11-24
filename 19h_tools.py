#!/usr/bin/env python3
"""Helpers for the FINAL species tree (19h_final_tree.sh, 19h_locus.sh, 19h_post.sh; 2026-10-19).
  cost W [C=1e-7] [MAXT=16]  expected MAFFT L-INS-i time per locus of W/loci, C * N^2 * L^2 seconds (N sequences,
                             L mean unaligned length; single-thread timing test on the loaded server 2026-10-19:
                             C = 4.6e-8 on a 332 x 195 locus, 7.6e-8 on 316 x 417, 2.0e-8 on 344 x 792, so 1e-7 is
                             an upper bound), and the MAFFT threads for it (ceil(time / 1 h), 1..MAXT, so no locus
                             is expected to run much over an hour) -> W/cost.tsv, most expensive first (the run
                             order; the long loci must not form the tail)
  times W                    observed per-locus MAFFT / IQ-TREE wall times, peak memory (W/log/*.time) and IQ-TREE
                             threads next to the expected cost -> W/locus_times.tsv"""
import math, os, sys


def cost(W, C=1e-7, maxt=16):
    rows = []
    for f in sorted(os.listdir(f"{W}/loci")):
        if not f.endswith(".faa"):
            continue
        n = tot = 0
        for l in open(f"{W}/loci/{f}"):
            if l.startswith(">"):
                n += 1
            else:
                tot += len(l.strip())
        L = tot / n if n else 0
        est = C * n * n * L * L
        rows.append((f[:-4], n, L, est, min(maxt, max(1, math.ceil(est / 3600)))))
    rows.sort(key=lambda r: -r[3])
    with open(f"{W}/cost.tsv.tmp", "w") as o:
        o.write("locus\tn_seq\tmean_len\test_s\tthreads\n")
        for r in rows:
            o.write(f"{r[0]}\t{r[1]}\t{r[2]:.0f}\t{r[3]:.0f}\t{r[4]}\n")
    os.replace(f"{W}/cost.tsv.tmp", f"{W}/cost.tsv")
    tot = sum(r[3] for r in rows)
    print(f"{len(rows)} loci; expected L-INS-i {tot / 3600:.0f} core-h (C = {C:g}); most expensive {rows[0][0]} "
          f"({rows[0][1]} x {rows[0][2]:.0f} aa, {rows[0][3] / 3600:.1f} h -> {rows[0][4]} threads); "
          f"{sum(r[4] > 1 for r in rows)} loci with > 1 thread")


def times(W):
    exp = {}
    for l in list(open(f"{W}/cost.tsv"))[1:]:
        f = l.rstrip("\n").split("\t")
        exp[f[0]] = f

    def rd(fn):
        try:
            v = open(fn).read().split()
            return v[-2], v[-1]          # /usr/bin/time -f "%e %M": wall seconds, peak KB (last line)
        except (OSError, IndexError):
            return "NA", "NA"

    def one(fn):
        try:
            return open(fn).read().strip() or "NA"
        except OSError:
            return "NA"
    with open(f"{W}/locus_times.tsv", "w") as o:
        o.write("locus\tn_seq\tmean_len\test_s\tmafft_threads\tmafft_s\tmafft_kb\tiqtree_threads\tiqtree_s\tiqtree_kb\n")
        for b, f in exp.items():
            ms, mk = rd(f"{W}/log/{b}.mafft.time")
            qs, qk = rd(f"{W}/log/{b}.iqtree.time")
            o.write("\t".join(f[:5] + [ms, mk, one(f"{W}/log/{b}.iqtree.threads"), qs, qk]) + "\n")
    print(f"{W}/locus_times.tsv: {len(exp)} loci")


if __name__ == "__main__":
    if sys.argv[1] == "cost":
        cost(sys.argv[2], *(float(x) for x in sys.argv[3:4]), *(int(x) for x in sys.argv[4:5]))
    elif sys.argv[1] == "times":
        times(sys.argv[2])
    else:
        sys.exit(f"unknown command {sys.argv[1]}")
