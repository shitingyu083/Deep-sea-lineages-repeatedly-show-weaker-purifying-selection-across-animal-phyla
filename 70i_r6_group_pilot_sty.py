#!/usr/bin/env python3
"""R6-2 timing pilot (PREREG_r6.md; NOT analysed): RELAX Group test mode with each deep lineage its own group
({L1}..{Lk}) and BG as reference, on production alignments of HOGs OUTSIDE the R6 sample (hylist minus nullhy).
Relabel: an {FG} node gets {Li} if all its tips lie in deep lineage i (lineages.tsv fg_lineages); BG/unlabelled kept.
Usage: python 70i_r6_group_pilot_sty.py CLADE N OUTDIR"""
import os, re, subprocess, sys, time
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv"); E = f"{R}/envs/dc/bin"
c, N, OUT = sys.argv[1], int(sys.argv[2]), sys.argv[3]
W = f"{R}/work/codon/{c}"
os.makedirs(OUT, exist_ok=True)
lin = {}
for l in list(open(f"{W}/lineages.tsv"))[1:]:
    p = l.rstrip("\n").split("\t")
    lin[p[0]] = [set(x.split(",")) for x in p[4].split(";") if x]
samp = {l.split("\t")[1] for l in list(open(f"{R}/work/addon/R6/sample.tsv"))[1:]}
null = {l.strip() for l in open(f"{W}/tiers/nullhy.txt")}
hogs = [h.strip() for h in open(f"{W}/tiers/hylist.txt") if h.strip() not in null and h.strip() not in samp]
strip = lambda n: re.sub(r"\{.*\}$", "", n)
for h in hogs[:N]:
    t = Tree(open(f"{W}/trees/{h}.hyphy.nwk").read().strip(), format=1)
    used = set()
    for n in t.traverse():
        if n.name.endswith("{FG}"):
            tips = {strip(x.name) for x in n.get_leaves()}
            k = [i for i, s in enumerate(lin.get(h, [])) if tips <= s]
            n.name = strip(n.name) + (f"{{L{k[0] + 1}}}" if k else "")
            if k:
                used.add(k[0] + 1)
    with open(f"{OUT}/{h}.grp.nwk", "w") as f:
        f.write(t.write(format=1) + "\n")
    t0 = time.time()
    r = subprocess.run([f"{E}/hyphy", "CPU=1", "relax", "--alignment", f"{W}/aln/{h}.fas", "--tree",
                        f"{OUT}/{h}.grp.nwk", "--mode", "Group mode", "--reference-group", "BG", "--models", "Minimal",
                        "--output", f"{OUT}/{h}.grp.relax.json"], capture_output=True, text=True, stdin=subprocess.DEVNULL)
    print(h, "groups", len(used), "of", len(lin.get(h, [])), "rc", r.returncode, "s", round(time.time() - t0),
          r.stdout[-300:].replace("\n", " | ") if r.returncode else "", flush=True)
