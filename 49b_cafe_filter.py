#!/usr/bin/env python3
"""CAFE5 rescue filter (2026-10-28, deviation logged in DEVIATIONS_addon.md): CAFE5 'Failed to initialize any
reasonable values' for Mollusca and Actinopteri_m. Following the CAFE5 advice, drop the families with the largest
max-min count range. Tier 1: drop the top 1% by range; tier 2: keep range <= 20.
Usage: python 49b_cafe_filter.py SRC_DIR DST_DIR TIER"""
import os, shutil, sys, glob

src, dst, tier = sys.argv[1], sys.argv[2], int(sys.argv[3])
os.makedirs(dst + "/logs", exist_ok=True)
lines = open(src + "/counts.tsv").read().splitlines()
head, body = lines[0], lines[1:]
rng = [(max(map(int, l.split("\t")[2:])) - min(map(int, l.split("\t")[2:])), l) for l in body]
if tier == 1:
    cut = sorted(r for r, _ in rng)[int(0.99 * len(rng))]
    keep = [l for r, l in rng if r < cut] if cut > 1 else [l for r, l in rng]
else:
    cut = 20
    keep = [l for r, l in rng if r <= cut]
open(dst + "/counts.tsv", "w").write("\n".join([head] + keep) + "\n")
for f in ["tree.nwk"] + [os.path.basename(x) for x in glob.glob(src + "/lambda_*.nwk")]:
    shutil.copy(src + "/" + f, dst + "/" + f)
print("tier %d: families %d -> %d (range cut %s, max range kept %d)" % (tier, len(body), len(keep), cut,
      max(r for r, l in rng if l in set(keep))))
