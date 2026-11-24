"""HyPhy tier lists on the FINAL-eligible frame after 43d_relabel.py --apply (2026-10-20).
43c drew hylist (N_HY HOGs) and nullhy (its first N_NULLHY) at random from the HOGs eligible on the DRAFT tree.
After relabelling on FINAL (43d diff.tsv: elig_draft, elig_final) the analysis frame is the FINAL-eligible HOGs:
a listed HOG that lost eligibility leaves the lists, and a HOG eligible only on FINAL (gained) had no chance to be
drawn. Top-up: each gained HOG enters with 43c's inclusion probabilities p_hy = |hylist| / |eligible| and
p_nh = |nullhy| / |eligible|, from one uniform draw u per HOG (random.Random(20261021 + crc32(hog))); u < p_nh
implies u < p_hy, so nullhy stays inside hylist as with 43c. Every FINAL-eligible HOG then has inclusion
probability p_hy in hylist_final and p_nh in nullhy_final.
Writes work/codon/<clade>/tiers/{hylist_final,nullhy_final}.txt (sorted) and tiers/topup43g.tsv (gained HOGs:
u, in hylist_final, in nullhy_final); refuses without relabel43d/APPLIED. Prints one summary line.
Usage: python 43g_topup.py CLADE [--wdir DIR]
"""
import csv
import os
import random
import sys
import zlib

R = os.path.expanduser("~/deepsea_conv")
clade = sys.argv[1]
W = sys.argv[sys.argv.index("--wdir") + 1] if "--wdir" in sys.argv else f"{R}/work/codon/{clade}"
T, O = f"{W}/tiers", f"{W}/relabel43d"
if not os.path.exists(f"{O}/APPLIED"):
    sys.exit(f"{O}/APPLIED missing: run 43d_relabel.py --apply first")


def rd(name):
    return set(open(f"{T}/{name}.txt").read().split())


def write(name, lines):
    tmp = f"{T}/.{name}.tmp"
    with open(tmp, "w", newline="\n") as f:
        f.write("".join(x + "\n" for x in lines))
    os.replace(tmp, f"{T}/{name}")


el, hy, nh = rd("eligible"), rd("hylist"), rd("nullhy")
assert nh <= hy <= el and hy, "43c lists are not nested (nullhy <= hylist <= eligible)"
rows = list(csv.DictReader(open(f"{O}/diff.tsv"), delimiter="\t"))
el_d = {r["hog"] for r in rows if r["elig_draft"] == "1"}
el_f = {r["hog"] for r in rows if r["elig_final"] == "1"}
assert el_d == el, f"diff.tsv elig_draft ({len(el_d)}) is not tiers/eligible.txt ({len(el)})"
p_hy, p_nh = len(hy) / len(el), len(nh) / len(el)
top = []
for h in sorted(el_f - el_d):
    u = random.Random(20261021 + zlib.crc32(h.encode())).random()
    top.append((h, u, u < p_hy, u < p_nh))
hyf = sorted((hy & el_f) | {h for h, _, a, _ in top if a})
nhf = sorted((nh & el_f) | {h for h, _, _, b in top if b})
write("topup43g.tsv", ["hog\tu\tin_hylist_final\tin_nullhy_final"] +
      [f"{h}\t{u:.6f}\t{int(a)}\t{int(b)}" for h, u, a, b in top])
write("hylist_final.txt", hyf)
write("nullhy_final.txt", nhf)
print(f"43g {clade}: FINAL-eligible {len(el_f)} (DRAFT {len(el)}: lost {len(el_d - el_f)}, gained {len(top)}); "
      f"p_hy {p_hy:.4f}, p_nh {p_nh:.4f}; hylist_final {len(hyf)} = {len(hy & el_f)} of hylist + "
      f"{sum(a for _, _, a, _ in top)} top-up; nullhy_final {len(nhf)} = {len(nh & el_f)} of nullhy + "
      f"{sum(b for _, _, _, b in top)} top-up")
