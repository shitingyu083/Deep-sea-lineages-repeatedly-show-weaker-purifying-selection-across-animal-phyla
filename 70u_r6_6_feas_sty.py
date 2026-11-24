# R6-6 feasibility counts only (no outcome column is read). Run from ~/deepsea_conv/work.
import csv, collections
rd = lambda p: csv.DictReader(open(p), delimiter="\t")
CL = ["Actinopteri_m", "Annelida", "Arthropoda", "Chondrichthyes", "Cnidaria", "Echinodermata", "Mollusca", "Porifera"]
og = {c: {r["hog"]: r["og"] for r in rd(f"func_v2/{c}/hog_annot.tsv") if r["og"]} for c in CL}
F = [(r["clade"], r["hog"], r["label"]) for r in rd("addon/R6/R6_3/r6_3_fits.tsv")]
real = {(c, h) for c, h, l in F if l == "real"}
null = {(c, h) for c, h, l in F if l == "null1"}
print("real", len(real), "null1", len(null), "null1 subset of real", len(null & real))
ogcl = collections.defaultdict(set)          # OG -> clades with a real fit (baseline donors)
for c, h in real:
    if h in og[c]: ogcl[og[c][h]].add(c)
def hasb(c, h): return h in og[c] and len(ogcl[og[c][h]] - {c}) > 0
pr = collections.Counter(); pa = collections.Counter(); bl = collections.Counter()
for c, h in real:
    if hasb(c, h): pr[c] += 1
    if (c, h) in null and hasb(c, h): pa[c] += 1
b = {(r["clade"], r["hog"]) for r in rd("addon/R6/R6_2b/sample_b.tsv")}
for c, h in b - real - null:
    if hasb(c, h): bl[c] += 1
print("clade real_withB paired_withB R6_2b_blind_withB(of blind total)")
for c in CL:
    print(c, pr[c], pa[c], bl[c], sum(1 for x in b - real - null if x[0] == c))
print("total", sum(pr.values()), sum(pa.values()), sum(bl.values()), len(b - real - null))
