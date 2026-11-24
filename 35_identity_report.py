"""Summarise work/identity/dist.tsv (mash all-vs-all) against panel taxonomy.

Outputs work/identity/
  duplicates.tsv   pairs of different species with Mash D < 0.01 (≈ >99 % ANI)
  nearest.tsv      per accession: nearest other panel genome, D, shared hashes,
                   same genus / family / order flags
  flags.tsv        accessions needing attention (duplicate, or nearest neighbour
                   outside own family although a confamilial genome with D < 0.25 exists)
Mash D saturates (~>0.3) between distant genomes, so only near neighbours are informative.
"""
import csv, os, re, collections
R = os.path.expanduser("~/deepsea_conv"); W = f"{R}/work/identity"
rows = list(csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"))
by_acc = collections.defaultdict(list)
for r in rows:
    by_acc[r["accession"]].append(r)
acc_of = lambda p: re.search(r"(GC[AF]_\d+\.\d+)", p).group(1)

D = collections.defaultdict(dict)
for l in open(f"{W}/dist.tsv"):
    a, b, d, p, sh = l.rstrip("\n").split("\t")
    a, b = acc_of(a), acc_of(b)
    if a != b:
        D[a][b] = (float(d), sh)

def tax(a):
    r = by_acc[a][0]
    return r["species"], r.get("family", ""), r.get("order", ""), r["species"].split()[0]

dup, near, flags = [], [], []
for a in sorted(D):
    sa, fa, oa, ga = tax(a)
    nb = sorted(D[a].items(), key=lambda x: x[1][0])
    b, (d, sh) = nb[0]
    sb, fb, ob, gb = tax(b)
    near.append([a, sa, fa, b, sb, fb, f"{d:.4f}", sh, int(ga == gb), int(fa == fb and fa != ""), int(oa == ob and oa != "")])
    for b2, (d2, sh2) in nb:
        if d2 >= 0.01:
            break
        if tax(b2)[0] != sa and a < b2:
            dup.append([a, sa, b2, tax(b2)[0], f"{d2:.5f}", sh2])
    conf = [d2 for b2, (d2, _) in nb if tax(b2)[1] == fa and fa]
    if fa and fb != fa and conf and min(conf) < 0.25:
        flags.append([a, sa, fa, "nearest_outside_family", f"{b} {sb} ({fb}) D={d:.4f}; best confamilial D={min(conf):.4f}"])
for x in dup:
    flags.append([x[0], x[1], tax(x[0])[1], "near_identical_to_other_species", f"{x[2]} {x[3]} D={x[4]} shared={x[5]}"])
    flags.append([x[2], x[3], tax(x[2])[1], "near_identical_to_other_species", f"{x[0]} {x[1]} D={x[4]} shared={x[5]}"])

def w(name, hdr, data):
    with open(f"{W}/{name}", "w") as f:
        f.write("\t".join(hdr) + "\n")
        for x in data:
            f.write("\t".join(map(str, x)) + "\n")
w("duplicates.tsv", ["acc1", "species1", "acc2", "species2", "mash_D", "shared"], dup)
w("nearest.tsv", ["acc", "species", "family", "nn_acc", "nn_species", "nn_family", "mash_D", "shared",
                  "same_genus", "same_family", "same_order"], near)
w("flags.tsv", ["acc", "species", "family", "flag", "detail"], flags)
print(f"genomes {len(D)}  duplicate pairs {len(dup)}  flags {len(flags)}")
