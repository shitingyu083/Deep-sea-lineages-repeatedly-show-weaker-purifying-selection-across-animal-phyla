#!/usr/bin/env python
"""R3 (PREREG_nature.md): collect RTF1 RELAX/BUSTED-PH (real + null1 labelling) and run the fish RTF1
amino-acid site scan (hadal/abyssal fishes vs all shallow fishes; 1000 random 4-shallow-species null)."""
import os, json, glob, random, csv
B = os.path.expanduser("~/deepsea_conv"); R = f"{B}/work/addon/R3/rtf1"
CODE = dict(zip([a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSSRRRRVVVVAAAADDEEGGGG"))

def g(d, *ks):
    for k in ks:
        if not isinstance(d, dict): return ""
        d = d.get(k, {})
    return "" if d == {} else d

rows = []
for f in sorted(glob.glob(f"{R}/*.json")):
    C, h1, h2, lab, test, _ = os.path.basename(f).split(".")
    try: j = json.load(open(f))
    except Exception: rows.append([C, f"{h1}.{h2}", lab, test, "PARSE_FAIL", "", "", "", ""]); continue
    if test == "relax":
        rd = g(j, "fits", "RELAX alternative", "Rate Distributions")
        rows.append([C, f"{h1}.{h2}", lab, test, g(j, "test results", "p-value"),
                     g(j, "test results", "relaxation or intensification parameter"), "", "", json.dumps(rd)[:300]])
    else:
        rows.append([C, f"{h1}.{h2}", lab, test, g(j, "test results", "p-value"), "",
                     g(j, "Background selection test results", "p-value"),
                     g(j, "Comparative selection test results", "p-value"), ""])
with open(f"{R}/rtf1_hyphy.tsv", "w") as fo:
    fo.write("clade\thog\tlabelling\ttest\tp\tK\tbg_p\tdiff_p\trate_dist\n")
    for r in rows: fo.write("\t".join(map(str, r)) + "\n")

# ---- fish site scan ----
lab = {r["code"]: r for r in csv.DictReader(open(f"{B}/meta/labels_v2.tsv"), delimiter="\t")}
aln = {}; k = None
for l in open(f"{B}/work/codon/Actinopteri/aln/N0.HOG0011415.fas"):
    l = l.strip()
    if l.startswith(">"): k = l[1:].split()[0]; aln[k] = []
    elif k: aln[k].append(l)
aa = {}
for t, v in aln.items():
    s = "".join(v).upper(); aa[t] = "".join(CODE.get(s[i:i + 3], "-") if "-" not in s[i:i + 3] else "-" for i in range(0, len(s) - 2, 3))
want = ["pseudoliparis", "notoliparis", "bassozetus", "bathysaurus"]
hadal = [t for t in aa if t in lab and any(w in lab[t]["species"].lower() for w in want)]
shal = [t for t in aa if t in lab and lab[t]["lab"] == "S"]
deep = [t for t in aa if t in lab and lab[t]["lab"] == "D"]
L = min(len(s) for s in aa.values())

def scan(focal, ref):
    hits = []
    for i in range(L):
        st = [aa[t][i] for t in focal if aa[t][i] not in "-X*"]
        refst = set(aa[t][i] for t in ref)
        for a in set(st):
            if st.count(a) >= 2 and a not in refst: hits.append((i + 1, a, st.count(a))); break
    return hits

obs = scan(hadal, shal)
rng = random.Random(20261028); null = []
for _ in range(1000):
    f4 = rng.sample(shal, 4); null.append(len(scan(f4, [t for t in shal if t not in f4])))
p = (1 + sum(x >= len(obs) for x in null)) / 1001.0
# depth-matched control: random 4 non-hadal deep fishes vs all shallow
dnull = []
dpool = [t for t in deep if t not in hadal]
for _ in range(1000 if len(dpool) >= 4 else 0):
    dnull.append(len(scan(rng.sample(dpool, 4), shal)))
with open(f"{R}/fish_sitescan.txt", "w") as fo:
    fo.write(f"aln_len_aa\t{L}\nhadal\t{','.join(hadal)}\t{[lab[t]['species'] for t in hadal]}\n")
    fo.write(f"n_shallow\t{len(shal)}\tn_deep\t{len(deep)}\n")
    fo.write(f"observed_sites\t{len(obs)}\t{obs}\n")
    fo.write(f"null_shallow4_mean\t{sum(null)/len(null):.2f}\tq95\t{sorted(null)[949]}\tp_ge\t{p:.4f}\n")
    if dnull: fo.write(f"null_deep4_mean\t{sum(dnull)/len(dnull):.2f}\tq95\t{sorted(dnull)[949]}\tp_ge\t{(1+sum(x>=len(obs) for x in dnull))/1001:.4f}\n")
    fo.write("note\tnull focal sets are compared with the remaining n_shallow-4 shallow fishes; observed with all n_shallow\n")
print(open(f"{R}/fish_sitescan.txt").read())
