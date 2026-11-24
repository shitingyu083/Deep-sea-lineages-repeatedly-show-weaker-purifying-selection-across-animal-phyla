#!/usr/bin/env python3
"""R6-2b independent replication prep (ADDENDUM_R6_4.md part 1). Same group/canonical-lineage construction as
70j_r6_group_sty.py (copied functions, unchanged), but on production original-route alignments and trees:
work/codon/<c>/aln/<h>.fas, trees/<h>.hyphy.nwk (real), trees/<h>.null1.hyphy.nwk (null1), fg.txt files.
Sample: eligible production HOGs not in R6 sample.tsv, sorted, random.Random(20261101).sample 200 per clade.
Writes OUT=work/addon/R6/R6_2b/{sample_b.tsv,canon.tsv,groups.tsv,jobs.txt,cut100.txt,trees/,PREP.DONE}.
collect: OUT/sel/<c>/<hog>.<set>.grp.relax.json -> OUT/r6_2_fits.tsv (same columns as R6-2).
Usage: envs/dc/bin/python 70m_r6_group_b_sty.py prep|collect"""
import csv, json, os, random, re, sys
from importlib import import_module
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv"); R6 = f"{R}/work/addon/R6"; OUT = f"{R6}/R6_2b"
sys.path.insert(0, f"{R}/scripts")
LAB = import_module("21_hog_select").labels()
STREE = f"{R}/work/final_v3/species_tree.nwk"
CLADES = ["Actinopteri_m", "Annelida", "Arthropoda", "Chondrichthyes", "Cnidaria", "Echinodermata", "Mollusca",
          "Porifera"]
strip = lambda n: re.sub(r"\{.*\}$", "", n)


def canon(clade):
    hc = "Actinopteri" if clade == "Actinopteri_m" else clade
    sp = set()
    for l in list(open(f"{R}/work/codon/{hc}/hogs.tsv"))[1:]:
        sp |= {x.split(":", 1)[0] for x in l.rstrip("\n").split("\t")[4].split(",") if x}
    t = Tree(STREE, format=1)
    t.prune([l for l in t.get_leaves() if l.name in sp], preserve_branch_length=True)
    st = {}
    for n in t.traverse("postorder"):
        st[n] = ({LAB.get(n.name)} - {"U", None}) if n.is_leaf() else set().union(*(st[c] for c in n.children))
    out = {"D": [], "S": []}
    for want in "DS":
        def rec(n):
            if st[n] == {want} and (n.up is None or st[n.up] != {want}):
                tips = frozenset(l.name for l in n.iter_leaves() if LAB.get(l.name) == want)
                if tips:
                    out[want].append(tips)
                return
            for c in n.children:
                rec(c)
        rec(t)
    return {w: {f"{clade[:4]}{w}{i + 1:03d}": s for i, s in enumerate(sorted(out[w], key=sorted))} for w in "DS"}


def fgsets(path, keep):
    g = {}
    for l in open(path):
        i, tip = l.rstrip("\n").split("\t")
        tip = tip.strip("^$")
        if tip in keep:
            g.setdefault(int(i), set()).add(tip)
    return g


def files(c, h):
    W = f"{R}/work/codon/{c}"
    return {"aln": f"{W}/aln/{h}.fas", "real": f"{W}/trees/{h}.hyphy.nwk", "null1": f"{W}/trees/{h}.null1.hyphy.nwk",
            "fg_real": f"{W}/trees/{h}.fg.txt", "fg_null1": f"{W}/trees/{h}.null1.fg.txt"}


def eligible(c, h):
    f = files(c, h)
    if not all(os.path.exists(p) for p in f.values()):
        return False
    for s in ("real", "null1"):
        txt = open(f[s]).read()
        tips = re.findall(r"[(,]([^(),:{}]+)(\{[A-Z]+\})?:", txt)
        if sum(1 for _, g in tips if g == "{FG}") < 2 or sum(1 for _, g in tips if g == "{BG}") < 2:
            return False
    return True


def prep():
    os.makedirs(OUT, exist_ok=True)
    used = {(r["clade"], r["hog"]) for r in csv.DictReader(open(f"{R6}/sample.tsv"), delimiter="\t")}
    per = {}
    with open(f"{OUT}/sample_b.tsv", "w") as S:
        S.write("clade\thog\tn_eligible\n")
        for c in CLADES:
            hs = sorted(f[:-len(".fas")] for f in os.listdir(f"{R}/work/codon/{c}/aln") if f.endswith(".fas"))
            el = [h for h in hs if (c, h) not in used and eligible(c, h)]
            per[c] = sorted(random.Random(20261101).sample(el, min(200, len(el))))
            for h in per[c]:
                S.write(f"{c}\t{h}\t{len(el)}\n")
            print(c, "aln", len(hs), "eligible new", len(el), "sampled", len(per[c]), flush=True)
    G = open(f"{OUT}/groups.tsv.tmp", "w"); J = open(f"{OUT}/jobs.txt.tmp", "w")
    G.write("clade\thog\tset\tgroup\tn_tips\tcanon\ttips\n")
    with open(f"{OUT}/canon.tsv", "w") as f:
        f.write("clade\tcanon\tn_species\tspecies\n")
    for c in CLADES:
        can = canon(c); os.makedirs(f"{OUT}/trees/{c}", exist_ok=True)
        with open(f"{OUT}/canon.tsv", "a") as f:
            for w in ("D", "S"):
                for cid, cs in sorted(can[w].items()):
                    f.write(f"{c}\t{cid}\t{len(cs)}\t{','.join(sorted(cs))}\n")
        for h in per[c]:
            F = files(c, h)
            for s, want in (("real", "D"), ("null1", "S")):
                t = Tree(open(F[s]).read().strip(), format=1)
                keep = {strip(l.name) for l in t.get_leaves()}
                lin = fgsets(F["fg_" + s], keep)
                usedg = {}
                for n in t.traverse():
                    if n.name.endswith("{FG}"):
                        tips = {strip(x.name) for x in n.get_leaves() if LAB.get(strip(x.name)) != "U"}
                        k = [i for i, v in lin.items() if tips and tips <= v]
                        n.name = strip(n.name) + (f"{{L{k[0]}}}" if k else "")
                        if k:
                            usedg[k[0]] = lin[k[0]]
                if not usedg:
                    continue
                with open(f"{OUT}/trees/{c}/{h}.{s}.grp.nwk", "w") as f:
                    f.write(t.write(format=1) + "\n")
                for i, v in sorted(usedg.items()):
                    m = [cid for cid, cs in can[want].items() if v <= cs]
                    G.write("\t".join([c, h, s, f"L{i}", str(len(v)), m[0] if m else "merged", ",".join(sorted(v))]) + "\n")
                J.write(f"{c} {h} {s}\n")
        print(c, len(per[c]), "canon D", len(can["D"]), "S", len(can["S"]), flush=True)
    G.close(); J.close()
    os.replace(f"{OUT}/groups.tsv.tmp", f"{OUT}/groups.tsv"); os.replace(f"{OUT}/jobs.txt.tmp", f"{OUT}/jobs.txt")
    rng = random.Random(20261031)
    with open(f"{OUT}/cut100.txt", "w") as f:
        for c in CLADES:
            for h in sorted(rng.sample(per[c], min(100, len(per[c])))):
                f.write(f"{c} {h}\n")
    open(f"{OUT}/PREP.DONE", "w").close()


def collect():
    rows = list(csv.DictReader(open(f"{OUT}/groups.tsv"), delimiter="\t"))
    cache = {}
    with open(f"{OUT}/r6_2_fits.tsv.tmp", "w") as f:
        f.write("clade\thog\tset\tgroup\tcanon\tn_tips\tK\tp_joint\tfit\n")
        for r in rows:
            p = f"{OUT}/sel/{r['clade']}/{r['hog']}.{r['set']}.grp.relax.json"
            if p not in cache:
                try:
                    tr = json.load(open(p))["test results"]
                    cache[p] = (tr["relaxation or intensification parameter"], tr["p-value"], "ok")
                except FileNotFoundError:
                    cache[p] = ({}, "NA", "missing")
                except Exception:
                    cache[p] = ({}, "NA", "unreadable")
            K, pj, st = cache[p]
            k = K.get(r["group"], "NA")
            f.write("\t".join(map(str, [r["clade"], r["hog"], r["set"], r["group"], r["canon"], r["n_tips"], k, pj,
                                        st if k != "NA" or st != "ok" else "nogroup"])) + "\n")
    os.replace(f"{OUT}/r6_2_fits.tsv.tmp", f"{OUT}/r6_2_fits.tsv")
    print("done")


{"prep": prep, "collect": collect}[sys.argv[1]]()
