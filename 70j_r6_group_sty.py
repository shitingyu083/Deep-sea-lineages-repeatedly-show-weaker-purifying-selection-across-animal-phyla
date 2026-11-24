#!/usr/bin/env python3
"""R6-2 per-origin RELAX (PREREG_r6.md + ADDENDUM_R6_2.md md5 0e9d4fe6).
prep:    HOGs of the R6-1 analysis set (all 4 R6-1 fits ok, r6_1_fits.tsv), <= 200 per clade. On the R6-1 pruned
         real / null1 trees an {FG} node becomes {Li} when its non-U tips all lie in lineage i of the production
         trees/<hog>{,.null1}.fg.txt (restricted to kept tips); other {FG} nodes lose the label; {BG} kept.
         Each group is mapped to a canonical lineage = maximal all-D (real) / all-S (null1) subtree of
         work/final_v3/species_tree.nwk pruned to the clade's species (hogs.tsv members), labels 21.labels();
         a group whose tips span > 1 canonical lineage is "merged".
         Writes OUT/trees/<c>/<hog>.<set>.grp.nwk, OUT/groups.tsv, OUT/jobs.txt, OUT/PREP.DONE.
collect: OUT/sel/<c>/<hog>.<set>.grp.relax.json -> OUT/r6_2_fits.tsv (one row per HOG x set x group).
Usage: python 70j_r6_group_sty.py prep|collect"""
import csv, json, os, re, sys
from importlib import import_module
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv"); R6 = f"{R}/work/addon/R6"; OUT = f"{R6}/R6_2"
sys.path.insert(0, f"{R}/scripts")
LAB = import_module("21_hog_select").labels()
STREE = f"{R}/work/final_v3/species_tree.nwk"
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


def prep(per=None):
    if per is None:
        fits = list(csv.DictReader(open(f"{R6}/r6_1_fits.tsv"), delimiter="\t"))
        ok = {}
        for r in fits:
            if r["arm"] in ("O", "U"):
                k = (r["clade"], r["hog"]); ok[k] = ok.get(k, True) and r["fit"] == "ok"
        hogs = sorted(k for k, v in ok.items() if v)
        per = {}
        for c, h in hogs:
            per.setdefault(c, []).append(h)
        per = {c: v[:200] for c, v in per.items()}
    os.makedirs(OUT, exist_ok=True)
    G = open(f"{OUT}/groups.tsv.tmp", "w"); J = open(f"{OUT}/jobs.txt.tmp", "w")
    G.write("clade\thog\tset\tgroup\tn_tips\tcanon\ttips\n")
    for c in sorted(per):
        can = canon(c); os.makedirs(f"{OUT}/trees/{c}", exist_ok=True)
        with open(f"{OUT}/canon.tsv", "a" if c != min(per) else "w") as f:
            if c == min(per):
                f.write("clade\tcanon\tn_species\tspecies\n")
            for w in ("D", "S"):
                for cid, cs in sorted(can[w].items()):
                    f.write(f"{c}\t{cid}\t{len(cs)}\t{','.join(sorted(cs))}\n")
        W = f"{R}/work/codon/{c}"
        for h in per[c]:
            for s, fg, want in (("real", f"{W}/trees/{h}.fg.txt", "D"), ("null1", f"{W}/trees/{h}.null1.fg.txt", "S")):
                t = Tree(open(f"{R6}/trees/{c}/{h}.{s}.nwk").read().strip(), format=1)
                keep = {strip(l.name) for l in t.get_leaves()}
                lin = fgsets(fg, keep)
                used = {}
                for n in t.traverse():
                    if n.name.endswith("{FG}"):
                        tips = {strip(x.name) for x in n.get_leaves() if LAB.get(strip(x.name)) != "U"}
                        k = [i for i, v in lin.items() if tips and tips <= v]
                        n.name = strip(n.name) + (f"{{L{k[0]}}}" if k else "")
                        if k:
                            used[k[0]] = lin[k[0]]
                if not used:
                    continue
                with open(f"{OUT}/trees/{c}/{h}.{s}.grp.nwk", "w") as f:
                    f.write(t.write(format=1) + "\n")
                for i, v in sorted(used.items()):
                    m = [cid for cid, cs in can[want].items() if v <= cs]
                    G.write("\t".join([c, h, s, f"L{i}", str(len(v)), m[0] if m else "merged", ",".join(sorted(v))]) + "\n")
                J.write(f"{c} {h} {s}\n")
        print(c, len(per[c]), "canon D", len(can["D"]), "S", len(can["S"]), flush=True)
    G.close(); J.close()
    os.replace(f"{OUT}/groups.tsv.tmp", f"{OUT}/groups.tsv"); os.replace(f"{OUT}/jobs.txt.tmp", f"{OUT}/jobs.txt")
    import random  # addendum budget cut: random 100 per clade, seed 20261031, analysis set sorted by clade, hog
    rng = random.Random(20261031)
    with open(f"{OUT}/cut100.txt", "w") as f:
        for c in sorted(per):
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


def test():  # dry run on 3 HOGs per clade with built R6-1 trees, output to /tmp (no fits needed)
    global OUT
    OUT = "/tmp/r6_2_test_sty"
    per = {}
    for c in sorted(os.listdir(f"{R6}/trees")):
        hs = sorted({f[:-len(".null1.nwk")] for f in os.listdir(f"{R6}/trees/{c}") if f.endswith(".null1.nwk")})
        per[c] = [h for h in hs if os.path.exists(f"{R6}/trees/{c}/{h}.real.nwk")][:3]
    prep(per)


{"prep": prep, "collect": collect, "test": test}[sys.argv[1]]()
