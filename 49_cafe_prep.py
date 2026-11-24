#!/usr/bin/env python3
"""CAFE5 inputs per clade (Fig. 4c; 2026-10-25).
Counts: genes per species of every N0 hierarchical orthogroup of the clade OrthoFinder run, outgroup species (21_hog_select.OUTG) excluded.
Filters (CAFE5 tutorial): present in >= 2 species and on both sides of the root (families absent at the root make the
birth-death likelihood degenerate); families with >= 100 copies in any species go to counts_large.tsv (fit with the
lambda of the main set, not used for the test).
Tree: the LSD2 time tree of the final species tree (work/final_v3/date/lsd.timetree.nwk) pruned to the clade species,
annotations stripped, branch lengths floored at 0.01 (CAFE5 requires positive lengths; ultrametric up to rounding).
Lambda trees (-y): deep = 2 for every branch whose labelled descendant tips are all deep (D) (the branch leading to and
within each maximal deep subtree), else 1. Null lambda trees null1..K: the same count of maximal all-shallow subtrees,
size-matched greedily to the deep subtrees (nearest tip count, without replacement, seed 20261025 + k), labelled 2.
The deep two-rate model is compared with the one-rate model (LRT, df 1) and its likelihood gain with the gains of the
null shallow labellings (empirical calibration, as the selection tests).
Usage: python 49_cafe_prep.py CLADE [--out DIR] [--tree FILE] [--nnull 5]"""
import os, random, sys
from importlib import import_module

R = os.path.expanduser("~/deepsea_conv")
sys.path.insert(0, f"{R}/scripts")
A = sys.argv
clade = A[1]
opt = lambda k, d: A[A.index(k) + 1] if k in A else d
OUT = opt("--out", f"{R}/work/cafe/{clade}")
TREE = opt("--tree", f"{R}/work/final_v3/date/lsd.timetree.nwk")
NNULL = int(opt("--nnull", "5"))
NPERM = int(opt("--nperm", "10"))
hs = import_module("21_hog_select")


class N:
    __slots__ = ("ch", "name", "bl", "up")

    def __init__(s):
        s.ch, s.name, s.bl, s.up = [], "", 0.0, None


def parse(t):
    t = t.strip().rstrip(";")
    out = []
    i = 0
    while i < len(t):                       # drop [...] annotations
        if t[i] == "[":
            j = t.index("]", i); i = j + 1; continue
        out.append(t[i]); i += 1
    t = "".join(out)
    root = N(); cur = root; i = 0; tok = ""

    def flush(node, tok):
        if ":" in tok:
            a, b = tok.split(":", 1); node.name = a.strip("'\""); node.bl = float(b) if b else 0.0
        elif tok:
            node.name = tok.strip("'\"")
    while i < len(t):
        c = t[i]
        if c == "(":
            n = N(); n.up = cur; cur.ch.append(n); cur = n
        elif c == ",":
            flush(cur, tok); tok = ""
            n = N(); n.up = cur.up; cur.up.ch.append(n); cur = n
        elif c == ")":
            flush(cur, tok); tok = ""; cur = cur.up
        else:
            tok += c
        i += 1
    flush(cur, tok)
    return root


def tips(n):
    return [n] if not n.ch else [x for c in n.ch for x in tips(c)]


def prune(n, keep):
    if not n.ch:
        return n if n.name in keep else None
    kids = [k for k in (prune(c, keep) for c in n.ch) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        kids[0].bl += n.bl; kids[0].up = n.up; return kids[0]
    n.ch = kids
    for k in kids:
        k.up = n
    return n


def nwk(n, val=None, internal=False):
    s = n.name if not n.ch else "(" + ",".join(nwk(c, val, internal) for c in n.ch) + ")"
    if n.up is None:
        return s
    if val is not None:
        return f"{s}:{val[id(n)]}"
    return f"{s}:{max(n.bl, 0.01):.4f}"


def nodes(n):
    yield n
    for c in n.ch:
        yield from nodes(c)


def maximal(root, lab, want):
    """Maximal subtrees (not the root) whose labelled tips are all `want` (>= 1 labelled tip)."""
    res = []

    def ok(n):
        L = [lab.get(t.name, "U") for t in tips(n)]
        L = [x for x in L if x in "DS"]
        return bool(L) and all(x == want for x in L)

    def walk(n):
        if n.up is not None and ok(n):
            res.append(n); return
        for c in n.ch:
            walk(c)
    walk(root)
    return res


def lam(root, subs):
    v = {id(n): 1 for n in nodes(root)}
    for s in subs:
        for n in nodes(s):
            v[id(n)] = 2
    return v


def main():
    os.makedirs(OUT, exist_ok=True)
    lab = hs.labels()
    res = hs.latest_results(clade.split("_")[0])
    n0 = f"{res}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
    with open(n0) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        sp_raw = hdr[3:]
        sp = [s.split(".")[0] for s in sp_raw]
        root = parse(open(TREE).read())
        # LSD2 time tree is rate-scaled (height = |tMRCA| x rate): rescale branch lengths to Ma
        lsdf = TREE.replace(".nwk", ".lsd")
        if os.path.exists(lsdf):
            import re
            tm = [float(m.group(1)) for line in open(lsdf) for m in [re.search(r"tMRCA\s+(-?[\d.]+)", line)] if m]

            def height(n):
                return 0.0 if not n.ch else max(c.bl + height(c) for c in n.ch)
            scale = abs(tm[-1]) / height(root)
            for n in nodes(root):
                n.bl *= scale
        else:
            scale = 1.0
        intree = {t.name for t in tips(root)}
        og = set(hs.OUTG.get(clade.split("_")[0], []))
        keep = [s for s in sp if s in intree and s not in og]       # outgroups out: the root split is the ingroup root
        miss = [s for s in sp if s not in intree and s not in og]
        root = prune(root, set(keep)); root.up = None; root.bl = 0.0
        order = [t.name for t in tips(root)]
        L, R2 = {t.name for t in tips(root.ch[0])}, {t.name for t in tips(root)} - {t.name for t in tips(root.ch[0])}
        col = {s: i for i, s in enumerate(sp)}
        nk = nl = n_low = n_root = 0
        fo, fl = open(f"{OUT}/counts.tsv", "w"), open(f"{OUT}/counts_large.tsv", "w")
        for f_ in (fo, fl):
            f_.write("Desc\tFamily ID\t" + "\t".join(order) + "\n")
        for line in f:
            r = line.rstrip("\n").split("\t")
            c = {s: (len(r[3 + col[s]].split(", ")) if 3 + col[s] < len(r) and r[3 + col[s]].strip() else 0) for s in order}
            pres = [s for s in order if c[s]]
            if len(pres) < 2:
                n_low += 1; continue
            if not (set(pres) & L and set(pres) & R2):
                n_root += 1; continue
            row = "(null)\t" + r[0] + "\t" + "\t".join(str(c[s]) for s in order) + "\n"
            if max(c.values()) >= 100:
                fl.write(row); nl += 1
            else:
                fo.write(row); nk += 1
        fo.close(); fl.close()
    open(f"{OUT}/tree.nwk", "w").write(nwk(root) + ";\n")
    deep = maximal(root, lab, "D")
    shal = maximal(root, lab, "S")
    open(f"{OUT}/lambda_deep.nwk", "w").write(nwk(root, lam(root, deep)) + ";\n")
    sz = lambda n: len(tips(n))
    made = 0
    for k in range(1, NNULL + 1):
        rng = random.Random(20261025 + k)
        pool = shal[:]; rng.shuffle(pool); pick = []
        for d in sorted(deep, key=sz, reverse=True):
            if not pool:
                break
            b = min(pool, key=lambda n: abs(sz(n) - sz(d)))
            pool.remove(b); pick.append(b)
        if len(pick) == len(deep):
            open(f"{OUT}/lambda_null{k}.nwk", "w").write(nwk(root, lam(root, pick)) + ";\n"); made += 1
    # tip-label permutation nulls (all clades; the only null where shallow subtrees are too few, e.g. the
    # deep-heavy fish): D/S labels shuffled among the labelled tips, maximal deep subtrees of the permuted labels
    labd = [s for s in order if lab.get(s) in ("D", "S")]
    nperm = 0
    for k in range(1, NPERM + 1):
        rng = random.Random(20261125 + k)
        v = [lab[s] for s in labd]; rng.shuffle(v)
        pl = dict(zip(labd, v))
        sub = maximal(root, pl, "D")
        open(f"{OUT}/lambda_nullp{k}.nwk", "w").write(nwk(root, lam(root, sub)) + ";\n"); nperm += 1
    with open(f"{OUT}/prep_report.txt", "w") as f:
        f.write(f"clade {clade}\tN0 {n0}\ntime scale x{scale:.1f} (Ma); outgroups excluded {len(og & set(sp))}\n"
                f"species {len(sp)} in_tree {len(keep)} missing {','.join(miss) or '-'}\n"
                f"families main {nk} large(>=100) {nl} dropped <2 species {n_low} dropped absent-at-root {n_root}\n"
                f"deep maximal subtrees {len(deep)} (tips {sum(map(sz, deep))}); shallow maximal subtrees {len(shal)}; "
                f"null lambda trees {made}/{NNULL}; permutation null trees {nperm}\n"
                f"labels D {sum(lab.get(s) == 'D' for s in order)} S {sum(lab.get(s) == 'S' for s in order)} "
                f"other {sum(lab.get(s) not in ('D', 'S') for s in order)}\n")
    print(open(f"{OUT}/prep_report.txt").read(), end="", flush=True)


if __name__ == "__main__":
    main()
