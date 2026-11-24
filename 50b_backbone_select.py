#!/usr/bin/env python3
"""Backbone taxon set for MCMCTree (2026-10-25; PREREG_addon A6, DEVIATIONS_fig4 2026-10-25b). The full-tree run
needs 2.1 s per iteration (~49 d per chain), so MCMCTree is run on N tips chosen such that every Benton calibration
of the full tree maps to the same node with the same applied / skipped status (50_mcmctree_prep.py logic), and the
ages are carried to all tips by 50c_bb_interp.R. Deterministic selection:
 1 both root children: their best tip;
 2 every calibration with both sides sampled (33_date_lsd2.CALS): an A-tip and a B-tip in different children of the
   full-tree MRCA m (tips already chosen first, then the highest locus occupancy), so MRCA(kept A + B) = m; if the
   calibration is skipped for intruders, also one intruder under m, so it is skipped again;
 3 orders (family, class if blank) not yet represented: greedy phylogenetic diversity on the time tree (gain = branch
   length from the tip to the chosen subtree; ties: occupancy), candidates with occupancy >= MINOCC;
 4 greedy phylogenetic-diversity fill to N (same rule).
Occupancy = loci present among the NLOC loci 50_mcmctree_prep.py uses (n_taxa desc, length desc, name; >= 4 taxa),
ties = non-gap residues. If steps 1-2 exceed N the set is kept and logged.
Usage: python 50b_backbone_select.py W [--out DIR] [--n 120] [--nloc 100] [--minocc 0.2] [--timetree TT]
       python 50b_backbone_select.py --check FULL/calibrations.tsv BB/calibrations.tsv   (exit 1 if they differ)
Writes OUT/{bb.timetree.nwk, bb.timetree.lsd (copy, root age), backbone.tsv, bb_calibrations.tsv, bb_report.txt}."""
import csv, glob, os, re, shutil, sys
from importlib import import_module

R = os.path.expanduser("~/deepsea_conv")
A = sys.argv


def check(full, bb):
    """applied / skipped status and root flag per calibration must be identical"""
    rd = lambda f: {r["node_no"]: (r["is_root"], r["applied"]) for r in csv.DictReader(open(f), delimiter="\t")}
    a, b = rd(full), rd(bb)
    bad = sorted(set(a) ^ set(b)) + sorted(k for k in set(a) & set(b) if a[k] != b[k])
    for k in bad:
        print(f"calibration {k}: full {a.get(k)} backbone {b.get(k)}")
    print(f"calibrations: full {len(a)}, backbone {len(b)}, applied {sum(v[1] == 'yes' for v in a.values())} / "
          f"{sum(v[1] == 'yes' for v in b.values())}, mismatches {len(bad)}")
    return 1 if bad else 0


def clean(s):
    return re.sub(r"\[[^\]]*\]", "", s.strip())


def read_fa(fn):
    seq, cur = {}, None
    for line in open(fn):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]; seq[cur] = []
        elif cur:
            seq[cur].append(line)
    return {k: "".join(v) for k, v in seq.items()}
def occupancy(W, tips, nloc):
    """loci present / non-gap residues per tip over the loci 50_mcmctree_prep.py would select"""
    st = list(csv.DictReader(open(f"{W}/locus_stats.tsv"), delimiter="\t"))
    key = list(st[0].keys())
    loc_col, nt_col = key[0], next(c for c in key if c.startswith("n_taxa") or c == "n")
    len_col = next((c for c in key if c.startswith("len")), None)
    st.sort(key=lambda r: (-int(r[nt_col]), -int(r[len_col]) if len_col else 0, r[loc_col]))
    ts = set(tips); occ = dict.fromkeys(tips, 0); res = dict.fromkeys(tips, 0); used = []
    for r in st:
        fs = [x for x in glob.glob(f"{W}/final_aln/{r[loc_col]}.*") if not x.endswith((".log", ".time"))]
        if not fs:
            continue
        sq = {k: v.upper().replace("*", "-").replace(".", "-") for k, v in read_fa(fs[0]).items() if k in ts}
        if len(sq) < 4:
            continue
        used.append(r[loc_col])
        for k, v in sq.items():
            n = sum(c not in "-?X" for c in v)
            if n:
                occ[k] += 1; res[k] += n
        if len(used) == nloc:
            break
    return occ, res, used
def main():
    from ete3 import Tree
    sys.path.insert(0, f"{R}/scripts")
    d33 = import_module("33_date_lsd2")
    W = os.path.abspath(A[1])
    opt = lambda k, d: A[A.index(k) + 1] if k in A else d
    OUT = opt("--out", f"{R}/work/mcmctree_bb")
    N, NLOC, MINF = int(opt("--n", "120")), int(opt("--nloc", "100")), float(opt("--minocc", "0.2"))
    TT = opt("--timetree", f"{R}/work/final_v3/date/lsd.timetree.nwk")
    os.makedirs(OUT, exist_ok=True)
    tax = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}
    t = Tree(clean(open(TT).read()), format=1, quoted_node_names=True)
    tips = t.get_leaf_names(); tp = [x for x in tips if x in tax]
    leaf = {l.name: l for l in t}
    occ, res, used = occupancy(W, tips, NLOC)
    minocc = MINF * len(used)
    q = lambda x: (-occ[x], -res[x], x)                             # tip quality (smaller = better)
    why, chosen = {}, []

    def add(x, w):
        if x not in why:
            chosen.append(x); why[x] = []
        why[x].append(w)

    for i, c in enumerate(t.children):                              # 1 root children
        add(min(c.get_leaf_names(), key=q), f"root{i}")
    cal_rows = []
    for k, (fa, fb) in d33.CALS.items():                            # 2 calibrations
        a = [x for x in tp if fa(tax[x])]; b = [x for x in tp if fb(tax[x])]
        if not a or not b:
            continue
        intr = d33.check(t, a, b, tax)[0]
        m = t.get_common_ancestor(a + b)
        if m.is_leaf():
            cal_rows.append((k, len(a), len(b), len(intr), m is t, "", "", "", "m is a leaf")); continue
        ch = {x: i for i, c in enumerate(m.children) for x in c.get_leaf_names()}
        S = set(chosen); pk = lambda xs: min(xs, key=lambda x: (x not in S,) + q(x))
        byA, byB = {}, {}
        for x in a:
            byA.setdefault(ch[x], []).append(x)
        for x in b:
            byB.setdefault(ch[x], []).append(x)
        pairs = [(pk(byA[i]), pk(byB[j])) for i in sorted(byA) for j in sorted(byB) if i != j]
        if not pairs:
            cal_rows.append((k, len(a), len(b), len(intr), m is t, "", "", "", "no cross-child pair")); continue
        pa, pb = min(pairs, key=lambda p: (sum(x not in S for x in p), -min(occ[p[0]], occ[p[1]]),
                                           -(occ[p[0]] + occ[p[1]]), p))
        add(pa, f"cal{k}A"); add(pb, f"cal{k}B"); ix = ""
        if intr:
            have = [x for x in intr if x in why]
            ix = have[0] if have else min(intr, key=q)
            add(ix, f"cal{k}intruder")
        cal_rows.append((k, len(a), len(b), len(intr), m is t, pa, pb, ix, "ok"))
    n_req = len(chosen)
    # 3-4 greedy phylogenetic diversity on the time tree (the root is covered by step 1)
    cov = set()

    def cover(x):
        n = leaf[x]
        while n is not None and id(n) not in cov:
            cov.add(id(n)); n = n.up

    def gain(x):
        n, g = leaf[x], 0.0
        while n.up is not None and id(n) not in cov:
            g += n.dist; n = n.up
        return round(g, 9)

    for x in chosen:
        cover(x)
    nb = lambda v: "" if (v or "").strip().lower() in ("", "na", "nan", "-", "none", "unclassified") else v.strip()
    grp = lambda x: (nb(tax[x].get("order")) or ("f:" + nb(tax[x].get("family")) if nb(tax[x].get("family")) else "")
                     or "c:" + nb(tax[x].get("class")))
    pool = [x for x in tp if x not in why and occ[x] >= minocc]
    n_low = sum(1 for x in tp if x not in why and occ[x] < minocc)
    gcov = {grp(x) for x in chosen if x in tax}
    n_grp = len({grp(x) for x in tp})
    while len(chosen) < N:
        unc = [x for x in pool if grp(x) not in gcov]
        if not unc:
            break
        x = min(unc, key=lambda x: (-gain(x),) + q(x))
        add(x, "group:" + grp(x)); cover(x); gcov.add(grp(x)); pool.remove(x)
    while len(chosen) < N and pool:
        x = min(pool, key=lambda x: (-gain(x),) + q(x))
        add(x, "pd"); cover(x); pool.remove(x)
    gcov = {grp(x) for x in chosen if x in tax}
    t2 = t.copy()
    t2.prune(chosen, preserve_branch_length=True)
    h_full = max(t.get_distance(l) for l in t); h_bb = max(t2.get_distance(l) for l in t2)
    open(f"{OUT}/bb.timetree.nwk", "w").write(t2.write(format=5, dist_formatter="%0.10g") + "\n")
    shutil.copy(TT.replace(".nwk", ".lsd"), f"{OUT}/bb.timetree.lsd")
    g = lambda x, c: tax[x].get(c, "") if x in tax else "outgroup"
    with open(f"{OUT}/backbone.tsv", "w") as f:
        f.write("tip\tphylum\tclass\torder\tfamily\tgroup\tn_loci\tresidues\treason\n")
        for x in chosen:
            f.write("\t".join([x, g(x, "phylum"), g(x, "class"), g(x, "order"), g(x, "family"),
                               grp(x) if x in tax else "", str(occ[x]), str(res[x]), ";".join(why[x])]) + "\n")
    with open(f"{OUT}/bb_calibrations.tsv", "w") as f:
        f.write("node_no\tn_sideA\tn_sideB\tn_intruders\tis_root\ttipA\ttipB\tintruder\tstatus\n")
        for r in cal_rows:
            f.write("\t".join(map(str, r)) + "\n")
    rep = (f"backbone {len(chosen)} of {len(tips)} tips (requested {N}; required by root + calibrations {n_req}"
           f"{' > N, kept' if n_req > N else ''}); groups (order / family / class) {len(gcov)} of {n_grp}; "
           f"occupancy floor {minocc:.0f} of {len(used)} loci ({n_low} tips below, not used for the fill); "
           f"backbone occupancy min {min(occ[x] for x in chosen)} median {sorted(occ[x] for x in chosen)[len(chosen) // 2]}\n"
           f"calibrations with both sides sampled {len(cal_rows)} ({sum(r[3] > 0 for r in cal_rows)} with intruders, "
           f"intruder kept); problems {[r[0] for r in cal_rows if r[-1] != 'ok'] or 'none'}; root height full "
           f"{h_full:.6f} backbone {h_bb:.6f} (rate-scaled)\n")
    open(f"{OUT}/bb_report.txt", "w").write(rep)
    print(rep, end="", flush=True)
    oc = [set(c.get_leaf_names()) for c in t.children]              # each pruned root child inside a distinct one
    same_root = len(t2.children) == len(t.children) == 2 and sorted(
        next((i for i, s in enumerate(oc) if set(c.get_leaf_names()) <= s), -1) for c in t2.children) == [0, 1]
    return 0 if all(r[-1] == "ok" for r in cal_rows) and same_root and abs(h_full - h_bb) < 1e-3 * h_full else 1


if __name__ == "__main__":
    sys.exit(check(A[2], A[3]) if len(A) > 1 and A[1] == "--check" else main())
