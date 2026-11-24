#!/usr/bin/env python3
"""G1 — repeated gene-family turnover across independent deep-sea origins.

Implements PREREG_G1_gene_family_convergence.md (frozen md5 8a19c101207c7016ee1affacb1ee8943).

Everything here is a RE-ANALYSIS of CAFE5 output that already exists on disk.
CAFE5 is not re-run, no new simulation is drawn, no family is filtered by result.

Usage:  ./envs/dc/bin/python scripts/g1_family_convergence.py [--bridge]
        --bridge  also run the cross-clade eggNOG bridge (slow, streams 1,111 files)
"""
import csv
import os
import sys
import collections

import numpy as np
import pandas as pd
from ete3 import Tree

WORK = os.path.expanduser("~/deepsea_conv/work")
CAFE = os.path.join(WORK, "cafe")
OUT = os.path.join(WORK, "addon", "G1")

# *_rescue clades were prepared as re-runs of a sibling clade with an identical
# species set and no prep_report.txt; their N0 is the sibling's (see n0_path).
SIBLING_N0 = {"Mollusca_rescue": "Mollusca", "Actinopteri_m_rescue": "Actinopteri_m"}

CLADES = ["Annelida", "Chondrichthyes", "Mollusca_rescue",
          "Cnidaria", "Echinodermata", "Actinopteri_m_rescue"]
# excluded by PREREG section 6, listed so the exclusion is visible in the log
EXCLUDED = {"Arthropoda": "8/27 tsa + 3 missing provenance; LRT gain 96, matched p 1.000",
            "Porifera": "deep model broken (gain -22889); 11/15 genome-derived",
            "Actinopteri_m": "entirely NA; superseded by _rescue",
            "Mollusca": "entirely NA; superseded by _rescue"}

os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "g1_run.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


# ---------------------------------------------------------------- tree handling
def build_map(tree_path):
    """Replay CAFE5 node numbering: tips 1..n left-to-right, internals n+1.. in
    DFS postorder. Validated against prep_report.txt for 5 clades and against the
    Base_change.tab column count for all 8."""
    t = Tree(tree_path, format=1)
    leaves = t.get_leaves()
    idmap = {}
    for i, nd in enumerate(leaves, 1):
        idmap[nd] = i
    k = len(leaves)
    for nd in t.traverse("postorder"):
        if not nd.is_leaf():
            k += 1
            idmap[nd] = k
    return t, idmap


def deep_components(idmap):
    """Maximal connected components of deep-labelled (dist==2) branches.

    A branch belongs to the node it leads to; two deep branches are in the same
    component when they are parent and child. Each component is one independent
    deep-sea origin for the clade (PREREG section 3)."""
    deep = {idmap[nd] for nd in idmap if int(nd.dist) == 2}
    parent = {x: x for x in deep}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for nd in idmap:
        a = idmap[nd]
        if a in deep and nd.up is not None and idmap[nd.up] in deep:
            ra, rb = find(a), find(idmap[nd.up])
            if ra != rb:
                parent[ra] = rb
    comp = collections.defaultdict(set)
    for x in deep:
        comp[find(x)].add(x)
    return list(comp.values())


def expected_cols(idmap):
    """Column names CAFE5 writes for this numbering."""
    return ["%s<%d>" % (nd.name, i) if nd.is_leaf() else "<%d>" % i
            for nd, i in idmap.items()]


# --------------------------------------------------------------- the statistic
def col_node(c):
    """'Alvcau<12>' or '<39>' -> 12 / 39."""
    return int(c[c.rindex("<") + 1:c.rindex(">")])


def recurrence(df, comps):
    """PREREG section 4.

    e(f) = #origins with subtree-summed Delta > 0
    c(f) = #origins with subtree-summed Delta < 0
    r(f) = max(e,c);   R = sum max(0, r-1)
    Returns R, expansion-only R, contraction-only R, per-family frame."""
    node_of = {c: col_node(c) for c in df.columns}
    cols_by_node = collections.defaultdict(list)
    for c, n in node_of.items():
        cols_by_node[n].append(c)
    sub = pd.DataFrame(index=df.index)
    for j, comp in enumerate(comps):
        cs = [c for n in comp for c in cols_by_node.get(n, [])]
        sub[j] = df[cs].sum(axis=1)
    e = (sub > 0).sum(axis=1)
    c_ = (sub < 0).sum(axis=1)
    r = np.maximum(e, c_)
    R = int(np.maximum(0, r - 1).sum())
    R_exp = int(np.maximum(0, e - 1).sum())
    R_con = int(np.maximum(0, c_ - 1).sum())
    fam = pd.DataFrame({"e": e, "c": c_, "r": r})
    return R, R_exp, R_con, fam


def load_change(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df


def analyse_clade(cl):
    d = os.path.join(CAFE, cl)
    say("\n" + "=" * 78)
    say("CLADE", cl)
    say("=" * 78)

    t, idmap = build_map(os.path.join(d, "lambda_deep.nwk"))
    comps = deep_components(idmap)
    K = len(comps)
    n_tips_deep = sum(1 for nd in idmap if int(nd.dist) == 2 and nd.is_leaf())

    obs_path = os.path.join(d, "deep", "Base_change.tab")
    df = load_change(obs_path)

    # ---- numbering check: every CAFE column must be reproduced by the replay
    exp = set(expected_cols(idmap))
    got = set(df.columns)
    if exp != got:
        say("  !! COLUMN MISMATCH  missing_from_tab=%s unexpected_in_tab=%s"
            % (sorted(exp - got)[:5], sorted(got - exp)[:5]))
        raise SystemExit("numbering does not reproduce the change table for " + cl)
    say("  numbering verified: %d columns reproduced exactly" % len(got))
    say("  families %d | deep origins K=%d | deep tips %d" % (len(df), K, n_tips_deep))
    if K < 2:
        say("  !! K<2, recurrence undefined; clade dropped")
        return None

    R, R_exp, R_con, fam = recurrence(df, comps)
    rate = R / (len(df) * (K - 1))
    say("  OBSERVED  R=%d  R_expansion=%d  R_contraction=%d  rate=%.4f"
        % (R, R_exp, R_con, rate))

    # ---- nulls, both batteries
    res = {"clade": cl, "K": K, "n_fam": len(df), "R": R, "R_exp": R_exp,
           "R_con": R_con, "rate": rate, "n_nullp": 0, "n_null": 0}
    for battery, tag in (("nullp", "permutation"), ("null", "structure")):
        Rs = []
        for i in range(1, 11 if battery == "nullp" else 6):
            np_path = os.path.join(d, "%s%d" % (battery, i))
            nwk = os.path.join(d, "lambda_%s%d.nwk" % (battery, i))
            if not (os.path.isdir(np_path) and os.path.exists(nwk)):
                continue
            ct = os.path.join(np_path, "Base_change.tab")
            if not os.path.exists(ct):
                continue
            _, im2 = build_map(nwk)
            c2 = deep_components(im2)
            if len(c2) < 2:
                continue
            dfn = load_change(ct)
            if set(dfn.columns) != set(expected_cols(im2)):
                say("    %s%d column mismatch, skipped" % (battery, i))
                continue
            Rn, _, _, _ = recurrence(dfn, c2)
            Rs.append(Rn)
        key = "nullp" if battery == "nullp" else "null"
        res["n_" + key] = len(Rs)
        res[key + "_R"] = Rs
        if Rs:
            p = (1 + sum(1 for x in Rs if x >= R)) / (1 + len(Rs))
            res["p_" + key] = p
            say("  %-11s n=%2d  R null %s" % (tag, len(Rs), sorted(Rs)))
            say("  %-11s p = %d/%d = %.3f%s"
                % (tag, 1 + sum(1 for x in Rs if x >= R), len(Rs) + 1, p,
                   "   <-- AT RESOLUTION FLOOR" if p <= 1.0 / (len(Rs) + 1) else ""))
        else:
            res["p_" + key] = np.nan
            say("  %-11s no usable nulls" % tag)

    # ---- recurring families, written out for the bridge / exploratory GO step
    rec = fam[fam.r >= 2].copy()
    rec["rate"] = rec.r / K
    rec.sort_values(["r", "e", "c"], ascending=False).to_csv(
        os.path.join(OUT, "%s_recurring.tsv" % cl), sep="\t")
    say("  families recurring in >=2 origins: %d (%.2f%%)"
        % (len(rec), 100.0 * len(rec) / len(df)))

    res["recurring_ids"] = list(rec.index)
    return res


# ------------------------------------------------------------ cross-clade bridge
def most_specific_og(field):
    """eggNOG_OGs -> the OG entry at the deepest taxid present, as (name, taxid)."""
    if not field or field == "-":
        return None
    best = None
    for part in field.split(","):
        if "@" not in part:
            continue
        name, rest = part.split("@", 1)
        taxid = rest.split("|")[0]
        try:
            tid = int(taxid)
        except ValueError:
            continue
        if best is None or tid > best[1]:
            best = (name, tid)
    return best


def bridge(clade_results):
    say("\n" + "=" * 78)
    say("CROSS-CLADE BRIDGE  (PREREG section 7: eggNOG_OGs only)")
    say("=" * 78)

    # family -> member genes, for recurring families only
    fam_genes = {}
    for cl in CLADES:
        r0 = clade_results.get(cl)
        if not r0 or not r0["recurring_ids"]:
            continue
        n0 = n0_path(cl)
        want = set(r0["recurring_ids"])
        got = {}
        with open(n0) as f:
            r = csv.reader(f, delimiter="\t")
            hdr = next(r)
            spcol = hdr[3:]
            for row in r:
                if row[0] in want:
                    genes = []
                    for j, v in enumerate(row[3:]):
                        if v.strip():
                            genes += [x.strip() for x in v.split(",") if x.strip()]
                    got[row[0]] = (row[1] if len(row) > 1 else "", genes)
        fam_genes[cl] = got
        say("  %-22s recurring families loaded %d" % (cl, len(got)))

    # which eggnog files we need
    need_sp = set()
    for cl, d in fam_genes.items():
        for _, (og, genes) in d.items():
            for g in genes:
                need_sp.add(g.split("_")[0])
    say("  species whose eggNOG annotation is needed: %d" % len(need_sp))

    want_genes = set()
    for d in fam_genes.values():
        for _, (og, genes) in d.items():
            want_genes.update(genes)

    gene_og = {}
    emap = os.path.join(WORK, "eggnog")
    nfile = 0
    for sp in sorted(need_sp):
        p = os.path.join(emap, sp + ".emapper.annotations")
        if not os.path.exists(p):
            continue
        nfile += 1
        with open(p) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                f2 = line.rstrip("\n").split("\t")
                if len(f2) < 5 or f2[0] not in want_genes:
                    continue
                gene_og[f2[0]] = most_specific_og(f2[4])
    say("  annotation files read %d ; genes resolved %d / %d"
        % (nfile, len(gene_og), len(want_genes)))

    # family -> OG by 50 % support rule
    fam_og = {}
    for cl, d in fam_genes.items():
        for fam, (_, genes) in d.items():
            ann = [gene_og[g] for g in genes if g in gene_og and gene_og[g]]
            if not ann:
                continue
            cnt = collections.Counter(n for n, _ in ann)
            tot = len(genes)
            cand = [(n, c) for n, c in cnt.items() if c / tot >= 0.5]
            if not cand:
                continue
            # tie-break on the deepest taxid at which that OG was called
            tid = {n: max(t for nn, t in ann if nn == n) for n, _ in cand}
            n, c = sorted(cand, key=lambda x: (-x[1], -tid[x[0]]))[0]
            fam_og[(cl, fam)] = (n, tid[n], c / tot)

    say("  families bridged to an eggNOG OG: %d / %d"
        % (len(fam_og), sum(len(v) for v in fam_genes.values())))

    # cross-clade recurrence
    og_clades = collections.defaultdict(set)
    og_rows = []
    for (cl, fam), (og, tid, frac) in fam_og.items():
        og_clades[og].add(cl)
        og_rows.append((og, tid, cl, fam, frac))
    pd.DataFrame(og_rows, columns=["eggNOG_OG", "taxid", "clade", "family", "support"]) \
        .sort_values(["eggNOG_OG", "clade"]).to_csv(
            os.path.join(OUT, "bridged_families.tsv"), sep="\t", index=False)

    multi = {og: cs for og, cs in og_clades.items() if len(cs) >= 2}
    n_og = len(og_clades)
    n_multi = len(multi)
    say("  distinct eggNOG OGs among recurring families: %d" % n_og)
    say("  of these, recurring in >=2 clades: %d (%.1f%%)"
        % (n_multi, 100.0 * n_multi / max(1, n_og)))

    # hypergeometric expectation: if clade membership were independent given the
    # per-clade counts, P(an OG appears in >=2 clades)
    sizes = []
    for cl in CLADES:
        k = len({og for (c2, f), (og, t, fr) in fam_og.items() if c2 == cl})
        sizes.append(k)
    N = n_og
    # probability a given OG is absent from a clade, assuming random assignment
    p_multi_exp = 1.0
    for k in sizes:
        p_multi_exp *= (N - k) / N
    p_multi_exp = 1 - p_multi_exp - sum(
        (k / N) * np.prod([(N - kk) / N for kk in sizes if kk != k]) for k in sizes)
    say("  per-clade bridged-OG counts: %s" % sizes)
    say("  expected share of OGs in >=2 clades under random assignment: %.1f%%"
        % (100.0 * max(0.0, p_multi_exp)))
    say("  note: this expectation is a rough independence benchmark, not a"
        " calibrated null; the calibrated null is the CAFE label permutation.")

    return {"n_og": n_og, "n_multi": n_multi, "exp_share": p_multi_exp,
            "multi": multi, "fam_og": fam_og}


def n0_path(cl):
    """Authoritative N0 path, taken from prep_report.txt when it exists."""
    pr = os.path.join(CAFE, cl, "prep_report.txt")
    if os.path.exists(pr):
        with open(pr) as f:
            first = f.readline()
        if "\t" in first:
            p = first.split("\t")[1].strip()
            # the field reads "N0 /abs/path/N0.tsv" -- drop the leading key
            if " " in p:
                p = p.split(" ", 1)[1].strip()
            if os.path.exists(p):
                return p
    # rescue clades have no prep_report: read the path CAFE actually logged
    rep = os.path.join(CAFE, cl, "deep", "Base_report.cafe")
    if os.path.exists(rep):
        with open(rep) as f:
            for line in f:
                if "N0" in line or "hierarchical" in line.lower():
                    for tok in line.replace("\t", " ").split():
                        if tok.endswith("N0.tsv") and os.path.exists(tok):
                            return tok
    # The *_rescue clades were prepared as re-runs of a sibling clade with a
    # byte-identical species set; they carry no prep_report of their own but
    # their counts.tsv family IDs are drawn verbatim from that sibling's N0
    # (verified: Mollusca_rescue 16712/16712 and Actinopteri_m_rescue
    # 17157/17157 family IDs are all present in the sibling N0).  Borrow it.
    # Recorded in DEVIATIONS_g1.md.
    sib = SIBLING_N0.get(cl)
    if sib:
        p = n0_path(sib)
        if os.path.exists(p):
            say("  [N0] %s has no prep_report; borrowing sibling %s N0: %s"
                % (cl, sib, p))
            return p
    raise SystemExit("cannot locate N0 for " + cl)


def main():
    say("G1 family convergence — re-analysis of existing CAFE5 output")
    say("PREREG md5 8a19c101207c7016ee1affacb1ee8943 (frozen 2026-11-03)")
    say("included: %s" % ", ".join(CLADES))
    for k, v in EXCLUDED.items():
        say("excluded: %-18s %s" % (k, v))

    results = {}
    for cl in CLADES:
        try:
            r = analyse_clade(cl)
            if r:
                results[cl] = r
        except SystemExit as e:
            say("  !! %s dropped: %s" % (cl, e))

    # ---- summary table
    rows = []
    for cl, r in results.items():
        rows.append({"clade": cl, "K": r["K"], "n_fam": r["n_fam"],
                     "R": r["R"], "R_expansion": r["R_exp"],
                     "R_contraction": r["R_con"], "rate": r["rate"],
                     "n_null_perm": r["n_nullp"], "p_perm": r["p_nullp"],
                     "n_null_struct": r["n_null"], "p_struct": r["p_null"],
                     "n_recurring": len(r["recurring_ids"])})
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(OUT, "G1_summary.tsv"), sep="\t", index=False)
    say("\n" + "=" * 78)
    say("SUMMARY (PREREG C1)")
    say("=" * 78)
    say(summ.to_string(index=False))

    if "--bridge" in sys.argv:
        b = bridge(results)
        with open(os.path.join(OUT, "G1_bridge_summary.txt"), "w") as f:
            f.write("distinct eggNOG OGs among recurring families\t%d\n" % b["n_og"])
            f.write("recurring in >=2 clades\t%d\n" % b["n_multi"])
            f.write("expected share under random assignment\t%.4f\n" % b["exp_share"])
            for og, cs in sorted(b["multi"].items(), key=lambda x: -len(x[1]))[:50]:
                f.write("%s\t%s\n" % (og, ",".join(sorted(cs))))
        say("\nwrote", os.path.join(OUT, "G1_bridge_summary.txt"))

    say("\nDONE")
    LOG.close()


if __name__ == "__main__":
    main()
