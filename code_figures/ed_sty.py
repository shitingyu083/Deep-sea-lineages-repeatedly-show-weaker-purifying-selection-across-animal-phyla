"""Extended Data figures 1-8 (Nature spec, 180 mm). Usage: python ed_sty.py [1 2 ...] (default all).

ED1 topology robustness; ED2 divergence times; ED3 depth labels and transitions; ED4 species-level PGLS controls (drawn by ed5());
ED5 preregistered combination sensitivities (drawn by ed4()); ED6 annotation/alignment confound checks (A4); ED7 candidate gene-set
tests (A5), competitive GO for CSUBST qualified pairs and CAFE5 gene-family gains (c, from Fig. 4d); ED8 timing tests and power (A6).
Every plotted value is also written to build/ed_facts.tsv.
"""
import math, os, re, sys
import numpy as np
from figcommon_sty import *
import matplotlib.colorbar as _mcb
_mcb.Colorbar.n_rasterize = 10 ** 9      # r22: colour bars drawn as vector cells, never rasterized


def vmat(ax, M, **kw):
    """r22: heat map as vector cells (pcolormesh) with the imshow geometry (cell centres on integers, row 0 at top)."""
    M = np.ma.masked_invalid(np.asarray(M, float)); nr, nc = M.shape
    q = ax.pcolormesh(np.arange(nc + 1) - .5, np.arange(nr + 1) - .5, M, edgecolors="face", linewidth=0.15, **kw)
    ax.set_xlim(-.5, nc - .5); ax.set_ylim(nr - .5, -.5)
    return q


facts = []
STAR = lambda q: "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def lab(c):
    return CLADE_LABEL.get(c, c)


def ed1():
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 55 * MM), gridspec_kw=dict(wspace=0.45, left=0.07, right=0.985, top=0.88, bottom=0.2))
    ax = axs[0]; panel(ax, "a", -0.25)
    bins = np.linspace(0, 100, 26)
    for f, nm, col in (("gcf", "Constrained", DEEP), ("gcf_ml", "Unconstrained ML", SHAL), ("gcf_unc", "wASTRAL", NULL)):
        v = np.array([fnum(r["gCF"]) for r in rtsv(W("final_v3", f + ".cf.stat"))]); v = v[~np.isnan(v)]
        ax.hist(v, bins, histtype="step", color=col, lw=0.9, label="%s (mean %.2f%%)" % (nm, v.mean()))
        facts.append(("ED1a", f, "n %d median %.2f mean %.3f" % (len(v), np.median(v), v.mean()), "work/final_v3/%s.cf.stat" % f))
    ax.set_xlabel("Gene concordance factor (%)"); ax.set_ylabel("Branches"); ax.set_ylim(0, 75); ax.legend(loc="upper left", fontsize=5)
    ax = axs[1]; panel(ax, "b", -0.25)
    s = np.array([fnum(r["sCF"]) for r in rtsv(W("final_v3", "scfl.cf.stat"))]); s = s[~np.isnan(s)]
    g = {r["ID"]: fnum(r["gCF"]) for r in rtsv(W("final_v3", "gcf.cf.stat"))}
    sc = {r["ID"]: fnum(r["sCF"]) for r in rtsv(W("final_v3", "scfl.cf.stat"))}
    k = [i for i in g if i in sc and not (math.isnan(g[i]) or math.isnan(sc[i]))]
    ax.scatter([g[i] for i in k], [sc[i] for i in k], s=3, color=DEEP, lw=0, alpha=0.7)
    ax.plot([0, 100], [33.3, 33.3], color=NULL, lw=0.5, ls=":"); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel("gCF (%)"); ax.set_ylabel("Site concordance factor, sCF (%)")
    ax.text(98, 4, "n = %d branches; median sCF %.1f%%" % (len(k), np.median(s)), ha="right", fontsize=5)
    facts.append(("ED1b", "scfl", "n %d median sCF %.2f" % (len(k), np.median(s)), "work/final_v3/scfl.cf.stat"))
    ax = axs[2]; panel(ax, "c", -0.32)
    txt = open(W("final_v3", "au.iqtree"), encoding="utf-8").read().split("USER TREES")[1]
    rows = [l.split() for l in txt.splitlines() if re.match(r"^\s+\d+\s+-?\d", l)]
    names = {r["tree"]: r["topology"] for r in rtsv(W("final_v3", "au_trees.tsv"))}
    names = {k: {"fixed": "constrained", "unconstrained wASTRAL": "wASTRAL", "concatenated ML": "unconstrained ML"}.get(v, v) for k, v in names.items()}  # r24 tree names
    for i, r in enumerate(rows):
        dl, pau = fnum(r[2]), fnum(r[-2] if r[-1] in "+-" else r[-1])
        ax.barh(i, dl, 0.6, color=DEEP if dl == 0 else LIGHT)
        # r24: P below 0.001 shown as "P < 0.001" (exact values in the text and Source Data; avoids sub-5-pt exponents)
        ps = "< 0.001" if pau < 1e-3 else "= %.2g" % pau if pau < 1 else "= 1"
        ax.text(dl + 200, i, "ΔlnL %.0f; AU test P %s" % (dl, ps), va="center", fontsize=5)
        facts.append(("ED1c", names.get(r[0], r[0]), "deltaL %s pAU %s" % (r[2], pau), "work/final_v3/au.iqtree"))
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([names.get(r[0], r[0]).replace(" ", "\n", 1) for r in rows])
    ax.set_ylim(len(rows) - 0.5, -0.5); ax.set_xlim(0, 22000); ax.set_xlabel("ΔlnL from best tree (concatenated loci)")
    save(fig, "extended_data", "ED_Fig1")


def ed2():
    na = rtsv(W("mcmctree_bb", "interp", "node_ages.tsv"))
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 72 * MM), gridspec_kw=dict(wspace=0.42, left=0.07, right=0.985, top=0.88, bottom=0.19))
    ax = axs[0]; panel(ax, "a", -0.25)
    for bb, col, nm in (("0", LIGHT, "Interpolated nodes"), ("1", DEEP, "MCMCTree backbone nodes")):
        rr = [r for r in na if r["backbone"] == bb]
        x = np.array([fnum(r["lsd2_Ma"]) for r in rr]); y = np.array([fnum(r["draw_med"]) for r in rr])
        lo = np.array([fnum(r["draw_lo"]) for r in rr]); hi = np.array([fnum(r["draw_hi"]) for r in rr])
        ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o", ms=1.6, color=col, ecolor=col, elinewidth=0.3, lw=0, label="%s (n=%d)" % (nm, len(rr)))
    ax.plot([0, 900], [0, 900], color=NULL, lw=0.5, ls=":"); ax.set_xlim(0, 900); ax.set_ylim(0, 900)
    ax.set_xlabel("LSD2 age (Ma)"); ax.set_ylabel("Backbone-dated age (Ma; median, 95% range)")
    rep = open(W("mcmctree_bb", "interp", "report.txt"), encoding="utf-8").read()
    r = re.search(r"r\(clock-2 posterior mean, LSD2\) ([0-9.]+)", rep).group(1)
    ax.text(880, 40, "backbone r = %s" % r, ha="right", fontsize=5); ax.legend(loc="upper left", fontsize=5, markerscale=1.5)
    facts.append(("ED2a", "r_backbone_lsd2", r, "work/mcmctree_bb/interp/report.txt"))
    ax = axs[1]; panel(ax, "b", -0.25)
    x = np.array([fnum(r["point_correlated_Ma"]) for r in na])
    for k, col, nm in (("point_relaxed_Ma", SHAL, "Independent rates"), ("point_clock_Ma", NULL, "Strict clock"), ("c3_correlated_Ma", INK, "Calibration set 3")):
        y = np.array([fnum(r[k]) for r in na]); ok = ~(np.isnan(x) | np.isnan(y))
        rr = np.corrcoef(x[ok], y[ok])[0, 1]
        ax.scatter(x[ok], y[ok], s=2, color=col, lw=0, label="%s (r = %.3f)" % (nm, rr))
        facts.append(("ED2b", k, "r vs point_correlated %.4f n %d" % (rr, ok.sum()), "work/mcmctree_bb/interp/node_ages.tsv"))
    ax.plot([0, 900], [0, 900], color=NULL, lw=0.5, ls=":"); ax.set_xlim(0, 900); ax.set_ylim(0, 900)
    ax.set_xlabel("Autocorrelated-rates age (Ma)"); ax.set_ylabel("Alternative interpolation (Ma)"); ax.legend(loc="upper left", fontsize=5)
    ax = axs[2]; panel(ax, "c", -0.42)
    cal = [c for c in rtsv(W("final_v3", "date", "calibrations_used.tsv")) if c["applied"] == "yes"]
    for i, c in enumerate(cal):
        ax.plot([fnum(c["min_Ma"]), fnum(c["max_Ma"])], [i, i], color=DEEP, lw=1.6, solid_capstyle="butt")
    ax.set_yticks(range(len(cal))); ax.set_yticklabels([c["node"].replace("CROWN ", "").title()[:24] for c in cal], fontsize=5)
    ax.set_ylim(len(cal) - 0.5, -0.5); ax.set_xlim(0, 900); ax.set_xlabel("Calibration bounds (Ma)")
    facts.append(("ED2c", "n_calibrations_applied", str(len(cal)), "work/final_v3/date/calibrations_used.tsv"))
    save(fig, "extended_data", "ED_Fig2")


def ed3():
    L = rtsv(os.path.join(DATA, "meta", "labels_v2.tsv"))
    tips = set(re.findall(r"[(,]([^(),:;]+)", open(W("final_v3", "species_tree.nwk")).read()))
    L = [r for r in L if r["code"] in tips]
    CL = [lab(c) for c in CLADES]
    grp = lambda r: r["clade"] if r["clade"] in CL else "Other"
    G = CL + ["Other"]
    fig = plt.figure(figsize=(180 * MM, 118 * MM))
    gs = fig.add_gridspec(2, 2, wspace=0.32, hspace=0.55, left=0.08, right=0.985, top=0.94, bottom=0.16)
    ax = fig.add_subplot(gs[0, :]); panel(ax, "a", -0.06)
    rng = np.random.default_rng(1)
    for i, gname in enumerate(G):
        for st, col, dx in (("D", DEEP, 0.17), ("S", SHAL, -0.17)):
            d = np.array([fnum(r["depth_median"]) for r in L if grp(r) == gname and r["lab"] == st]); d = d[~np.isnan(d)]
            ax.scatter(i + dx + rng.uniform(-0.09, 0.09, len(d)), np.log10(np.clip(d, 0.1, None)), s=2.5, color=col, lw=0)
            facts.append(("ED3a", gname, "%s n %d median depth %.0f m" % (st, len(d), np.median(d) if len(d) else float("nan")), "meta/labels_v2.tsv"))
    ax.axhline(np.log10(200), color=NULL, lw=0.5, ls=":"); ax.text(-0.45, np.log10(200) + 0.05, "200 m", ha="left", fontsize=5, color=NULL)
    ax.set_xticks(range(len(G))); ax.set_xticklabels(G); ax.set_ylabel("Median occurrence depth (log10 m)")
    ax.scatter([], [], s=6, color=DEEP, label="Deep (lab = D)"); ax.scatter([], [], s=6, color=SHAL, label="Shallow (lab = S)")
    ax.legend(loc="upper right", ncol=2, fontsize=5)
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "b", -0.14)
    w = 0.2
    for j, (k, nm, col) in enumerate((("lab", "All deep labels", DEEP), ("lab_strict", "Strict", "#5b84ad"), ("lab_1000", "≥1,000 m", "#8fb3d4"))):
        n = [sum(1 for r in L if grp(r) == g and r[k] == "D") for g in G]
        ax.bar(np.arange(len(G)) + (j - 1.5) * w, n, w, color=col, label=nm)
        facts.append(("ED3b", k, " ".join("%s:%d" % (g, x) for g, x in zip(G, n)), "meta/labels_v2.tsv"))
    n = [sum(1 for r in L if grp(r) == g and r["lab"] == "S") for g in G]
    ax.bar(np.arange(len(G)) + 1.5 * w, n, w, color=SHAL, label="Shallow")
    ax.set_xticks(range(len(G))); ax.set_xticklabels(G, rotation=45, ha="right"); ax.set_ylabel("Species")
    ax.legend(loc="upper right", fontsize=5, ncol=2)
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "c", -0.14)
    tr = rtsv(W("final_v3", "asr", "transitions.tsv"))
    for k, col, nm in (("S,D", DEEP, "Shallow → deep"), ("D,S", SHAL, "Deep → shallow")):
        v = np.array([fnum(r[k]) for r in tr])
        ax.hist(v, np.arange(150, 271, 4), color=col, alpha=0.75, lw=0,
                label="%s: mean %.1f (95%% %.0f–%.0f)" % (nm, v.mean(), np.percentile(v, 2.5), np.percentile(v, 97.5)))
        facts.append(("ED3c", k, "mean %.1f 95%% %.0f-%.0f n %d" % (v.mean(), np.percentile(v, 2.5), np.percentile(v, 97.5), len(v)), "work/final_v3/asr/transitions.tsv"))
    ax.set_xlabel("Transitions per stochastic map (ARD Mk, 1,000 maps)"); ax.set_ylabel("Maps"); ax.set_ylim(0, 200); ax.legend(loc="upper left", fontsize=5, frameon=False)
    save(fig, "extended_data", "ED_Fig3")


T25 = ["RELAX p<0.05", "BUSTED-PH FG p<0.05", "BUSTED-PH strict (FG p<0.05 and diff p<0.05)", "CSUBST perm p<0.05",
       "CSUBST qualified pairs", "CSUBST nochemo qualified pairs"]
T25L = ["RELAX", "BUSTED-PH", "BUSTED-PH strict", "CSUBST permutation", "CSUBST qualified", "CSUBST qualified,\nno chemosynthetic"]


def ed4():
    fig = plt.figure(figsize=(180 * MM, 82 * MM))
    gs = fig.add_gridspec(1, 3, wspace=0.6, left=0.155, right=0.95, top=0.86, bottom=0.27, width_ratios=[1.1, 1, 1])
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", -0.62)
    V = [("final_main", "Main (rescued RELAX, matched fish)", DEEP, "o"), ("si_norescue", "No RELAX rescue", SHAL, "s"),
         ("si_fishold", "Original fish design", NULL, "^")]
    for j, (f, nm, col, mk) in enumerate(V):
        cb = {r["test"]: r for r in rtsv(W("combine25c", f + ".combined.tsv"))}
        z = [fnum(cb[t]["stouffer_z"]) for t in T25]
        ax.plot(z, np.arange(len(T25)) + (j - 1) * 0.2, mk, ms=3, color=col, lw=0, label=nm)
        for t in T25:
            facts.append(("ED5a", f, t, "z %s p %s holm %s" % (cb[t]["stouffer_z"], cb[t]["p_stouffer"], cb[t]["p_holm"])))
    der = rtsv(os.path.join(BUILD, "src_25c_derived_sensitivity.tsv"))
    for i, t in enumerate(T25):
        lo = [fnum(r["z"]) for r in der if r["test"] == t and r["variant"] == "LOCO"]
        if lo:
            ax.plot([min(lo), max(lo)], [i + 0.3, i + 0.3], color=INK, lw=0.8)
            facts.append(("ED5a", "LOCO", t, "z range %.3f-%.3f" % (min(lo), max(lo))))
    ax.plot([], [], color=INK, lw=0.8, label="Leave-one-clade-out range (main)")
    ax.axvline(1.645, color=NULL, lw=0.5, ls=":"); ax.set_yticks(range(len(T25))); ax.set_yticklabels(T25L, fontsize=5.5)
    ax.set_ylim(len(T25) - 0.4, -0.6); ax.set_xlabel("Weighted Stouffer z (w = √n)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.35, -0.2), fontsize=5, ncol=2, columnspacing=0.6, handletextpad=0.2, frameon=False)
    for k, (f, letter, nm) in enumerate((("final_main", "b", "Main"), ("si_norescue", "c", "No RELAX rescue"))):
        ax = fig.add_subplot(gs[0, k + 1]); panel(ax, letter, -0.42)
        pc = rtsv(W("combine25c", f + ".perclade.tsv"))
        M = np.full((len(CLADES), len(T25)), np.nan); S = [[""] * len(T25) for _ in CLADES]
        for r in pc:
            c = r["clade"] if r["clade"] in CLADES else ("Actinopteri_m" if r["clade"].startswith("Actinopteri") else None)
            if c is None or r["test"] not in T25:
                continue
            i, j = CLADES.index(c), T25.index(r["test"])
            M[i, j] = -math.log10(max(fnum(r["paired_perm_p"]), 1e-4)); S[i][j] = STAR(fnum(r["p_holm_within_test"]))
            facts.append(("ED5" + letter, c, r["test"], "p %s holm %s n %s" % (r["paired_perm_p"], r["p_holm_within_test"], r["n_paired"])))
        im = vmat(ax, M, cmap="Blues", vmin=0, vmax=4)
        for i in range(len(CLADES)):
            for j in range(len(T25)):
                if np.isnan(M[i, j]):
                    ax.text(j, i, "–", ha="center", va="center", fontsize=5, color=NULL)
                elif S[i][j]:
                    ax.text(j, i, S[i][j], ha="center", va="center", fontsize=5, color="white" if M[i, j] > 2.5 else INK)
        ax.set_xticks(range(len(T25))); ax.set_xticklabels(["RELAX", "BPH", "BPH strict", "CS perm.", "CS qual.", "CS qual. no chemo."], rotation=55, ha="right", fontsize=5)
        ax.set_yticks(range(len(CLADES))); ax.set_yticklabels([lab(c) for c in CLADES], fontsize=5.5); ax.set_title(nm, fontsize=6)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=fig.axes[-2:], fraction=0.03, pad=0.02); cb.set_label("−log10 paired permutation P", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    save(fig, "extended_data", "ED_Fig5")


PG_COLS = [("final_v3/pgls_controls", "all", "Deep label"), ("final_v3/pgls_controls", "ref", "Reference annotations"),
           ("final_v3/pgls_controls", "deepT", "+ habitat T"), ("final_v3/pgls_controls", "deepTGC", "+ T + GC"),
           ("final_v3/pgls_controls", "nochemo", "No chemosynthetic"), ("final_v3/pgls_controls", "ldepthT", "log depth + T"),
           ("final_v3/pgls_controls", "deepdose", "log depth, deep only"), ("final_v3/pgls_controls_strict", "all", "Strict labels"),
           ("final_v3/pgls_controls_1000", "all", "≥1,000 m labels"), ("final_v3_ml/pgls_controls", "all", "Unconstrained ML tree"),
           ("final_v3_unc/pgls_controls", "all", "wASTRAL tree")]


def ed5():
    tabs = {d: {r["trait"]: r for r in rtsv(W(*d.split("/"), "controls_summary.tsv"))} for d in set(c[0] for c in PG_COLS)}
    main = tabs["final_v3/pgls_controls"]
    tz = lambda r, m: fnum(r.get(m + "_est")) / fnum(r.get(m + "_se")) if r else float("nan")
    traits = sorted(main, key=lambda t: -abs(tz(main[t], "all")))
    M = np.array([[tz(tabs[d].get(t), m) for d, m, _ in PG_COLS] for t in traits])
    Q = [[fnum(tabs[d].get(t, {}).get(m + "_q")) for d, m, _ in PG_COLS] for t in traits]
    fig, ax = plt.subplots(figsize=(120 * MM, 150 * MM)); fig.subplots_adjust(left=0.3, right=0.88, top=0.97, bottom=0.17)
    im = vmat(ax, M, cmap="RdBu_r", vmin=-6, vmax=6)
    for i, t in enumerate(traits):
        for j, (d, m, nm) in enumerate(PG_COLS):
            s = STAR(Q[i][j]) if not math.isnan(Q[i][j]) else ""
            if s:
                ax.text(j, i, s, ha="center", va="center", fontsize=5, color="white" if abs(M[i, j]) > 4 else INK)
            facts.append(("ED4", t, "%s/%s" % (d, m), "t %.3f q %s" % (M[i, j], Q[i][j])))
    AA3 = dict(A="Ala", C="Cys", D="Asp", E="Glu", F="Phe", G="Gly", H="His", I="Ile", K="Lys", L="Leu", M="Met", N="Asn",
               P="Pro", Q="Gln", R="Arg", S="Ser", T="Thr", V="Val", W="Trp", Y="Tyr")
    PROP = {"pI": "Isoelectric point", "R_RK": "Arg/(Arg+Lys)", "acid_base": "Acidic/basic ratio", "net_charge7_per_res": "Net charge (pH 7)",
            "charged_DEKR": "Charged (DEKR)", "ENC": "Effective no. codons", "gravy": "GRAVY", "ST": "Ser+Thr", "IVYWREL": "IVYWREL",
            "PG": "Pro+Gly", "aliphatic": "Aliphatic index", "aromatic": "Aromatic", "GC3": "GC3", "GC4": "GC4 (fourfold)", "cds_GC": "CDS GC"}

    def nice(t):
        scope, prop = ("BUSCO", t[6:]) if t.startswith("busco_") else ("proteome", t[6:])
        lab = AA3[prop[2:]] + " frequency" if prop.startswith("f_") else PROP.get(prop, prop)
        return "%s (%s)" % (lab, scope)
    ax.set_yticks(range(len(traits))); ax.set_yticklabels([nice(t) for t in traits], fontsize=5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if math.isnan(M[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color=LIGHT, lw=0))
                ax.text(j, i, "n.e.", ha="center", va="center", fontsize=5, color=INK)
    ax.set_xticks(range(len(PG_COLS))); ax.set_xticklabels([c[2] for c in PG_COLS], rotation=55, ha="right", fontsize=5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("PGLS t (deep − shallow; positive = higher in deep)", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    save(fig, "extended_data", "ED_Fig4")


A4I = [("relax", "RELAX any"), ("relax_relaxed", "RELAX relaxed"), ("relax_intens", "RELAX intensified"), ("bph", "BUSTED-PH"),
       ("bph_strict", "BUSTED-PH strict"), ("cs_perm", "CSUBST perm."), ("cs_qual", "CSUBST qualified")]


def ed6():
    pc = [r for r in rtsv(W("addon", "A4", "perclade.tsv")) if r["term"] == "is_real"]
    nf = rtsv(W("addon", "A4", "null_fpr.tsv"))
    fig, axs = plt.subplots(1, 2, figsize=(180 * MM, 75 * MM), gridspec_kw=dict(wspace=0.3, left=0.1, right=0.985, top=0.86, bottom=0.14))
    cols = [DEEP, "#5b84ad", "#8fb3d4", SHAL, "#a07020", NULL, INK]
    for ax, letter, title in ((axs[0], "a", "Clade-level M1 (HOG-clustered s.e.)"), (axs[1], "b", "Real rate vs matched-null 95% interval")):
        panel(ax, letter, -0.17); ax.set_title(title, fontsize=6)
    ax = axs[0]
    for j, (ind, nm) in enumerate(A4I):
        for i, c in enumerate(CLADES):
            r = [x for x in pc if x["ind"] == ind and x["clade"] == c]
            if not r:
                continue
            e, s = fnum(r[0]["est"]), fnum(r[0]["se_cl"]); y = i + (j - 3) * 0.11
            if abs(e) > 8 or math.isnan(s):
                continue
            ax.plot([e - 1.96 * s, e + 1.96 * s], [y, y], color=cols[j], lw=0.5); ax.plot(e, y, "o", ms=1.8, color=cols[j])
            facts.append(("ED6a", c, ind, "logOR %.3f se_cl %.3f p_cl %s" % (e, s, r[0]["p_cl"])))
    ax.axvline(0, color=NULL, lw=0.5, ls=":"); ax.set_xlim(-4, 4)
    ax.set_yticks(range(len(CLADES))); ax.set_yticklabels([lab(c) for c in CLADES]); ax.set_ylim(len(CLADES) - 0.5, -0.6)
    ax.set_xlabel("log odds ratio, real deep vs matched shallow foreground")
    ax = axs[1]
    for j, (ind, nm) in enumerate(A4I):
        for i, c in enumerate(CLADES):
            r = [x for x in nf if x["ind"] == ind and x["clade"] == c]
            if not r:
                continue
            y = i + (j - 3) * 0.11; lo, hi, rr = fnum(r[0]["lo"]), fnum(r[0]["hi"]), fnum(r[0]["real_rate"])
            ax.plot([lo, hi], [y, y], color=LIGHT, lw=0.9); ax.plot(rr, y, "o", ms=1.8, color=cols[j])
            facts.append(("ED6b", c, ind, "real %.3f null %s (%.3f-%.3f)" % (rr, r[0]["rate"], lo, hi)))
        ax.plot([], [], "o", ms=2.5, color=cols[j], label=nm)
    ax.set_yticks(range(len(CLADES))); ax.set_yticklabels([lab(c) for c in CLADES]); ax.set_ylim(len(CLADES) - 0.5, -0.6)
    ax.set_xlabel("Fraction of orthogroups called"); ax.set_xscale("symlog", linthresh=0.05); ax.set_xlim(0, 0.8)
    ax.set_xticks([0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8]); ax.set_xticklabels(["0", "0.02", "0.05", "0.1", "0.2", "0.4", "0.8"]); ax.minorticks_off()
    fig.legend(*axs[1].get_legend_handles_labels(), loc="upper center", ncol=7, fontsize=5, bbox_to_anchor=(0.55, 1.0), handletextpad=0.1)
    save(fig, "extended_data", "ED_Fig6")


def ed7():
    a5 = rtsv(W("addon", "A5", "combined.tsv"))
    sets = []
    for r in a5:
        if (r["set"], r["label"]) not in sets:
            sets.append((r["set"], r["label"]))
    TT = [("relax_relaxed", "RELAX relaxed"), ("relax_intens", "RELAX intensified"), ("bph_strict", "BUSTED-PH strict"), ("cs_qual", "CSUBST qualified")]
    fig = plt.figure(figsize=(180 * MM, 135 * MM))
    gs = fig.add_gridspec(2, 2, wspace=1.05, hspace=0.55, left=0.23, right=0.985, top=0.94, bottom=0.08, width_ratios=[0.8, 1], height_ratios=[1, 0.62])
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", -0.95)
    M = np.full((len(sets), len(TT)), np.nan)
    for r in a5:
        if r["test"] in dict(TT):
            i, j = sets.index((r["set"], r["label"])), [t[0] for t in TT].index(r["test"]); M[i, j] = fnum(r["Z"])
            if math.isnan(M[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color=LIGHT, lw=0, zorder=3))
                ax.text(j, i, "nt", ha="center", va="center", fontsize=5, color=INK, zorder=4)
            elif r["call"]:
                ax.text(j, i, {"n.s.": "", "pass": "+", "FAIL (departs from null)": "×", "not testable": "nt"}.get(r["call"], r["call"]),
                        ha="center", va="center", fontsize=5.5, color="white" if abs(M[i, j]) > 2 else INK, zorder=4)
            facts.append(("ED7a", r["set"], r["test"], "Z %s p_adj %s call %s n_clades %s" % (r["Z"], r["p_adj"], r["call"], r["n_clades"])))
    im = vmat(ax, M, cmap="RdBu_r", vmin=-3, vmax=3)
    ax.set_yticks(range(len(sets))); ax.set_yticklabels(["%s %s" % s for s in sets], fontsize=5)
    ax.set_xticks(range(len(TT))); ax.set_xticklabels([t[1] for t in TT], rotation=50, ha="right", fontsize=5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03); cb.set_label("Stouffer z (deep vs matched null)", fontsize=5); cb.ax.tick_params(labelsize=5)
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", -0.75)
    cq = sorted([r for r in rtsv(W("addon", "D3_func_competitive", "competitive_cs_qual.tsv")) if fnum(r["q"]) < 0.05],
                key=lambda r: fnum(r["p"]))
    for i, r in enumerate(cq):
        q = fnum(r["q"])
        ax.barh(i, -math.log10(fnum(r["p"])), 0.7, color=DEEP if q < 0.05 else LIGHT)
        facts.append(("ED7b", r["term"], r["name"], "Z %s p %s q %s clades %s" % (r["Z"], r["p"], r["q"], r["n_clades"])))
    ax.set_yticks(range(len(cq))); ax.set_yticklabels([((r["name"][:40] + "…") if len(r["name"]) > 41 else r["name"]) + " (%s)" % r["namespace"].replace("GO_", "") for r in cq], fontsize=5)
    ax.set_ylim(len(cq) - 0.5, -0.5); ax.set_xlabel("−log10 P, competitive GO test\n(CSUBST qualified pairs)")
    ax.text(0.98, 0.04, "All %d of 3,289 terms with q < 0.05\nCSUBST qualified" % len(cq) + " is not robust to the\nannotation confound model (ED Fig. 7)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5, color=INK)
    # c: CAFE5 (moved from Fig. 4d on 2026-10-30; failed fits are indeterminate)
    ax = fig.add_subplot(gs[1, :]); panel(ax, "c", -0.27)
    # 2026-10-29: 40 tip permutations per clade (R2, PREREG_nature.md; work/addon/R2/cafe_perm/cafe_gains_long.tsv,
    # cafe_rank40.tsv); Mollusca and Actinopteri from the family-range rescue; Porifera has its original 10 permutations.
    OR = ["Mollusca", "Cnidaria", "Echinodermata", "Annelida", "Arthropoda", "Actinopteri_m", "Chondrichthyes", "Porifera"]
    LONG = rtsv(W("addon", "R2", "cafe_perm", "cafe_gains_long.tsv"))
    RK = {r["clade"]: r for r in rtsv(W("addon", "R2", "cafe_perm", "cafe_rank40.tsv"))}
    for i, c in enumerate(OR):
        key = c + "_rescue" if any(r["clade"] == c + "_rescue" for r in LONG) else c
        rows = {r["model"]: fnum(r["gain"]) for r in LONG if r["clade"] == key and not math.isnan(fnum(r["gain"]))}
        dg = rows.get("deep", float("nan"))
        mat = [v for m, v in rows.items() if re.fullmatch(r"null\d", m)]; per = [v for m, v in rows.items() if m.startswith("nullp")]
        ax.scatter(per, [i + 0.17] * len(per), s=4, marker="x", color=NULL, lw=0.4, alpha=0.8)
        ax.scatter(mat, [i - 0.17] * len(mat), s=5, color=SHAL, lw=0)
        ax.scatter([dg], [i], s=14, color=DEEP, lw=0, zorder=3)
        rp = (1 + sum(v >= dg for v in per)) / (1 + len(per))
        if key in RK:
            rp = fnum(RK[key]["rank_p"])
        lbl = "P = %.3f" % rp + (" (fit failed; indeterminate)" if dg < 0 else "") + ("" if len(per) == 40 else " (%d permutations)" % len(per))
        if not per:
            rp, lbl = float("nan"), "indeterminate (no fit converged)"
        ax.text(0.995, i, lbl, transform=ax.get_yaxis_transform(), va="center", ha="right", fontsize=5, color=DEEP if rp < 0.05 else INK)
        facts.append(("ED7c", c, "deep %.0f perm rank p %.4f K %d matched %d" % (dg, rp, len(per), len(mat)), key))
    ax.set_xscale("symlog", linthresh=100); ax.axvline(0, color=NULL, lw=0.5, ls=":")
    ax.set_xlim(-1.2e5, 1e6); ax.set_xticks([-100000, -10000, -1000, -100, 0, 100, 1000, 10000]); ax.set_xticklabels(["−10⁵", "−10⁴", "−10³", "−10²", "0", "10²", "10³", "10⁴"])
    ax.minorticks_off()
    ax.set_yticks(range(len(OR))); ax.set_yticklabels([CLADE_LABEL.get(c, c) for c in OR]); ax.set_ylim(len(OR) - 0.5, -0.5)
    ax.set_xlabel("CAFE5 log-likelihood gain of a separate foreground rate")
    ax.scatter([], [], s=14, color=DEEP, label="Deep foreground"); ax.scatter([], [], s=5, color=SHAL, label="Matched shallow (5)")
    ax.scatter([], [], s=4, marker="x", color=NULL, label="Permuted foreground (40)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=5, handletextpad=0.1, borderaxespad=0.2, columnspacing=1)
    save(fig, "extended_data", "ED_Fig7")


def ed8():
    d = rtsv(W("addon", "A6", "decision.tsv")); dl = {r["event"]: r for r in rtsv(W("addon", "A6_lsd", "decision.tsv"))}
    ps = rtsv(W("addon", "A6", "power_summary.tsv"))
    fig, axs = plt.subplots(1, 2, figsize=(180 * MM, 62 * MM), gridspec_kw=dict(wspace=0.5, left=0.17, right=0.80, top=0.88, bottom=0.17))
    ax = axs[0]; panel(ax, "a", -0.32)
    for i, r in enumerate(d):
        ax.plot(fnum(r["z"]), i - 0.12, "o", ms=2.6, color=DEEP); ax.plot(fnum(dl[r["event"]]["z"]), i + 0.12, "s", ms=2.2, color=SHAL)
        facts.append(("ED8a", r["event"], r["family"], "z %s holm %s call %s; lsd z %s" % (r["z"], r["p_holm"], r["call"], dl[r["event"]]["z"])))
    ax.plot([], [], "o", ms=2.6, color=DEEP, label="MCMCTree-backbone ages"); ax.plot([], [], "s", ms=2.2, color=SHAL, label="LSD2 ages")
    ax.axvline(0, color=NULL, lw=0.5, ls=":"); ax.axvline(1.645, color=NULL, lw=0.5, ls="--")
    ax.set_yticks(range(len(d))); ax.set_yticklabels(["%s (%s Ma)" % (r["event"], r["E"]) for r in d], fontsize=5)
    ax.set_ylim(len(d) - 0.5, -0.5); ax.set_xlabel("z, excess deep-colonisation rate in interval"); ax.legend(loc="lower left", fontsize=5, frameon=False); ax.text(1.62, len(d) - 0.7, "one-sided\nP = 0.05", ha="right", va="bottom", fontsize=5, color=NULL)
    ax = axs[1]; panel(ax, "b", -0.2)
    for fam, mk, ls in (("window", "o", "-"), ("epoch", "s", "--")):
        for k, E in enumerate(sorted(set(fnum(r["E"]) for r in ps if r["family"] == fam))):
            rr = sorted([r for r in ps if r["family"] == fam and fnum(r["E"]) == E], key=lambda r: fnum(r["mult"]))
            col = plt.cm.viridis(k / 8 * 0.85) if fam == "window" else ("#a07020", "#d4a03a", "#888888")[k]
            ax.plot([fnum(r["mult"]) for r in rr], [fnum(r["power"]) for r in rr], ls, marker=mk, ms=2, lw=0.6, color=col,
                    label="%s %s Ma" % ("10-Myr window at" if fam == "window" else "after", ("%g" % E)))
            for r in rr:
                facts.append(("ED8b", fam, "E %s mult %s" % (r["E"], r["mult"]), "power %s n %s" % (r["power"], r["n"])))
    ax.axhline(0.8, color=NULL, lw=0.5, ls=":"); ax.set_xscale("log"); ax.set_xticks([2, 4, 5, 8, 10, 30]); ax.set_xticklabels(["2", "4", "5", "8", "10", "30"])
    ax.minorticks_off(); ax.set_ylim(0, 1); ax.set_xlabel("Simulated rate multiplier in interval"); ax.set_ylabel("Power (50 simulations)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=5, ncol=1, handlelength=2, frameon=False)
    save(fig, "extended_data", "ED_Fig8")


if __name__ == "__main__":
    which = sys.argv[1:] or [str(i) for i in range(1, 9)]
    for w in which:
        globals()["ed" + w]()
    p = os.path.join(BUILD, "ed_facts.tsv")
    old = [l for l in (open(p, encoding="utf-8").read().splitlines()[1:] if os.path.exists(p) else []) if not any(l.startswith("ED%s" % w) for w in which)]
    with open(p, "w", encoding="utf-8") as f:
        f.write("panel\tkey\titem\tvalue\n")
        for l in old + ["\t".join(str(x) for x in r) for r in facts]:
            f.write(l + "\n")
