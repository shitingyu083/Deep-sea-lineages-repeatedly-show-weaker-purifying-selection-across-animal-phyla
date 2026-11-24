"""Fig. 4 v2 - what relaxes (R6-3a omega-class decomposition) and whether it recurs across independent deep-sea
origins (R6-2 per-origin Group-mode RELAX, ADDENDUM_R6_2 0e9d4fe6).

Forms: a, b, d, e = Sackton et al. 2019 Fig. 4A (dots, grey capped whiskers, labels at 45 degrees) with the dashed
reference line of Sackton Fig. 3C and the dark condition box of Cell 2025 Fig. 4G; c = terr Extended Data Fig. 8a
(binned distribution, white bin edges, horizontal grid, full-height red line at the observed value with red text),
two stacked panels on one x axis (deep lineages; matched null pseudo-lineages in grey).
Data (read-only): R6/R6_3/r6_3_stats.tsv ca2600c9 (a, b); R6/R6_2/r6_2_summary.tsv acb54400 (d, e: confirmatory
numbers as written by the frozen analysis); R6/R6_2/r6_2_lineages.tsv 997e9842 (c top, e lineage counts);
R6/R6_2/r6_2_fits.tsv 97078a7f (c bottom: null pseudo-lineage effects by the ADDENDUM lineage rule, descriptive only;
the same code reproduces all 65 deep-lineage effects of r6_2_lineages.tsv, asserted below). No bootstrap is redone.
Outputs in viz_v2_sty/out only.
"""
from figcommon_v2_sty import *  # noqa: F401,F403
from collections import defaultdict

mt_arial()
S3 = {(r["analysis"], r["metric"], r["clade"]): r for r in rtsv(W("addon", "R6", "R6_3", "r6_3_stats.tsv"))}
g3 = lambda m, c: tuple(fnum(S3[("R6-3a_paired", m, c)][k]) for k in ("estimate", "lo", "hi"))
assert abs(g3("Dpur", "ALL")[0] - .077440575607305) < 1e-12 and S3[("R6-3a_paired", "Dpur", "ALL")]["n_clades_pos"] == "6"
SM = {r["stat"]: tuple(fnum(r[k]) for k in ("est", "lo", "hi")) for r in rtsv(W("addon", "R6", "R6_2", "r6_2_summary.tsv"))}
assert abs(SM["pooled"][0] + .0512255539301645) < 1e-12 and abs(SM["share_lt0"][0] - .707692307692308) < 1e-12
LIN = rtsv(W("addon", "R6", "R6_2", "r6_2_lineages.tsv"))
assert len(LIN) == 65
# null pseudo-lineage effects (ADDENDUM: logK clamped, lineage median over HOGs where present, >= 5 HOGs, minus the
# clade null median); merged groups excluded
lk = lambda K: np.log10(min(max(K, .001), 50))
val = defaultdict(list)
for r in rtsv(W("addon", "R6", "R6_2", "r6_2_fits.tsv")):
    if r["fit"] == "ok" and r["canon"] not in ("merged", "", "NA") and np.isfinite(fnum(r["K"])):
        val[(r["clade"], r["set"], r["canon"])].append(lk(fnum(r["K"])))
NMED = {c: np.median([x for (cc, s, _), v in val.items() if cc == c and s == "null1" for x in v])
        for c in {k[0] for k in val}}
EFF = {(c, s, cn): (np.median(v) - NMED[c], len(v)) for (c, s, cn), v in val.items() if len(v) >= 5}
for r in LIN:
    e, n = EFF[(r["clade"], "real", r["canon"])]
    assert abs(e - fnum(r["eff"])) < 1e-12 and n == int(r["n_hog"]), ascii(r)
assert sum(1 for k in EFF if k[1] == "real") == 65
REAL = np.array([fnum(r["eff"]) for r in LIN])
NULLE = np.array([e for (c, s, _), (e, _) in EFF.items() if s == "null1"])
assert len(NULLE) == 83 and abs((REAL < 0).mean() - SM["share_lt0"][0]) < 1e-12
POOLC = ["Porifera"] + [c for c in PHYLO if c != "Porifera"]
NL = {c: (sum(1 for r in LIN if r["clade"] == c and fnum(r["eff"]) < 0), sum(1 for r in LIN if r["clade"] == c))
      for c in PHYLO}
assert NL["Porifera"][1] < 2 and all(NL[c][1] >= 2 for c in PHYLO if c != "Porifera"), NL

Wf, Hf = 180.0, 96.0
T1, H1, T2, H2 = 5.0, 26.0, 50.0, 31.0
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
f6 = lambda v: "%.6f" % v
src = [("panel", "group", "series", "estimate", "lo95", "hi95", "n", "note")]
facts = [("key", "value", "source")]
LAB = [clab(c) for c in PHYLO] + [ALL_LAB]
xs = np.r_[np.arange(8), 8.8]
for k, (m, nm, x0, lim, tk) in enumerate((("Dpur", "Purifying classes ($\\Delta_{\\mathrm{pur}}$)", 11.0, (-.3, .45), .15),
                                          ("Dpos", "Positive classes ($\\Delta_{\\mathrm{pos}}$)", 98.0, (-.06, .06), .03))):
    ax = axes_mm(fig, Wf, Hf, x0, T1, 70.0, H1)
    v = np.array([g3(m, c) for c in PHYLO] + [g3(m, "ALL")])
    assert v[:, 1].min() > lim[0] and v[:, 2].max() < lim[1] - .35 * (lim[1] - lim[0]) / 2, ascii((m, v[:, 1].min(), v[:, 2].max()))
    ref_line(ax, 0)
    dot_whisker(ax, xs, v[:, 0], v[:, 1], v[:, 2], LAB, lab_fs=5)
    ax.set_ylim(*lim); ax.set_yticks(np.round(np.arange(lim[0], lim[1] + tk / 2, tk), 3))
    ax.set_ylabel("Excess over matched null", fontsize=5.5)
    cond_box(ax, nm, fs=5)
    letter(fig, Wf, Hf, "ab"[k], 1.0 if k == 0 else 88.0, T1 - 1.9)
    r = S3[("R6-3a_paired", m, "ALL")]
    for c, row in zip(PHYLO + ["ALL"], v):
        src.append(("ab"[k], c, m, f6(row[0]), f6(row[1]), f6(row[2]), S3[("R6-3a_paired", m, c)]["n"],
                    "r6_3_stats.tsv; p_boot %s" % S3[("R6-3a_paired", m, c)]["p_boot"]))
    facts.append(("R6_3a_%s" % m, "%.3f (%.3f to %.3f), p_boot %s, %s of 8 clades > 0" % (
        tuple(v[-1]) + (r["p_boot"], r["n_clades_pos"])), "r6_3_stats.tsv"))
# c: deep-lineage effects (top) and null pseudo-lineage effects (bottom), one x axis
bins = np.arange(-1.1, 1.11, .1)
assert REAL.min() > bins[0] and REAL.max() < bins[-1] and NULLE.min() > bins[0] and NULLE.max() < bins[-1], \
    (REAL.min(), REAL.max(), NULLE.min(), NULLE.max())
hh = (H2 - 3.0) / 2
for j, (v, col, at, lines) in enumerate((
        (REAL, None, SM["pooled"][0], ["Pooled %s" % ("%.3f" % SM["pooled"][0]).replace("-", "−"), "%.0f%% < 0" % (100 * SM["share_lt0"][0])]),
        (NULLE, LIGHT, 0.0, ["Null", "%.0f%% < 0" % (100 * (NULLE < 0).mean())]))):
    ax = axes_mm(fig, Wf, Hf, 11.0, T2 + j * (hh + 3.0), 46.0, hh)
    ed8_hist(ax, v, bins, at, lines, xlabel="Lineage effect, Δ median log$_{10}$ $K$" if j else None,
             ylabel="Lineages", color=col, head=1.6)
    ax.set_xlim(bins[0], bins[-1])
    if not j:
        ax.set_xticklabels([])
    cond_box(ax, ("Deep lineages (%d)" % len(REAL), "Null pseudo-lineages (%d)" % len(NULLE))[j], fs=5)
letter(fig, Wf, Hf, "c", 1.0, T2 - 1.9)
for r in LIN:
    src.append(("c", "%s|%s" % (r["clade"], r["canon"]), "deep_lineage", f6(fnum(r["eff"])), "", "", r["n_hog"],
                "r6_2_lineages.tsv"))
for (c, s, cn), (e, n) in sorted(EFF.items()):
    if s == "null1":
        src.append(("c", "%s|%s" % (c, cn), "null_pseudo_lineage", f6(e), "", "", n, "r6_2_fits.tsv, descriptive"))
facts.append(("lineage_effects", "deep %d lineages, %.1f%% < 0 (%.3f to %.3f); null %d pseudo-lineages, %.1f%% < 0" % (
    len(REAL), 100 * (REAL < 0).mean(), REAL.min(), REAL.max(), len(NULLE), 100 * (NULLE < 0).mean()),
    "r6_2_lineages.tsv; r6_2_fits.tsv [null new, descriptive]"))
# d: pooled effect and the pre-declared leave-outs
DS = [("pooled", "Pooled"), ("lineage_weighted", "Lineage weighted"), ("without_top2_clades", "−Top 2 clades"),
      ("without_largest_lineage", "−Largest lineage")] + [("loco_" + c, "−" + clab(c)) for c in PHYLO
                                                                 if "loco_" + c in SM]
ax = axes_mm(fig, Wf, Hf, 73.0, T2, 50.0, H2)
v = np.array([SM[k] for k, _ in DS])
xd = np.r_[np.arange(4), 4.8 + np.arange(len(DS) - 4)]
ref_line(ax, 0)
dot_whisker(ax, xd, v[:, 0], v[:, 1], v[:, 2], [l for _, l in DS], lab_fs=5)
assert v[:, 1].min() > -.2 and v[:, 2].max() < .05
ax.set_ylim(-.2, .05); ax.set_yticks([-.2, -.15, -.1, -.05, 0, .05])
ax.set_ylabel("Pooled lineage effect", fontsize=5.5)
letter(fig, Wf, Hf, "d", 63.0, T2 - 1.9)
for (k, _), row in zip(DS, v):
    src.append(("d", k, "R6-2", f6(row[0]), f6(row[1]), f6(row[2]), "", "r6_2_summary.tsv"))
facts.append(("R6_2_claim", "pooled %.4f (%.4f to %.4f); share<0 %.3f (%.3f to %.3f); without top2 %.4f (%.4f to %.4f)"
              "; SUPPORTED" % (SM["pooled"] + SM["share_lt0"] + SM["without_top2_clades"]), "r6_2_report.txt"))
# e: per-clade mean lineage effect, lineages < 0 / summarised lineages in the label
ax = axes_mm(fig, Wf, Hf, 134.0, T2, 44.0, H2)
EC = [c for c in PHYLO if "clade_" + c in SM]
v = np.array([SM["clade_" + c] for c in EC])
ref_line(ax, 0)
dot_whisker(ax, np.arange(len(EC)), v[:, 0], v[:, 1], v[:, 2], ["%s (%d/%d)" % ((clab(c),) + NL[c]) for c in EC],
            lab_fs=5)
assert v[:, 1].min() > -.7 and v[:, 2].max() < .3
ax.set_ylim(-.7, .3); ax.set_yticks([-.6, -.4, -.2, 0, .2])
ax.set_ylabel("Clade mean lineage effect", fontsize=5.5)
letter(fig, Wf, Hf, "e", 124.0, T2 - 1.9)
for c, row in zip(EC, v):
    src.append(("e", c, "clade_mean", f6(row[0]), f6(row[1]), f6(row[2]), "%d/%d" % NL[c],
                "r6_2_summary.tsv; lineages <0 / summarised"))
    facts.append(("R6_2_clade_%s" % c, "%.4f (%.4f to %.4f); %d of %d lineages < 0" % (tuple(row) + NL[c]),
                  "r6_2_summary.tsv; r6_2_lineages.tsv"))
res = qa_v2(fig)
clean = qa_report(res, "Fig4")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig4_qa_v2.txt"))
write_facts("fig4_source_v2.tsv", src)
write_facts("fig4_facts_v2.tsv", facts)
save_v2(fig, "Fig4")
print("Fig4 written; QA clean:", clean)
