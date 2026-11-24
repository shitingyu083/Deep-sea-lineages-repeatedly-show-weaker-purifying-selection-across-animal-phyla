"""Fig. 3 v2 - the relaxation contrast under one uniform annotation route (R6-1) and the functional-concentration test
(R6-3b).

Forms: a, b, c, f = Sackton et al. 2019 Fig. 4A (dots, grey capped whiskers, labels at 45 degrees) with the dashed
reference line of Sackton Fig. 3C and the dark condition box of Cell 2025 Fig. 4G; d = Cell 2025 Fig. 4A (framed
scatter, two series with OLS lines and 95% bands, statistics in the legend at top left); e = Cell 2025 Fig. 4G/H
(adjacent bars, capped error bars, square-patch legend outside right).
Data (read-only): R6/r6_1_stats.tsv f4760668 (all estimates and CIs of a-c, as written by the frozen R6-1 analysis),
R6/r6_1_fits.tsv 722b737e (d, e; analysis set, kappa and diagnostics re-derived here and asserted against
r6_1_report.txt), R6/R6_3/r6_3_stats.tsv ca2600c9 (f). K is clamped to [0.001, 50] before log10 (PREREG_r6).
e error bars are Wilson 95% intervals (descriptive; new). Outputs in viz_v2_sty/out only.
"""
from figcommon_v2_sty import *  # noqa: F401,F403

mt_arial()
ST = {(r["outcome"], r["estimand"], r["clade"]): r for r in rtsv(W("addon", "R6", "r6_1_stats.tsv"))}
g3 = lambda o, e, c="ALL": tuple(fnum(ST[(o, e, c)][k]) for k in ("estimate", "lo", "hi"))
assert abs(g3("relaxed", "U")[0] - 1.43786515203715) < 1e-12 and abs(g3("relaxed", "R")[0] - .945825951330721) < 1e-12
rows = rtsv(W("addon", "R6", "r6_1_fits.tsv"))
assert len(rows) == 5886, len(rows)
F = {(r["clade"], r["hog"], r["set"], r["arm"]): r for r in rows}
HOGS = sorted({(r["clade"], r["hog"]) for r in rows})
ok = lambda r: r is not None and np.isfinite(fnum(r["relax_p"])) and np.isfinite(fnum(r["relax_K"]))
A = [h for h in HOGS if all(ok(F.get(h + (s, a))) for s in ("real", "null1") for a in ("O", "U"))]
assert len(HOGS) == 981 and len(A) == 844, (len(HOGS), len(A))
lk = lambda r: np.log10(min(max(fnum(r["relax_K"]), .001), 50))
rel = lambda r: int(fnum(r["relax_p"]) < .05 and fnum(r["relax_K"]) < 1)
LK, KAP, R_ = {}, {}, {}
for s in ("real", "null1"):
    LK[s] = np.array([[lk(F[h + (s, a)]) for a in ("O", "U")] for h in A])
    y = np.array([[rel(F[h + (s, a)]) for a in ("O", "U")] for h in A])
    po = (y[:, 0] == y[:, 1]).mean(); p1, p2 = y.mean(0)
    KAP[s] = (po - (p1 * p2 + (1 - p1) * (1 - p2))) / (1 - (p1 * p2 + (1 - p1) * (1 - p2)))
    R_[s] = np.corrcoef(LK[s][:, 0], LK[s][:, 1])[0, 1]
assert abs(KAP["real"] - .8350809) < 1e-7 and abs(KAP["null1"] - .8103821) < 1e-7, KAP   # r6_1_report.txt
DMED = {s: float(np.median(LK[s][:, 1] - LK[s][:, 0])) for s in LK}
DIAG = {}   # (set, arm) -> (fail events, K-at-bound events, n built)
REP = {("null1", "O"): (.0540, .0071), ("real", "O"): (.0754, .0204), ("null1", "U"): (.0510, .0051),
       ("real", "U"): (.0663, .0153)}
for s in ("real", "null1"):
    for a in ("O", "U"):
        R = [F[h + (s, a)] for h in HOGS]
        ev = (sum(r["fit"] != "ok" for r in R), sum(ok(r) and not .01 <= fnum(r["relax_K"]) <= 49 for r in R), len(R))
        assert abs(ev[0] / ev[2] - REP[(s, a)][0]) < 5e-5 and abs(ev[1] / ev[2] - REP[(s, a)][1]) < 5e-5, (s, a, ev)
        DIAG[(s, a)] = ev


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * (c - h), 100 * (c + h)


S3 = {(r["analysis"], r["metric"]): r for r in rtsv(W("addon", "R6", "R6_3", "r6_3_stats.tsv"))}
SETS = [("LIGHT", "Light"), ("CIRCADIAN", "Circadian"), ("PIGMENT", "Pigment")]
NSET = {"LIGHT": 450, "CIRCADIAN": 124, "PIGMENT": 119}
for k, _ in SETS:
    r = S3[("R6-3b_set_interaction", k)]
    assert int(r["n"]) == NSET[k] and fnum(r["p_holm"]) == 1, ascii(k)

# ---- layout (mm) ---------------------------------------------------------------------------------------------------
Wf, Hf = 180.0, 103.0
T1, H1, T2, H2 = 5.0, 30.0, 59.0, 32.0
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
EST = [("O", "Original"), ("U", "Uniform"), ("P", "Pooled"), ("inter", "Interaction")]
facts = [("key", "value", "source"),
         ("analysis_set", "844 of 981 built orthogroups (all four real/null1 x original/uniform fits return P and K)",
          "r6_1_fits.tsv; r6_1_report.txt")]
src = [("panel", "group", "series", "estimate", "lo95", "hi95", "n", "note")]
f6 = lambda v: "%.6f" % v
for k, (o, nm, x0) in enumerate((("relaxed", "Relaxed selection", 11.0), ("intens", "Intensified selection", 46.0))):
    ax = axes_mm(fig, Wf, Hf, x0, T1, 25.0, H1)
    v = np.array([g3(o, e) for e, _ in EST])
    ref_line(ax, 0)
    dot_whisker(ax, np.arange(4), v[:, 0], v[:, 1], v[:, 2], [l for _, l in EST], lab_fs=5)
    ax.set_ylim(*((-1, 3) if k == 0 else (-2, 2))); ax.set_yticks(np.arange(*((-1, 3.1) if k == 0 else (-2, 2.1))))
    if k == 0:
        ax.set_ylabel("Log odds ratio, deep vs null", fontsize=5.5)
    cond_box(ax, nm, fs=5)
    rr = g3(o, "R")
    ax.annotate("$R$ = %.2f (%.2f–%.2f)" % rr, xy=(.04, .04), xycoords="axes fraction", fontsize=5,
                ha="left", va="bottom")
    for (e, l), row in zip(EST + [("R", "R")], np.vstack([v, rr])):
        src.append(("a", o, e, f6(row[0]), f6(row[1]), f6(row[2]), 844, "r6_1_stats.tsv"))
    facts.append(("%s_U_O_P_inter_R" % o, "; ".join("%s %.3f (%.3f to %.3f)" % ((e,) + tuple(g3(o, e)))
                  for e in ("U", "O", "P", "inter", "R")), "r6_1_stats.tsv"))
letter(fig, Wf, Hf, "a", 1.0, T1 - 1.9)
# b: per-clade uniform-route log OR, n orthogroups in the label, dagger for clades with < 20 orthogroups
ax = axes_mm(fig, Wf, Hf, 85.0, T1, 56.0, H1)
cl = [(c, int(ST[("relaxed", "U_clade", c)]["n_hog"])) for c in PHYLO]
v = np.array([g3("relaxed", "U_clade", c) for c in PHYLO] + [g3("relaxed", "U")])
labs = ["%s (%d)%s" % (clab(c), n, "†" if n < 20 else "") for c, n in cl] + ["%s (844)" % ALL_LAB]
xs = np.r_[np.arange(8), 8.8]
ref_line(ax, 0)
dot_whisker(ax, xs, v[:, 0], v[:, 1], v[:, 2], labs, lab_fs=5)
assert v[:, 1].min() > -4 and v[:, 2].max() < 4
ax.set_ylim(-4, 4); ax.set_yticks(np.arange(-4, 4.1, 2))
ax.set_ylabel("Log odds ratio, uniform route", fontsize=5.5)
letter(fig, Wf, Hf, "b", 76.0, T1 - 1.9)
for (c, n), row in zip(cl, v):
    src.append(("b", c, "U_clade", f6(row[0]), f6(row[1]), f6(row[2]), n, "r6_1_stats.tsv"))
    facts.append(("U_clade_%s" % c, "%.3f (%.3f to %.3f), %d orthogroups" % (tuple(row) + (n,)), "r6_1_stats.tsv"))
# c: leave-one-clade-out uniform-route log OR (point estimates only; the frozen analysis wrote no CI)
ax = axes_mm(fig, Wf, Hf, 152.0, T1, 26.0, H1)
lo = np.array([fnum(ST[("relaxed", "U_leave_one_clade_out", c)]["estimate"]) for c in PHYLO])
ref_line(ax, g3("relaxed", "U")[0])
ax.annotate("All clades", xy=(3, g3("relaxed", "U")[0]), xytext=(0, -1.2), textcoords="offset points",
            ha="center", va="top", fontsize=5, color=NULL)
dot_whisker(ax, np.arange(8), lo, lo, lo, [clab(c) for c in PHYLO], lab_fs=5)
ax.set_ylim(1.0, 1.8); ax.set_yticks([1.0, 1.2, 1.4, 1.6, 1.8])
ax.set_ylabel("Log odds ratio, clade left out", fontsize=5.5)
letter(fig, Wf, Hf, "c", 143.0, T1 - 1.9)
for c, e in zip(PHYLO, lo):
    src.append(("c", c, "U_leave_one_clade_out", f6(e), "", "", ST[("relaxed", "U_leave_one_clade_out", c)]["n_hog"],
                "r6_1_stats.tsv"))
facts.append(("U_loco_range", "%.3f to %.3f" % (lo.min(), lo.max()), "r6_1_stats.tsv"))
# d: per-orthogroup log10 K, uniform vs original route
ax = axes_mm(fig, Wf, Hf, 11.0, T2, 40.0, H2)
cell_frame(ax)
leg = []
for s, col, lc, nm in (("null1", LIGHT, NULL, "Null branches"), ("real", DEEP3, DEEP, "Deep branches")):
    n, r, b1, b0 = scatter_fit(ax, LK[s][:, 0], LK[s][:, 1], col, lc, np.ones(844, bool), alpha=.7)
    leg.append(("%s ($r$ = %.2f, κ = %.2f)" % (nm, r, KAP[s]), col))
    facts.append(("logK_%s" % s, "r %.4f, slope %.4f, kappa %.7f, median U-O %.4f (clamped)" % (r, b1, KAP[s], DMED[s]),
                  "r6_1_fits.tsv [r new]"))
    for h, (x, y) in zip(A, LK[s]):
        src.append(("d", "%s|%s" % h, s, f6(y), "", "", "", "log10 K uniform; original %s" % f6(x)))
ax.set_xlim(-3.2, 2.0); ax.set_ylim(-3.2, 3.2)  # headroom for the legend (Cell Fig. 4A)
ax.set_xticks([-3, -2, -1, 0, 1]); ax.set_yticks([-3, -2, -1, 0, 1, 2])
ax.set_xlabel("log$_{10}$ $K$, original annotation", fontsize=5.5, labelpad=1.5)
ax.set_ylabel("log$_{10}$ $K$, uniform annotation", fontsize=5.5, labelpad=1.5)
cell_legend(ax, [l for l, _ in leg[::-1]], [c for _, c in leg[::-1]], fs=5)
letter(fig, Wf, Hf, "d", 1.0, T2 - 1.9)
# e: technical diagnostics by label x route (all 981 built orthogroups)
ax = axes_mm(fig, Wf, Hf, 62.0, T2, 52.0, H2)
SER = [("real", "O", "Deep, original", DEEP), ("real", "U", "Deep, uniform", DEEP2),
       ("null1", "O", "Null, original", NULL), ("null1", "U", "Null, uniform", LIGHT)]
V = np.array([[wilson(DIAG[(s, a)][j], DIAG[(s, a)][2]) for s, a, _, _ in SER] for j in (0, 1)])
grouped_bars(ax, ["Fit failure", "$K$ at bound"], V[..., 0], V[..., 1], V[..., 2], [c for *_, c in SER], gap=1.0,
             lab_fs=5.5)
assert V[..., 2].max() < 10
ax.set_ylim(0, 10); ax.set_yticks([0, 2, 4, 6, 8, 10])
ax.set_ylabel("Fits (%)", fontsize=5.5)
sq_legend(fig, 116.0 / Wf, 1 - (T2 + H2 / 2) / Hf, [l for _, _, l, _ in SER], [c for *_, c in SER], fs=5)
letter(fig, Wf, Hf, "e", 53.0, T2 - 1.9)
for j, gname in enumerate(("fit_failure", "K_at_bound")):
    for (s, a, l, _), row in zip(SER, V[j]):
        src.append(("e", gname, "%s_%s" % (s, a), f6(row[0]), f6(row[1]), f6(row[2]), DIAG[(s, a)][2],
                    "%% of built fits; events %d; Wilson 95%%" % DIAG[(s, a)][j]))
facts.append(("diagnostics", "; ".join("%s %s fail %d/981 bound %d/981" % (s, a, DIAG[(s, a)][0], DIAG[(s, a)][1])
                                       for s, a, _, _ in SER), "r6_1_fits.tsv = r6_1_report.txt"))
# f: R6-3b gene-set x deep interaction
ax = axes_mm(fig, Wf, Hf, 153.0, T2, 25.0, H2)
v = np.array([[fnum(S3[("R6-3b_set_interaction", k)][c]) for c in ("estimate", "lo", "hi")] for k, _ in SETS])
ref_line(ax, 0)
dot_whisker(ax, np.arange(3), v[:, 0], v[:, 1], v[:, 2], ["%s (%d)" % (l, NSET[k]) for k, l in SETS], lab_fs=5)
ax.set_ylim(-2, 1.5); ax.set_yticks([-2, -1, 0, 1])
ax.set_ylabel("Interaction, log odds ratio", fontsize=5.5)
ax.annotate("$P_{\\mathrm{Holm}}$ = 1", xy=(.97, .97), xycoords="axes fraction", fontsize=5, ha="right", va="top")
letter(fig, Wf, Hf, "f", 143.0, T2 - 1.9)
for (k, _), row in zip(SETS, v):
    src.append(("f", k, "set_x_deep", f6(row[0]), f6(row[1]), f6(row[2]), NSET[k], "r6_3_stats.tsv; Holm P 1"))
    facts.append(("R6_3b_%s" % k, "%.3f (%.3f to %.3f), n %d, Holm P 1" % (tuple(row) + (NSET[k],)), "r6_3_stats.tsv"))
res = qa_v2(fig)
clean = qa_report(res, "Fig3")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig3_qa_v2.txt"))
write_facts("fig3_source_v2.tsv", src)
write_facts("fig3_facts_v2.tsv", facts)
save_v2(fig, "Fig3")
print("Fig3 written; QA clean:", clean)
