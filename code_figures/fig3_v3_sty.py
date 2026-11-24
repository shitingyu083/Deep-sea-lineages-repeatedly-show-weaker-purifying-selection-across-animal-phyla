"""Fig. 3 v3 - the relaxation contrast under one uniform annotation route (R6-1) and the functional-concentration
test (R6-3b). Every panel copies the form of one viewed published panel (panel map: paper/figure_benchmark_map.md,
"Fig. 3 v3"); no form is invented and no panel mixes two benchmarks.
  a     = Bergeron et al. 2023 Nature Fig. 1c data column, drawn twice (relaxed above, intensified below):
          half-step light stripes, thin reference line, two-line header with one double-headed arrow, black bars
          from the no-difference value with grey capped 95% intervals; rows = four estimands, upright names with
          dotted leaders, each column with its own ticks and x label
  b     = Bergeron Fig. 1c in full (time tree, upright names, silhouettes, dotted leaders, one data column)
  c     = Bergeron Fig. 1c data column sharing b's rows, bars from the all-clade estimate, no intervals
  d     = Bergeron Fig. 2a,b: one scatter per branch set in its own axes, grey thick dashed least-squares line across
          the axes, centred multi-line statistics, one-row dot legend below
  e     = Wei et al. 2026 Nature Extended Data Fig. 4: annotated heatmap, integer counts in the tiles, zero cells
          grey, red-white-blue log scale about the geometric midpoint, row names right-aligned with the 981-fit
          totals in parentheses, column names at 35 degrees, colour bar at the right
  f     = Bergeron Fig. 1c data column: the three gene sets (Holm P = 1 goes to the caption and Source Data)
Data (read-only): R6/r6_1_stats.tsv f4760668, R6/r6_1_fits.tsv 722b737e, R6/R6_3/r6_3_stats.tsv ca2600c9. Every v2
assert is kept; e is new and its eight clade columns are asserted to sum to the v2 totals 74/20, 53/7, 65/15, 50/5.
Outputs in viz_v3_sty/out only.
"""
from figcommon_v3_sty import *  # noqa: F401,F403
from matplotlib.colors import LinearSegmentedColormap

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
# e: the same two diagnostics split into the eight clade columns; the columns must add up to the v2 totals
CLN = {c: sum(1 for (k, _) in HOGS if k == c) for c in PHYLO}
DGC = {}
for cl in PHYLO:
    assert CLN[cl] > 0
    for s in ("real", "null1"):
        for a in ("O", "U"):
            R = [F[(k, h, s, a)] for (k, h) in HOGS if k == cl]
            DGC[(s, a, cl)] = (sum(r["fit"] != "ok" for r in R),
                               sum(ok(r) and not .01 <= fnum(r["relax_K"]) <= 49 for r in R), len(R))
assert sum(CLN.values()) == 981
for s in ("real", "null1"):
    for a in ("O", "U"):
        for j in (0, 1):
            assert sum(DGC[(s, a, cl)][j] for cl in PHYLO) == DIAG[(s, a)][j], (s, a, j)


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * p, 100 * (c - h), 100 * (c + h)


S3 = {(r["analysis"], r["metric"]): r for r in rtsv(W("addon", "R6", "R6_3", "r6_3_stats.tsv"))}
SETS = [("LIGHT", "Light"), ("CIRCADIAN", "Circadian"), ("PIGMENT", "Pigment")]
NSET = {"LIGHT": 450, "CIRCADIAN": 124, "PIGMENT": 119}
for k, _ in SETS:
    r = S3[("R6-3b_set_interaction", k)]
    assert int(r["n"]) == NSET[k] and fnum(r["p_holm"]) == 1, ascii(k)

# ---- layout (mm) ---------------------------------------------------------------------------------------------------
Wf, Hf = 180.0, 144.0
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
ax = mm_axes(fig, Wf, Hf)
f6 = lambda v: "%.6f" % v
src = [("panel", "group", "series", "estimate", "lo95", "hi95", "n", "note")]
facts = [("key", "value", "source"),
         ("analysis_set", "844 of 981 built orthogroups (all four real/null1 x original/uniform fits return P and K)",
          "r6_1_fits.tsv; r6_1_report.txt")]
EST = [("O", "Original"), ("U", "Uniform"), ("P", "Pooled"), ("inter", "Interaction")]
OUT = [("relaxed", "Relaxed call", -1, 3, 6.0), ("intens", "Intensified call", -2, 2, 6.0)]
COLA = {"relaxed": DEEP, "intens": DEEP2}            # matches Fig. 2a: relaxed DEEP, intensified DEEP2

# ---- a: Berg 1c data column, twice (relaxed above, intensified below) ----------------------------------------------
PA = 5.2
XS_A = 22.0                                      # spine (axis start) of both columns
XVA = 49.5                                       # right-aligned value column, as the berg 1a count boxes
YA = {"relaxed": 20.0, "intens": 56.0}
YAX_A = {k: v + 3 * PA + 2.6 for k, v in YA.items()}
for (o, nm, lo, hi, unit) in OUT:
    y0 = YA[o]
    ya = YAX_A[o]
    ys = [y0 + r * PA for r in range(4)]
    u = berg_col(ax, XS_A, lo, hi, unit, (y0 - 0.5 * PA, ya), ya, y0 - 8.2, y0 - 4.9,
                 head1="No difference", head2=("Null higher", "Deep higher"), xlabel="%s, log odds ratio" % nm,
                 ticks=np.arange(lo, hi + .01, 1))
    for r, (e, lab) in enumerate(EST):
        name_row(ax, 3.0, ys[r], lab, XS_A)
    v = np.array([g3(o, e) for e, _ in EST])
    assert v[:, 1].min() >= lo and v[:, 2].max() <= hi, (o, v)
    berg_bars_c(ax, u, ys, v[:, 0], v[:, 1], v[:, 2], COLA[o], pitch=PA)
    for y, row in zip(ys, v):
        ax.text(XVA, y, "%.2f" % row[0], ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
    for (e, lab), row in zip(EST, v):
        src.append(("a", o, e, f6(row[0]), f6(row[1]), f6(row[2]), 844, "r6_1_stats.tsv; bars from 0"))
    facts.append(("%s_U_O_P_inter_R" % o, "; ".join("%s %.3f (%.3f to %.3f)" % ((e,) + tuple(g3(o, e)))
                  for e in ("U", "O", "P", "inter", "R")), "r6_1_stats.tsv"))
letter(fig, Wf, Hf, "a", 1.0, 6.5)

# ---- b: Berg 1c in full - tree, names, silhouettes, one data column -----------------------------------------------
PB = 5.6
YB0 = 18.0
yb = lambda r: YB0 + r * PB
YB_ALL = yb(7) + 0.6 * PB
YAX_B = YB_ALL + 0.5 * PB
rows_b = (YB0 - 0.5 * PB, YAX_B)
TB = focal_tree(PHYLO)
YB = {**{c: yb(r) for r, c in enumerate(PHYLO)}, **{clab(c): yb(r) for r, c in enumerate(PHYLO)}}
SIL = {clab(c): min(0.9 * PB, 6.0 / aspect(PNG[clab(c)])) for c in PHYLO}
draw_tree(ax, TB, 52.0, 59.5, YB)
cl = [(c, int(ST[("relaxed", "U_clade", c)]["n_hog"])) for c in PHYLO]
XS_B = 88.0
XVB = 117.0
COLS = [ccol(c) for c in PHYLO] + [ccol(ALL_LAB)]
for r, (c, n) in enumerate(cl):
    name_row_c(ax, 60.5, yb(r), "%s (%d)%s" % (clab(c), n, "†" if n < 20 else ""), XS_B, clab(c),
               PNG[clab(c)], h=SIL[clab(c)])
name_row_c(ax, 60.5, YB_ALL, "%s (844)" % ALL_LAB, XS_B, ALL_LAB, None)
UB = berg_col(ax, XS_B, -4, 4, 2.9, rows_b, YAX_B, YB0 - 8.2, YB0 - 4.9, head1="No difference",
              head2=("Null higher", "Deep higher"), arrow_to=0.4, xlabel="Log odds ratio, uniform route",
              ticks=[-4, -2, 0, 2, 4])
v = np.array([g3("relaxed", "U_clade", c) for c in PHYLO] + [g3("relaxed", "U")])
assert v[:, 1].min() > -4 and v[:, 2].max() < 4
ys_b = [yb(r) for r in range(8)] + [YB_ALL]
berg_bars_c(ax, UB, ys_b, v[:, 0], v[:, 1], v[:, 2], COLS, pitch=PB)
for y, row in zip(ys_b, v):
    ax.text(XVB, y, "%.2f" % row[0], ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "b", 46.5, 6.5)
for (c, n), row in zip(cl, v):
    src.append(("b", c, "U_clade", f6(row[0]), f6(row[1]), f6(row[2]), n, "r6_1_stats.tsv"))
    facts.append(("U_clade_%s" % c, "%.3f (%.3f to %.3f), %d orthogroups" % (tuple(row) + (n,)), "r6_1_stats.tsv"))
src.append(("b", "ALL", "U", f6(v[8, 0]), f6(v[8, 1]), f6(v[8, 2]), 844, "r6_1_stats.tsv"))

# ---- c: Berg 1c data column, bars from the all-clade estimate, no intervals ---------------------------------------
LOCO = np.array([fnum(ST[("relaxed", "U_leave_one_clade_out", c)]["estimate"]) for c in PHYLO])
REF = g3("relaxed", "U")[0]
XVC = 140.5
UC = berg_col(ax, 122.0, 1.0, 1.6, 22.0, rows_b, YAX_B, YB0 - 8.2, YB0 - 4.9, ref=REF, head1="All clades",
              head2=("Lower", "Higher"), arrow_to=1.45, xlabel="Log odds ratio, clade left out",
              ticks=[1.0, 1.2, 1.4, 1.6])
assert LOCO.min() >= 1.0 and LOCO.max() <= 1.9, LOCO
berg_bars_c(ax, UC, ys_b[:8], LOCO, [None] * 8, [None] * 8, COLS[:8], base=REF, pitch=PB)
for y, e in zip(ys_b[:8], LOCO):
    ax.text(XVC, y, "%.2f" % e, ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "c", 110.0, 6.5)
for c, e in zip(PHYLO, LOCO):
    src.append(("c", c, "U_leave_one_clade_out", f6(e), "", "", ST[("relaxed", "U_leave_one_clade_out", c)]["n_hog"],
                "r6_1_stats.tsv; bars from %.6f; no CI in the frozen analysis" % REF))
facts.append(("U_loco_range", "%.3f to %.3f (all clades %.3f)" % (LOCO.min(), LOCO.max(), REF), "r6_1_stats.tsv"))

# ---- f: Berg 1c data column, the three gene sets -------------------------------------------------------------------
YF0 = 34.0
XVF = 178.0
ys_f = [YF0 + r * PB for r in range(3)]
UF = berg_col(ax, 162.0, -2, 1, 3.2, (YF0 - 0.5 * PB, YF0 + 2 * PB + 2.6), YF0 + 2 * PB + 2.6, YF0 - 8.2,
              YF0 - 4.9, head1="No difference", head2=("Lower", "Higher"), arrow_to=0.4,
              xlabel="Interaction, log odds ratio", ticks=[-2, -1, 0, 1])
VF = np.array([[fnum(S3[("R6-3b_set_interaction", k)][c]) for c in ("estimate", "lo", "hi")] for k, _ in SETS])
assert VF[:, 1].min() >= -2 and VF[:, 2].max() <= 1, VF
berg_bars_c(ax, UF, ys_f, VF[:, 0], VF[:, 1], VF[:, 2], DEEP, pitch=PB)
for y, row in zip(ys_f, VF):
    ax.text(XVF, y, "%.2f" % row[0], ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
for (k, lab), y, row in zip(SETS, ys_f, VF):
    name_row(ax, 143.5, y, "%s (%d)" % (lab, NSET[k]), 162.0)
    src.append(("f", k, "set_x_deep", f6(row[0]), f6(row[1]), f6(row[2]), NSET[k], "r6_3_stats.tsv; Holm P 1"))
    facts.append(("R6_3b_%s" % k, "%.3f (%.3f to %.3f), n %d, Holm P 1" % (tuple(row) + (NSET[k],)),
                  "r6_3_stats.tsv"))
letter(fig, Wf, Hf, "f", 141.0, 22.5)

# ---- d: Berg 2a,b - one scatter per branch set, own axes, dashed fit, centred statistics, dot legend --------------
DX, DW = (13.0, 49.0), 30.0
YD0, YBOT = 86.0, 124.0
for i, (s, nm, col) in enumerate((("real", "Deep branches", DEEP2), ("null1", "Null branches", NULL))):
    axd = axes_mm(fig, Wf, Hf, DX[i], YD0, DW, YBOT - YD0)
    x, y = LK[s][:, 0], LK[s][:, 1]
    axd.scatter(x, y, s=1.7, c=col, lw=0, alpha=.62 if s == "real" else .55, zorder=2)
    axd.set_xlim(-3.2, 1.8); axd.set_ylim(-3.2, 1.8)
    axd.set_xticks([-3, -2, -1, 0, 1]); axd.set_yticks([-3, -2, -1, 0, 1])
    axd.tick_params(labelsize=5, pad=1.5)
    b1, b0 = np.polyfit(x, y, 1)
    xx = np.array([-3.2, 1.8])
    axd.plot(xx, b1 * xx + b0, color=NULL, lw=1.3, ls=(0, (4, 2)), solid_capstyle="butt", zorder=3)
    axd.text(.32, .965, "%s\n$r$ = %.2f\n$\\kappa$ = %.2f" % (nm, R_[s], KAP[s]), transform=axd.transAxes,
             ha="center", va="top", fontsize=5.5, color=INK, linespacing=1.35, zorder=10)
    axd.set_xlabel("log$_{10}$ $K$, original", fontsize=5.5, labelpad=1.5)
    axd.set_ylabel("log$_{10}$ $K$, uniform", fontsize=5.5, labelpad=1.5)
    facts.append(("logK_%s" % s, "r %.4f, slope %.4f, kappa %.7f, median U-O %.4f (clamped)" % (R_[s], b1, KAP[s],
                  DMED[s]), "r6_1_fits.tsv [r new]"))
    for h, (xv, yv) in zip(A, LK[s]):
        src.append(("d", "%s|%s" % h, s, f6(yv), "", "", "", "log10 K uniform; original %s" % f6(xv)))
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=2.6, mfc=c, mec="none", label=l)
                   for l, c in (("Deep branches", DEEP2), ("Null branches", NULL))],
          loc="upper center", bbox_to_anchor=(31.0 / Wf, 1 - 132.0 / Hf), ncol=2, frameon=False, fontsize=5.5,
          handletextpad=.3, columnspacing=1.2, borderpad=0)
letter(fig, Wf, Hf, "d", 1.0, 82.0)

# ---- e: terr ED Fig. 4 annotated heatmap ---------------------------------------------------------------------------
SER = [("real", "O", "deep, original"), ("real", "U", "deep, uniform"),
       ("null1", "O", "null, original"), ("null1", "U", "null, uniform")]
DNAME = ("Fit failure", "$K$ at bound")
EZ = np.array([[DGC[(s, a, c)][j] for c in PHYLO] for j in (0, 1) for s, a, _ in SER], float)   # 8 rows x 8 clades
assert EZ.shape == (8, 8), EZ.shape
ER = ["%s, %s (%d)" % (DNAME[j], lab, DIAG[(s, a)][j]) for j in (0, 1) for s, a, lab in SER]
EC = [clab(c) for c in PHYLO]
nz = EZ[EZ > 0]
EMID = float(np.sqrt(nz.min() * nz.max()))
EFILL = np.where(EZ > 0, np.log10(np.maximum(EZ, 1) / EMID), np.nan)
EV = float(np.nanmax(np.abs(EFILL)))
assert np.isfinite(EV) and EV > 0
CM = LinearSegmentedColormap.from_list("ed9z", ED9_STOPS).copy()
CM.set_bad("#b4b4b4")
axe, stub = ed9_axes(fig, Wf, Hf, 118.0, 84.0, 8, 8, 5.6, 5.6, stub=.12)
ed9_heat(axe, EFILL, [["%d" % EZ[i, j] for j in range(8)] for i in range(8)], ER, EC, EV, stub, fs=5.5, lab_fs=5.0,
         rot=35, cmap=CM)
cax = axes_mm(fig, Wf, Hf, 166.5, 96.0, 4.6, 22.0)
ed9_cbar(cax, EV, [np.log10(1 / EMID), 0.0, np.log10(nz.max() / EMID)], ["%d" % nz.min(), "%d" % round(EMID),
         "%d" % nz.max()], "Number of fits", fs=5, title_fs=5.5, cmap=CM)
letter(fig, Wf, Hf, "e", 86.0, 82.0)
for j, dn in enumerate(("fit_failure", "K_at_bound")):
    for i, (s, a, lab) in enumerate(SER):
        for c in PHYLO:
            k_, n_ = DGC[(s, a, c)][j], DGC[(s, a, c)][2]
            assert EZ[4 * j + i, PHYLO.index(c)] == k_
            src.append(("e", "%s|%s" % (dn, lab.replace(", ", "_")), c, "%d" % k_, "", "", n_,
                        "fits with that event of %d built fits" % n_))
            if k_:
                q, l, h = wilson(k_, n_)
                src.append(("e_pct", "%s|%s" % (dn, lab.replace(", ", "_")), c, "%.6f" % q, "%.6f" % l,
                            "%.6f" % h, n_, "percent of built fits; Wilson 95%"))
    facts.append(("diag_%s_by_clade" % dn, "; ".join("%s %d of %d" % (clab(c), DGC[("real", "U", c)][j], CLN[c])
                                                     for c in PHYLO), "r6_1_fits.tsv [clade split, new in v3]"))
facts.append(("diagnostics", "; ".join("%s %s fail %d/981 bound %d/981" % (s, a, DIAG[(s, a)][0], DIAG[(s, a)][1])
                                       for s, a, _ in SER), "r6_1_fits.tsv = r6_1_report.txt"))

res = qa_v2(fig)
clean = qa_report(res, "Fig3")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig3_qa_v2.txt"))
write_facts("fig3_source_v3.tsv", src)
write_facts("fig3_facts_v3.tsv", facts)
save_v2(fig, "Fig3")
print("Fig3 written; QA clean:", clean)
print("logK r", {k: round(v, 4) for k, v in R_.items()}, "kappa", {k: round(v, 7) for k, v in KAP.items()})
print("loco", np.round(LOCO, 3), "ref", round(REF, 3))
print("E mid", round(EMID, 2), "range", nz.min(), nz.max(), "EV", round(EV, 3))
print("E rows", EZ.astype(int).tolist())
