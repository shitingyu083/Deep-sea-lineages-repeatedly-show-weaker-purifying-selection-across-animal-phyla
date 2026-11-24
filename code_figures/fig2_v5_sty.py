"""Fig. 2 v5 - core relaxation effect in the R6-3 sample (production fits; 1,000 deep and 400 matched null branch fits
per clade). Statistics block = v4 (unchanged, asserted against R); drawing redone in denser published forms
(panel map: paper/figure_benchmark_map.md, "Fig. 2-4 v5"):
  a = lam6d, Nature 2022 Extended Data Fig. 6d: faceted volcano grid; columns = clades, rows = deep / matched null
      foreground; RELAX log10 K against -log10 P of every fit, calls (P < 0.05) in the clade colour, other fits grey.
  b = stil2b, Stiller et al. 2024 Nature Fig. 2b: time tree with clade-coloured terminal branches and names under them,
      one light violin per clade of the paired per-orthogroup difference log10 K(deep) - log10 K(matched null) with a
      darker clade-median bar, silhouettes to the right, leader-line labels.
  c = berg 1c data columns (v3/v4 Fig. 2b,c): relaxed and intensified log odds ratios, deep vs matched null.
Data: R6/R6_3/r6_3_fits.tsv (read-only). Outputs: viz_v5_sty/out only."""
from figcommon_v5_sty import *  # noqa: F401,F403

mt_arial()
rows = rtsv(W("addon", "R6", "R6_3", "r6_3_fits.tsv"))
assert len(rows) == 11200, len(rows)
FIT, KV = {}, {}
for r in rows:
    p, K = fnum(r["relax_p"]), fnum(r["relax_K"]); ok = bool(np.isfinite(p) and np.isfinite(K))
    k = (r["clade"], r["label"], r["hog"])
    assert k not in FIT and r["label"] in ("real", "null1"), ascii(k)
    FIT[k] = (int(ok and p < .05 and K < 1), int(ok and p < .05 and K > 1), fnum(r["has_go"]) == 1, ok)
    if ok:
        KV.setdefault((r["clade"], r["label"]), []).append(K)
KLO, KHI = -3.0, np.log10(50.0)   # axis floor (K <= 0.001, incl. three K = 0 fits, drawn at the floor) and RELAX bound
LK = {k: np.clip(np.log10(np.maximum(np.asarray(v), 1e-300)), KLO, KHI) for k, v in KV.items()}
assert sum(len(v) for v in LK.values()) == 11200 - 10, sum(len(v) for v in LK.values())
NFLOOR = sum(int((v <= KLO).sum()) for v in LK.values())
EXP = {  # deep relaxed, deep intensified, null relaxed, null intensified (r6g rates in stats_nohup.txt x 1000 / 400)
    "Porifera": (120, 147, 54, 39), "Cnidaria": (530, 18, 62, 54), "Annelida": (315, 117, 56, 78),
    "Mollusca": (641, 11, 70, 62), "Arthropoda": (218, 134, 81, 65), "Echinodermata": (430, 13, 49, 46),
    "Chondrichthyes": (148, 153, 53, 51), "Actinopteri_m": (206, 271, 53, 81)}
B, C = 2000, len(PHYLO)
rng = np.random.default_rng(20261031)
E1, E0 = np.zeros((B + 1, C, 2)), np.zeros((B + 1, C, 2))   # row 0 = observed, rows 1..B = resamples
N1, N0 = np.zeros((B + 1, C)), np.zeros((B + 1, C))
for c, cl in enumerate(PHYLO):
    U = sorted(h for (k, l, h) in FIT if k == cl and l == "real")
    Nh = {h for (k, l, h) in FIT if k == cl and l == "null1"}
    assert len(U) == 1000 and len(Nh) == 400 and Nh <= set(U), ascii(cl)
    y1 = np.array([FIT[(cl, "real", h)][:2] for h in U], float)
    y0 = np.array([FIT[(cl, "null1", h)][:2] if h in Nh else (0, 0) for h in U], float)
    m0 = np.array([h in Nh for h in U], float)
    got = tuple(int(v) for v in np.r_[y1.sum(0), y0.sum(0)])
    assert got == EXP[cl], ascii((cl, got))
    w = np.vstack([np.ones(len(U)), boot_w(rng, len(U), B)])
    E1[:, c], N1[:, c], E0[:, c], N0[:, c] = w @ y1, w.sum(1), w @ y0, w @ m0
R1, R0 = 100 * E1 / N1[..., None], 100 * E0 / N0[..., None]   # % of orthogroup fits, (B + 1, C, outcome)
RA1, RA0 = R1.mean(1), R0.mean(1)                             # All clades: mean of the eight clade percentages
LC = logor2(E1, N1[..., None], E0, N0[..., None])             # per-clade log OR, deep vs matched null
LA = np.stack([fe_logit(E1[..., o], N1, E0[..., o], N0) for o in (0, 1)], -1)   # clade-adjusted log OR (B + 1, 2)
for o, ref in enumerate((1.1004475620, -0.3834714319)):        # R glm on all 11,200 fits (xcheck_fig2_sty.txt)
    assert abs(LA[0, o] - ref) < 1e-6, ascii((o, float(LA[0, o])))
REF_LC = {"Porifera": (-0.1349754362, 0.4669893515), "Cnidaria": (1.8160558223, -2.1417648220),
          "Annelida": (1.0384437672, -0.6033085470), "Mollusca": (2.1303044808, -2.8028875484),
          "Arthropoda": (0.0933822704, -0.2263018467), "Echinodermata": (1.6871147732, -2.2890651653),
          "Chondrichthyes": (0.1286586133, 0.2119835162), "Actinopteri_m": (0.5298255739, 0.3811870370)}
for c, cl in enumerate(PHYLO):
    assert np.allclose(LC[0, c], REF_LC[cl], rtol=0, atol=1e-8), ascii(cl)
HG = np.zeros((C, 2, 2, 2))   # clade x (deep, null) x outcome x (events, fits): fits of has_go orthogroups with P and K
for (cl, l, h), (yr, yi, g, ok) in FIT.items():
    if g and ok:
        s = HG[PHYLO.index(cl), int(l == "null1")]
        s[:, 0] += (yr, yi); s[:, 1] += 1
assert HG[:, :, 0, 1].sum() == 9208, HG[:, :, 0, 1].sum()
LG = [float(fe_logit(*(HG[:, k, o, j] for k in (0, 1) for j in (0, 1)))) for o in (0, 1)]
for o, ref in enumerate((1.0497951406, -0.3496607684)):        # relaxed = r6g stats_nohup.txt 1.049795
    assert abs(LG[o] - ref) < 1e-6, ascii((o, LG[o]))


def summ(A):  # observed, 95% percentile CI and two-sided bootstrap P of an array whose row 0 is observed
    lo, hi = pci(A[1:])
    return A[0], lo, hi, bp(A[1:])


r1, r1l, r1h, _ = summ(R1); r0, r0l, r0h, _ = summ(R0)
a1, a1l, a1h, _ = summ(RA1); a0, a0l, a0h, _ = summ(RA0)
lc, lcl, lch, lcp = summ(LC); la, lal, lah, lap = summ(LA)
cells = np.stack([E1[1:], N1[1:, :, None] - E1[1:], E0[1:], N0[1:, :, None] - E0[1:]])
NZ = (cells.min(0) <= 0).sum(0)                                # resamples with a zero cell (Haldane +0.5), (C, 2)
NOPK = {l: sum(1 for (k, ll, h), v in FIT.items() if ll == l and not v[3]) for l in ("real", "null1")}
for o, nm in enumerate(("relaxed", "intensified")):
    print(nm, "pooled log OR %.4f [%.4f, %.4f] P=%.4f; has_go %.4f; zero-cell resamples %d" % (
        la[o], lal[o], lah[o], lap[o], LG[o], NZ[:, o].sum()))
print("fits without P or K:", NOPK)
OUT, LAB = ("relaxed", "intensified"), [clab(c) for c in PHYLO] + [ALL_LAB]
f6 = lambda v: "%.6f" % v
src = [("panel", "outcome", "group", "series", "events", "fits", "estimate", "lo95", "hi95", "p_boot", "note")]
for o, nm in enumerate(OUT):
    pa = "a"                                                  # v3: both percentages in one Cell 2B panel
    for c, cl in enumerate(PHYLO):
        src.append((pa, nm, cl, "deep", int(E1[0, c, o]), int(N1[0, c]), f6(r1[c, o]), f6(r1l[c, o]), f6(r1h[c, o]),
                    "", "% of fits"))
        src.append((pa, nm, cl, "matched_null", int(E0[0, c, o]), int(N0[0, c]), f6(r0[c, o]), f6(r0l[c, o]),
                    f6(r0h[c, o]), "", "% of fits"))
    src.append((pa, nm, "ALL", "deep", int(E1[0, :, o].sum()), int(N1[0].sum()), f6(a1[o]), f6(a1l[o]), f6(a1h[o]),
                "", "mean of 8 clade %"))
    src.append((pa, nm, "ALL", "matched_null", int(E0[0, :, o].sum()), int(N0[0].sum()), f6(a0[o]), f6(a0l[o]),
                f6(a0h[o]), "", "mean of 8 clade %"))
    for c, cl in enumerate(PHYLO):
        src.append((("bc"[o]), nm, cl, "log_OR", "", "", f6(lc[c, o]), f6(lcl[c, o]), f6(lch[c, o]),
                    "%.4f" % lcp[c, o], "2x2 deep vs matched null; zero-cell resamples %d" % NZ[c, o]))
    src.append((("bc"[o]), nm, "ALL", "log_OR", "", "", f6(la[o]), f6(lal[o]), f6(lah[o]), "%.4f" % lap[o],
                "glm y ~ deep + clade, all 11200 fits; has_go subset (9208 fits) %.6f" % LG[o]))
SRC2 = src  # written at the end (v5 adds rows)
pct = lambda v, l, h: "%.1f%% (%.1f-%.1f)" % (v, l, h)
facts = [("key", "value", "source"), ("fits_per_clade", "1,000 deep and 400 matched null (null1) branch fits; "
          "8 clades; 11,200 fits", "r6_3_fits.tsv"),
         ("no_P_or_K", "deep %d, null %d (counted as neither)" % (NOPK["real"], NOPK["null1"]), "r6_3_fits.tsv"),
         ("bootstrap", "2,000 resamples of orthogroups within clade; numpy default_rng(20261031) [new: numpy stream]",
          "fig2_v5_sty.py (statistics block = v4)")]
for o, nm in enumerate(OUT):
    facts.append(("%s_all_clades" % nm, "deep %s vs null %s" % (pct(a1[o], a1l[o], a1h[o]), pct(a0[o], a0l[o],
                  a0h[o])), "fig2_source_v5.tsv [new CIs]"))
    for c, cl in enumerate(PHYLO):
        facts.append(("%s_%s" % (nm, cl), "deep %s vs null %s; log OR %.2f (%.2f to %.2f), %s" % (
            pct(r1[c, o], r1l[c, o], r1h[c, o]), pct(r0[c, o], r0l[c, o], r0h[c, o]), lc[c, o], lcl[c, o],
            lch[c, o], pfmt(lcp[c, o]).replace("$", "")), "fig2_source_v5.tsv [new CIs, P]"))
    facts.append(("%s_pooled_log_OR" % nm, "%.2f (%.2f to %.2f), %s; has_go subset %.4f" % (
        la[o], lal[o], lah[o], pfmt(lap[o]).replace("$", ""), LG[o]), "fig2_source_v5.tsv [new CI, P]"))

from scipy.stats import gaussian_kde


def sv(v, d=1):
    """Value with a true minus sign (U+2212) at `d` decimals."""
    return (MINUS if v < 0 else "") + ("%." + str(d) + "f") % abs(v)


def star(p):
    """Bootstrap significance keyed to the manuscript: *** P<0.01, ** P<0.05, * P<0.1."""
    return "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else ""


# ---- per-fit values for a and paired orthogroup differences for b --------------------------------------------------
LP_CAP = 10.0                                    # -log10 P drawn at most 10 (P <= 1e-10, incl. P = 0, at the top)
VOL, LKH = {}, {}
for r in rows:
    p, K = fnum(r["relax_p"]), fnum(r["relax_K"])
    if not (np.isfinite(p) and np.isfinite(K)):
        continue
    lk = float(np.log10(min(max(K, 1e-3), 50.0)))
    VOL.setdefault((r["clade"], r["label"]), []).append((lk, min(-np.log10(max(p, 1e-300)), LP_CAP), p < 0.05, K))
    LKH[(r["clade"], r["label"], r["hog"])] = lk
for c in PHYLO:
    for l, k in (("real", 0), ("null1", 2)):
        v = VOL[(c, l)]
        nr = sum(1 for x in v if x[2] and x[3] < 1); ni = sum(1 for x in v if x[2] and x[3] > 1)
        assert (nr, ni) == EXP[c][k:k + 2], (c, l, nr, ni)
PAIR = {}
for c in PHYLO:
    hs = sorted(h for (cc, l, h) in LKH if cc == c and l == "null1" and (c, "real", h) in LKH)
    PAIR[c] = np.array([LKH[(c, "real", h)] - LKH[(c, "null1", h)] for h in hs])
assert sum(len(v) for v in PAIR.values()) == 3196, sum(len(v) for v in PAIR.values())
for c in PHYLO:                                  # share of significant deep fits that lie at K < 1 (Fig. 2a)
    nr, ni = EXP[c][0], EXP[c][1]
    facts.append(("sigK<1_%s" % c, "%d of %d significant deep fits at K < 1 (%.1f%%)" % (nr, nr + ni, 100 * nr / (nr + ni)),
                  "r6_3_fits.tsv (P < 0.05)"))

Wf, Hf = 180.0, 133.0
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
ax = mm_axes(fig, Wf, Hf)
src = SRC2

# ---- a: lam6d faceted volcano grid --------------------------------------------------------------------------------
VX0, VY0, PW, PH, GX, GY, STRIP = 14.0, 10.6, 19.0, 16.4, 1.05, 1.2, 3.4
XLO, XHI, YLO, YHI = -3.2, np.log10(50.0) + 0.2, -0.7, 12.9
ROWS = (("real", "Deep"), ("null1", "Matched null"))
NS_COL = "#c4c4c4"
for j, c in enumerate(PHYLO):
    x0 = VX0 + j * (PW + GX)
    ax.add_patch(Rectangle((x0, VY0 - 0.7 - STRIP), PW, STRIP, fc="white", ec="k", lw=0.6, zorder=8))
    ax.text(x0 + PW / 2, VY0 - 0.7 - STRIP / 2, clab(c), ha="center", va="center", fontsize=5.5, color=INK,
            zorder=9)._qa_skip = True
    ci = PHYLO.index(c)
    for i, (lab, nm) in enumerate(ROWS):
        y0 = VY0 + i * (PH + GY)
        a = axes_mm(fig, Wf, Hf, x0, y0, PW, PH)
        v = np.array([x[:3] for x in VOL[(c, lab)]], float)
        sig = v[:, 2] > 0
        a.scatter(v[~sig, 0], v[~sig, 1], s=0.8, c=NS_COL, alpha=0.65, lw=0, edgecolors="none", rasterized=True,
                  zorder=2)
        a.scatter(v[sig, 0], v[sig, 1], s=0.8, c=ccol(c), alpha=0.85, lw=0, edgecolors="none", rasterized=True,
                  zorder=3)
        a.set_xlim(XLO, XHI); a.set_ylim(YLO, YHI); a.set_facecolor("none")
        for sp in a.spines.values():
            sp.set_visible(False)
        if j == 0:
            a.spines["left"].set_visible(True); a.spines["left"].set_linewidth(0.5)
            a.spines["left"].set_bounds(0, LP_CAP)
            a.set_yticks([0, 5, 10]); a.set_yticklabels(["0", "5", "≥10"])
            a.tick_params(axis="y", length=2, width=0.5, pad=1.0, labelsize=5.0, colors=INK)
        else:
            a.set_yticks([])
        if i == 1:
            a.spines["bottom"].set_visible(True); a.spines["bottom"].set_linewidth(0.5)
            a.spines["bottom"].set_bounds(-3, 1)
            a.set_xticks([-3, -2, -1, 0, 1]); a.set_xticklabels(["", "0.01", "", "1", "10"])
            a.tick_params(axis="x", length=2, width=0.5, pad=1.0, labelsize=5.0, colors=INK)
        else:
            a.set_xticks([])
        pr = (r1 if lab == "real" else r0)[ci, 0]
        pi = (r1 if lab == "real" else r0)[ci, 1]
        a.text(XLO + 0.08, YHI - 0.3, "%.1f%%" % pr, ha="left", va="top", fontsize=5.0, color=ccol(c), zorder=10)
        a.text(XHI - 0.08, YHI - 0.3, "%.1f%%" % pi, ha="right", va="top", fontsize=5.0, color=INK, zorder=10)
        ncap = int((v[:, 1] >= LP_CAP).sum())
        src.append(("a", "volcano", c, "deep" if lab == "real" else "matched_null", int(sig.sum()), len(v),
                    "%.1f" % pr, "%.1f" % pi, "", "", "fits with P and K; P<0.05 coloured; %d at the -log10 P cap 10; "
                    "relaxed %% left, intensified %% right" % ncap))
XS_R = VX0 + 8 * PW + 7 * GX + 0.7
for i, (lab, nm) in enumerate(ROWS):
    y0 = VY0 + i * (PH + GY)
    ax.add_patch(Rectangle((XS_R, y0), STRIP, PH, fc="white", ec="k", lw=0.6, zorder=8))
    ax.text(XS_R + STRIP / 2, y0 + PH / 2, nm, rotation=270, ha="center", va="center", fontsize=5.5, color=INK,
            zorder=9)._qa_skip = True
YVB = VY0 + 2 * PH + GY
ax.text(4.2, VY0 + PH + GY / 2, "$-\\log_{10}$ RELAX $P$", rotation=90, ha="center", va="center", fontsize=6.0,
        color=INK, zorder=10)
ax.text(VX0 + (8 * PW + 7 * GX) / 2, YVB + 4.6, "RELAX selection-intensity parameter $K$ (log scale)", ha="center",
        va="top", fontsize=6.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "a", 1.0, 4.4)

# ---- b: stil2b, time tree + one violin per clade of the paired orthogroup difference ---------------------------------
PB, YB0 = 6.9, 65.4
yb = lambda r: YB0 + r * PB
YB_ALL = yb(7) + 0.62 * PB
YAX = YB_ALL + 0.55 * PB                       # shared baseline of b (time and violin axes) and c
TB = focal_tree(PHYLO)
TX0, TX1 = 2.0, 30.0                           # root, present


def light(hexcol, f=0.45):
    """Stiller's pale terminal-branch tint: the clade colour mixed towards white."""
    c = np.asarray(to_rgb(hexcol))
    return tuple(c + (1 - c) * f)


def draw_tree_stil(root, x0, x1, ymap):
    A = root.age
    X = lambda a_: x0 + (A - a_) / A * (x1 - x0)

    def sety(q):
        if not q.kids:
            q.y = ymap[q.clade]
        else:
            for k in q.kids:
                sety(k)
            q.y = 0.5 * (q.kids[0].y + q.kids[-1].y)

    def dr(q):
        xq = X(q.age)
        ax.plot([xq, xq], [q.kids[0].y, q.kids[-1].y], color="k", lw=0.8, solid_capstyle="projecting", zorder=4)
        for k in q.kids:
            if k.kids:
                ax.plot([xq, X(k.age)], [k.y, k.y], color="k", lw=0.8, solid_capstyle="butt", zorder=4)
                dr(k)
            else:          # black to the clade's crown, then the pale clade line to the present
                xc = X(k.age)
                if xc - xq > 0.05:
                    ax.plot([xq, xc], [k.y, k.y], color="k", lw=0.8, solid_capstyle="butt", zorder=4)
                ax.plot([xc, x1], [k.y, k.y], color=light(ccol(k.clade)), lw=1.3, solid_capstyle="butt", zorder=4)
    sety(root)
    dr(root)
    return X


YM = {c: yb(r) for r, c in enumerate(PHYLO)}
XT = draw_tree_stil(TB, TX0, TX1, YM)
for r, c in enumerate(PHYLO):
    ax.text(TX1 - 0.4, yb(r) + 0.75, clab(c), ha="right", va="top", fontsize=5.5, color=INK, zorder=10)
ax.plot([XT(TB.age), TX1], [YAX, YAX], color="k", lw=0.6, solid_capstyle="butt", zorder=4)
for t in (800, 600, 400, 200, 0):
    if t <= TB.age:
        ax.plot([XT(t)] * 2, [YAX, YAX + 1.0], color="k", lw=0.5, solid_capstyle="butt", zorder=4)
        ax.text(XT(t), YAX + 1.4, "%d" % t, ha="center", va="top", fontsize=5.0, color=INK, zorder=10)
ax.text(0.5 * (TX0 + TX1), YAX + 4.2, "Time (Ma)", ha="center", va="top", fontsize=6.0, color=INK, zorder=10)
XV0, XV1, VL, VH = 33.5, 74.0, -1.5, 1.5
uv = lambda v: XV0 + (v - VL) / (VH - VL) * (XV1 - XV0)
GRID = np.linspace(VL, VH, 481)
ax.plot([uv(0)] * 2, [YB0 - 0.6 * PB, YAX], color=NULL, lw=0.5, ls=(0, (3, 2)), zorder=1)
for r, c in enumerate(PHYLO):
    d = PAIR[c]; y = yb(r)
    dens = gaussian_kde(d)(GRID)                  # Scott bandwidth on all values; drawn inside the axis window
    hw = dens / dens.max() * 0.40 * PB
    ax.fill_between([uv(g) for g in GRID], y - hw, y + hw, color=light(ccol(c), 0.35), lw=0, alpha=0.8, zorder=2)
    med = float(np.median(d))
    ax.plot([uv(med)] * 2, [y - 0.40 * PB, y + 0.40 * PB], color=ccol(c), lw=1.4, solid_capstyle="butt", zorder=3)
    sh = min(0.80 * PB, 7.5 / aspect(PNG[clab(c)]))
    silhouette_c(ax, PNG[clab(c)], XV1 + 1.6, y, ccol(c), h=sh, z=6)
    q1, q3 = np.percentile(d, [25, 75])
    nb = int(((d < VL) | (d > VH)).sum())
    src.append(("b", "paired_dlog10K", c, "deep_minus_null", "", len(d), "%.6f" % med, "%.6f" % q1, "%.6f" % q3, "",
                "median, quartiles; mean %.6f; %d beyond +-1.5 (violin truncated)" % (float(d.mean()), nb)))
    facts.append(("paired_dlog10K_%s" % c, "n %d median %.3f (IQR %.3f to %.3f) mean %.3f, beyond window %d" % (
        len(d), med, q1, q3, float(d.mean()), nb), "r6_3_fits.tsv, orthogroups with both fits"))
ax.plot([uv(VL), uv(VH)], [YAX, YAX], color="k", lw=0.6, solid_capstyle="butt", zorder=4)
for t in (-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5):
    ax.plot([uv(t)] * 2, [YAX, YAX + 1.0], color="k", lw=0.5, solid_capstyle="butt", zorder=4)
    ax.text(uv(t), YAX + 1.4, sv(t, 1) if t else "0", ha="center", va="top", fontsize=5.0, color=INK, zorder=10)
ax.text(0.5 * (XV0 + XV1), YAX + 4.2, "Paired difference per orthogroup, $\\log_{10}K$ deep $-$ matched null",
        ha="center", va="top", fontsize=6.0, color=INK, zorder=10)
# leader-line labels (stil2b "Mean across species" / "Clade mean"), on the top row
y0r, c0 = yb(0), PHYLO[0]
m0 = float(np.median(PAIR[c0]))
t_ = ax.text(uv(-1.45), y0r - 0.62 * PB, "Orthogroups", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.plot([uv(-1.45) + text_mm(ax, t_)[0] + 0.6, uv(-0.30)], [y0r - 0.62 * PB] * 2, color=INK, lw=0.4, zorder=5)
ax.plot([uv(-0.30), uv(-0.30)], [y0r - 0.62 * PB, y0r - 0.12 * PB], color=INK, lw=0.4, zorder=5)
ax.text(uv(0.62), y0r - 0.62 * PB, "Clade median", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.plot([uv(m0), uv(0.62) - 0.6], [y0r - 0.62 * PB] * 2, color=INK, lw=0.4, zorder=5)
ax.plot([uv(m0), uv(m0)], [y0r - 0.62 * PB, y0r - 0.44 * PB], color=INK, lw=0.4, zorder=5)
letter(fig, Wf, Hf, "b", 1.0, YB0 - 9.4)

# ---- c: berg 1c data columns, relaxed and intensified log odds ratios ----------------------------------------------
YH1, YH2 = YB0 - 6.3, YB0 - 3.4
rows_c = (YB0 - 0.5 * PB, YAX)
ys_c = [yb(r) for r in range(8)] + [YB_ALL]
XN, XS1, XS2 = 86.0, 104.0, 141.5
for r, c in enumerate(PHYLO):
    name_row_c(ax, XN, yb(r), clab(c), XS1, c, None, fs=5.5)
name_row_c(ax, XN, YB_ALL, ALL_LAB, XS1, ALL_LAB, None, fs=5.5)
UB = berg_col(ax, XS1, -1, 3, 6.0, rows_c, YAX, YH1, YH2, head1="No difference", head2=("Null higher", "Deep higher"),
              arrow_to=0.9, xlabel="Relaxed call, log odds ratio")
UC = berg_col(ax, XS2, -4, 1, 5.0, rows_c, YAX, YH1, YH2, head1="No difference", head2=("Null higher", "Deep higher"),
              arrow_to=0.3, xlabel="Intensified call, log odds ratio")
COLS = [ccol(c) for c in PHYLO] + [ccol(ALL_LAB)]
EB, LB, HB = np.r_[lc[:, 0], la[0]], np.r_[lcl[:, 0], lal[0]], np.r_[lch[:, 0], lah[0]]
PB_ = np.r_[lcp[:, 0], lap[0]]
EC, LCC, HCC = np.r_[lc[:, 1], la[1]], np.r_[lcl[:, 1], lal[1]], np.r_[lch[:, 1], lah[1]]
PC_ = np.r_[lcp[:, 1], lap[1]]
berg_bars_c(ax, UB, ys_c, EB, LB, HB, COLS, pitch=PB)
berg_bars_c(ax, UC, ys_c, EC, LCC, HCC, COLS, pitch=PB)
for xcol, E, P_ in ((137.6, EB, PB_), (176.5, EC, PC_)):
    for y, e, p_ in zip(ys_c, E, P_):
        ax.text(xcol, y, "%s%s" % (sv(e, 2), star(p_)), ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "c", 84.0, YB0 - 9.4)

res = qa_v2(fig)
clean = qa_report(res, "Fig2")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig2_qa_v2.txt"))
write_facts("fig2_source_v5.tsv", src)
write_facts("fig2_facts_v5.tsv", facts)
save_v2(fig, "Fig2")
print("Fig2 written; QA clean:", clean)
