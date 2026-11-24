"""Fig. 2 v3 - core relaxation effect in the R6-3 sample (production fits; 1,000 deep and 400 matched null branch fits
per clade). Forms copied from two viewed published panels (panel map: paper/figure_benchmark_map.md, "Fig. 2-4 v3"):
  a     = Xu et al. 2025 Cell Fig. 2B: time tree, upright names, two thin-framed stacked-bar panels (deep, matched
          null) with the stacked box legend inside each frame
  b, c  = Bergeron et al. 2023 Nature Fig. 1c: time tree + upright names + dotted leaders + silhouettes, black bars
          from the no-difference value with grey capped 95% intervals, half-step grey stripes, thin reference line
          and a two-line header with one double-headed arrow
Data: R6/R6_3/r6_3_fits.tsv (read-only). Relaxed = P < 0.05 and K < 1, intensified = P < 0.05 and K > 1 (PREREG_r6);
a fit without P or K counts as neither (as the r6g rates). 95% CIs: 2,000 resamples of orthogroups within clade (the
deep and matched null fits of an orthogroup stay together, as r6g/r6h), numpy default_rng(20261031) over PHYLO: the r6
scheme, but a numpy, not an R, random stream. All clades: a = mean of the eight clade percentages; b, c = clade-adjusted
log odds ratio of glm(y ~ deep + clade) (fe_logit; checked against R in _scratch_sty/xcheck_fig2_sty.R).
Outputs (viz_v3_sty/out only): Fig2.pdf/png/tif, fig2_source_v3.tsv, fig2_facts_v3.tsv, Fig2_qa_v2.txt.
"""
from figcommon_v3_sty import *  # noqa: F401,F403

mt_arial()
rows = rtsv(W("addon", "R6", "R6_3", "r6_3_fits.tsv"))
assert len(rows) == 11200, len(rows)
FIT = {}
for r in rows:
    p, K = fnum(r["relax_p"]), fnum(r["relax_K"]); ok = bool(np.isfinite(p) and np.isfinite(K))
    k = (r["clade"], r["label"], r["hog"])
    assert k not in FIT and r["label"] in ("real", "null1"), ascii(k)
    FIT[k] = (int(ok and p < .05 and K < 1), int(ok and p < .05 and K > 1), fnum(r["has_go"]) == 1, ok)
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
write_facts("fig2_source_v3.tsv", src)
pct = lambda v, l, h: "%.1f%% (%.1f-%.1f)" % (v, l, h)
facts = [("key", "value", "source"), ("fits_per_clade", "1,000 deep and 400 matched null (null1) branch fits; "
          "8 clades; 11,200 fits", "r6_3_fits.tsv"),
         ("no_P_or_K", "deep %d, null %d (counted as neither)" % (NOPK["real"], NOPK["null1"]), "r6_3_fits.tsv"),
         ("bootstrap", "2,000 resamples of orthogroups within clade; numpy default_rng(20261031) [new: numpy stream]",
          "fig2_v3_sty.py")]
for o, nm in enumerate(OUT):
    facts.append(("%s_all_clades" % nm, "deep %s vs null %s" % (pct(a1[o], a1l[o], a1h[o]), pct(a0[o], a0l[o],
                  a0h[o])), "fig2_source_v3.tsv [new CIs]"))
    for c, cl in enumerate(PHYLO):
        facts.append(("%s_%s" % (nm, cl), "deep %s vs null %s; log OR %.2f (%.2f to %.2f), %s" % (
            pct(r1[c, o], r1l[c, o], r1h[c, o]), pct(r0[c, o], r0l[c, o], r0h[c, o]), lc[c, o], lcl[c, o],
            lch[c, o], pfmt(lcp[c, o]).replace("$", "")), "fig2_source_v3.tsv [new CIs, P]"))
    facts.append(("%s_pooled_log_OR" % nm, "%.2f (%.2f to %.2f), %s; has_go subset %.4f" % (
        la[o], lal[o], lah[o], pfmt(lap[o]).replace("$", ""), LG[o]), "fig2_source_v3.tsv [new CI, P]"))

def txt_on(hexcol):
    """White on dark fills, ink on light fills (Cell 2B legend boxes)."""
    r, g, b = (int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5))
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "w" if lum < 0.45 else INK


# ---------------------------------------------------------------------------------------------------------------
Wf, Hf = 180.0, 128.0
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
ax = mm_axes(fig, Wf, Hf)


def sv(v, d=1):
    """Value with a true minus sign (U+2212, not a hyphen) at `d` decimals."""
    return (MINUS if v < 0 else "") + ("%." + str(d) + "f") % abs(v)


def star(p):
    """Bootstrap significance keyed to the manuscript: *** P<0.01, ** P<0.05, * P<0.1."""
    return "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.1 else ""


def chip_legend(xc, y, entries, h=2.4, gap=1.2, fs=5.2):
    """One legend row centred on xc: a filled chip plus its label per entry (Cell 2B legend boxes, lifted out of
    the data area so that no bar runs underneath them)."""
    ts = [ax.text(0, 0, lab, fontsize=fs, ha="left", va="center") for _, lab in entries]
    ws = [h + 0.9 + text_mm(ax, t)[0] for t in ts]
    for t in ts:
        t.remove()
    x = xc - (sum(ws) + gap * (len(ws) - 1)) / 2
    for (fill, lab), w in zip(entries, ws):
        ax.add_patch(Rectangle((x, y - h / 2), h, h, fc=fill, ec="none", lw=0, zorder=6))
        ax.text(x + h + 0.9, y, lab, ha="left", va="center", fontsize=fs, color=INK, zorder=10)
        x += w + gap


# ---- a: Cell 2B, two framed stacked-bar panels off one tree -----------------------------------------------------
YTA, PA, GAP = 11.0, 4.6, 0.6                     # frame top, row pitch, row gap before All clades
yc = lambda r: YTA + (r + 0.5) * PA
YALL = YTA + (8 + 0.5 + GAP) * PA
YBA = YALL + 0.5 * PA
ya_ax = [yc(r) - YTA for r in range(8)] + [YALL - YTA]
ymap = lambda yf: {**{c: yf(r) for r, c in enumerate(PHYLO)}, **{clab(c): yf(r) for r, c in enumerate(PHYLO)}}
            # focal_tree keeps the raw clade code (Actinopteri_m); names use the label (Actinopteri)
draw_tree(ax, focal_tree(PHYLO), 2.0, 18.0, ymap(yc))
for r, c in enumerate(PHYLO):
    chip_row(ax, 19.0, yc(r), clab(c), c, w=1.5, h=1.5, gap=0.9, fs=6.0)
chip_row(ax, 19.0, YALL, ALL_LAB, ALL_LAB, w=1.5, h=1.5, gap=0.9, fs=6.0)
FR = (40.0, 67.5), (108.3, 67.7)                  # left edge and width of the deep and matched-null frames
axd = cell2b_frame(fig, Wf, Hf, FR[0][0], YTA, FR[0][1], YBA - YTA, (-2, 82), (0, 20, 40, 60, 80), "Deep branches")
axn = cell2b_frame(fig, Wf, Hf, FR[1][0], YTA, FR[1][1], YBA - YTA, (-2, 82), (0, 20, 40, 60, 80),
                   "Matched null branches")
DREL, DINT = np.r_[r1[:, 0], a1[0]], np.r_[r1[:, 1], a1[1]]
NREL, NINT = np.r_[r0[:, 0], a0[0]], np.r_[r0[:, 1], a0[1]]
for axf, rel, itn, c1, c2 in ((axd, DREL, DINT, DEEP, DEEP2), (axn, NREL, NINT, NULL, LIGHT)):
    axf.barh(ya_ax, rel, height=0.8 * PA, color=c1, lw=0, zorder=2)
    axf.barh(ya_ax, itn, left=rel, height=0.8 * PA, color=c2, lw=0, zorder=2)
    for y, a, b in zip(ya_ax, rel, itn):          # bar-end total, i.e. orthogroups with a significant RELAX test
        end_val(axf, lambda v: v, y, "%.1f" % (a + b), a + b, base=0.0, fs=5.0, gap=1.4)
chip_legend(0.5 * (FR[0][0] + FR[0][1]), 4.4, ((DEEP, "Relaxed"), (DEEP2, "Intensified")))
chip_legend(0.5 * (FR[1][0] + FR[1][1]), 4.4, ((NULL, "Relaxed"), (LIGHT, "Intensified")))
ax.text(0.5 * (FR[0][0] + FR[1][0] + FR[1][1]), YBA + 3.4, "Orthogroups with a significant RELAX test (%)",
        ha="center", va="top", fontsize=6, color=INK, zorder=10)
letter(fig, Wf, Hf, "a", 1.0, 3.6)

# ---- b, c: Berg 1c, same tree and name column, one column per outcome -------------------------------------------
YH1, YH2, PB, YB0 = 63.6, 66.8, 5.8, 70.2
yb = lambda r: YB0 + r * PB
YB_ALL = yb(7) + 0.6 * PB
YAX = YB_ALL + 0.5 * PB
rows_bc = (YB0 - 0.5 * PB, YAX)
ys_b = [yb(r) for r in range(8)] + [YB_ALL]
TB = focal_tree(PHYLO); YB = ymap(yb)
SIL = {clab(c): min(0.9 * PB, 6.0 / aspect(PNG[clab(c)])) for c in PHYLO}
COLS = [ccol(c) for c in PHYLO] + [ccol(ALL_LAB)]
for x0, x1, xname, xsp in ((2.0, 16.0, 17.0, 42.0), (82.0, 96.0, 97.0, 122.0)):
    draw_tree(ax, TB, x0, x1, YB)
    for r, c in enumerate(PHYLO):
        name_row_c(ax, xname, yb(r), clab(c), xsp, c, PNG[clab(c)], h=SIL[clab(c)])
    name_row_c(ax, xname, YB_ALL, ALL_LAB, xsp, ALL_LAB, None)
UB = berg_col(ax, 42.0, -1, 3, 8.5, rows_bc, YAX, YH1, YH2, head1="No difference",
              head2=("Null higher", "Deep higher"), xlabel="Relaxed call, log odds ratio")
UC = berg_col(ax, 122.0, -4, 1, 8.5, rows_bc, YAX, YH1, YH2, head1="No difference",
              head2=("Null higher", "Deep higher"), xlabel="Intensified call, log odds ratio")
EB, LB, HB = np.r_[lc[:, 0], la[0]], np.r_[lcl[:, 0], lal[0]], np.r_[lch[:, 0], lah[0]]
PB_ = np.r_[lcp[:, 0], lap[0]]
EC, LCC, HCC = np.r_[lc[:, 1], la[1]], np.r_[lcl[:, 1], lal[1]], np.r_[lch[:, 1], lah[1]]
PC_ = np.r_[lcp[:, 1], lap[1]]
berg_bars_c(ax, UB, ys_b, EB, LB, HB, COLS, pitch=PB)
berg_bars_c(ax, UC, ys_b, EC, LCC, HCC, COLS, pitch=PB)
XBCOL, XCCOL = 79.5, 171.5                            # value column at the frame end (berg 1a count boxes)
for xcol, ysc, E, P in ((XBCOL, ys_b, EB, PB_), (XCCOL, ys_b, EC, PC_)):
    for y, e, p in zip(ysc, E, P):
        ax.text(xcol, y, "%.2f%s" % (e, star(p)), ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(XBCOL, rows_bc[0] - 0.75 * PB, "log OR", ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(XCCOL, rows_bc[0] - 0.75 * PB, "log OR", ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "b", 1.0, 59.0); letter(fig, Wf, Hf, "c", 80.0, 59.0)

res = qa_v2(fig)
clean = qa_report(res, "Fig2")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig2_qa_v2.txt"))
write_facts("fig2_facts_v3.tsv", facts)
save_v2(fig, "Fig2")
print("Fig2 written; QA clean:", clean)
