"""Fig. 2 v2 - core relaxation effect in the R6-3 sample (production fits; 1,000 deep and 400 matched null branch fits
per clade).

Forms: a, b = Cell 2025 Fig. 4G/H (adjacent edge-less bars, capped error bars, dark condition box, square-patch legend
outside right of the lower panel); c = Sackton et al. 2019 Fig. 4A (dots, grey capped whiskers, labels at 45 degrees)
with the dashed reference line of Sackton Fig. 3C (no label at 0).
Data: R6/R6_3/r6_3_fits.tsv (read-only). Relaxed = P < 0.05 and K < 1, intensified = P < 0.05 and K > 1 (PREREG_r6);
a fit without P or K counts as neither (as the r6g rates). 95% CIs: 2,000 resamples of orthogroups within clade (the
deep and matched null fits of an orthogroup stay together, as r6g/r6h), numpy default_rng(20261031) over PHYLO: the r6
scheme, but a numpy, not an R, random stream. All clades: a, b = mean of the eight clade fractions; c = clade-adjusted
log odds ratio of glm(y ~ deep + clade) (fe_logit; checked against R in _scratch_sty/xcheck_fig2_sty.R).
Outputs (viz_v2_sty/out only): Fig2.pdf/png/tif, fig2_source_v2.tsv, fig2_facts_v2.tsv, Fig2_qa_v2.txt.
"""
from figcommon_v2_sty import *  # noqa: F401,F403

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
    pa = "ab"[o]
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
        src.append(("c", nm, cl, "log_OR", "", "", f6(lc[c, o]), f6(lcl[c, o]), f6(lch[c, o]), "%.4f" % lcp[c, o],
                    "2x2 deep vs matched null; zero-cell resamples %d" % NZ[c, o]))
    src.append(("c", nm, "ALL", "log_OR", "", "", f6(la[o]), f6(lal[o]), f6(lah[o]), "%.4f" % lap[o],
                "glm y ~ deep + clade, all 11200 fits; has_go subset (9208 fits) %.6f" % LG[o]))
write_facts("fig2_source_v2.tsv", src)
pct = lambda v, l, h: "%.1f%% (%.1f-%.1f)" % (v, l, h)
facts = [("key", "value", "source"), ("fits_per_clade", "1,000 deep and 400 matched null (null1) branch fits; "
          "8 clades; 11,200 fits", "r6_3_fits.tsv"),
         ("no_P_or_K", "deep %d, null %d (counted as neither)" % (NOPK["real"], NOPK["null1"]), "r6_3_fits.tsv"),
         ("bootstrap", "2,000 resamples of orthogroups within clade; numpy default_rng(20261031) [new: numpy stream]",
          "fig2_v2_sty.py")]
for o, nm in enumerate(OUT):
    facts.append(("%s_all_clades" % nm, "deep %s vs null %s" % (pct(a1[o], a1l[o], a1h[o]), pct(a0[o], a0l[o],
                  a0h[o])), "fig2_source_v2.tsv [new CIs]"))
    for c, cl in enumerate(PHYLO):
        facts.append(("%s_%s" % (nm, cl), "deep %s vs null %s; log OR %.2f (%.2f to %.2f), %s" % (
            pct(r1[c, o], r1l[c, o], r1h[c, o]), pct(r0[c, o], r0l[c, o], r0h[c, o]), lc[c, o], lcl[c, o],
            lch[c, o], pfmt(lcp[c, o]).replace("$", "")), "fig2_source_v2.tsv [new CIs, P]"))
    facts.append(("%s_pooled_log_OR" % nm, "%.2f (%.2f to %.2f), %s; has_go subset %.4f" % (
        la[o], lal[o], lah[o], pfmt(lap[o]).replace("$", ""), LG[o]), "fig2_source_v2.tsv [new CI, P]"))
Wf, Hf = 180.0, 125.0
X0, AW, AH, CW = 11.0, 137.0, 30.0, 62.0                       # mm: left edge, bar panel width, panel height
TOP = {"a": 5.0, "b": 43.4, "c": 81.9}                         # axes tops; letter baselines 1.9 mm above
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
COND = ("Relaxed selection ($P$ < 0.05, $K$ < 1)", "Intensified selection ($P$ < 0.05, $K$ > 1)")
for o, (pa, lim, tk) in enumerate((("a", 80, 20), ("b", 40, 10))):
    ax = axes_mm(fig, Wf, Hf, X0, TOP[pa], AW, AH)
    v = np.c_[np.r_[r1[:, o], a1[o]], np.r_[r0[:, o], a0[o]]]
    lo = np.c_[np.r_[r1l[:, o], a1l[o]], np.r_[r0l[:, o], a0l[o]]]
    hi = np.c_[np.r_[r1h[:, o], a1h[o]], np.r_[r0h[:, o], a0h[o]]]
    assert hi.max() < lim, ascii((pa, hi.max()))
    grouped_bars(ax, LAB, v, lo, hi, (DEEP, NULL), extra=(8,), lab_fs=5)
    ax.set_ylim(0, lim); ax.set_yticks(np.arange(0, lim + 1, tk))
    ax.set_ylabel("Orthogroups (%)", fontsize=5.5)
    cond_box(ax, COND[o])
    letter(fig, Wf, Hf, pa, 1.0, TOP[pa] - 1.9)
sq_legend(fig, (X0 + AW + 3) / Wf, 1 - (TOP["b"] + AH / 2) / Hf, ["Deep branches", "Matched null branches"],
          (DEEP, NULL))
letter(fig, Wf, Hf, "c", 1.0, TOP["c"] - 1.9)
xs = np.r_[np.arange(C), C + .8]
for o, x0 in enumerate((X0, X0 + AW - CW)):
    ax = axes_mm(fig, Wf, Hf, x0, TOP["c"], CW, AH)
    e, l, h = np.r_[lc[:, o], la[o]], np.r_[lcl[:, o], lal[o]], np.r_[lch[:, o], lah[o]]
    ref_line(ax, 0)
    dot_whisker(ax, xs, e, l, h, LAB, lab_fs=5)
    y0, y1 = np.floor(l.min() - .05), np.ceil(h.max() + .05)
    ax.set_ylim(y0, y1); ax.set_yticks(np.arange(y0, y1 + .5, 1))
    ax.set_ylabel("Log odds ratio (%s)" % OUT[o], fontsize=5.5)
res = qa_v2(fig)
clean = qa_report(res, "Fig2")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig2_qa_v2.txt"))
write_facts("fig2_facts_v2.tsv", facts)
save_v2(fig, "Fig2")
print("Fig2 written; QA clean:", clean)
