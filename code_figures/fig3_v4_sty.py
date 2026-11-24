"""Fig. 3 v4 (= v3 Fig. 4 plus the blind replicate) -- every panel copies the form of one viewed published panel
(paper/figure_benchmark_map.md, "Fig. 2-4 v4"). v4 additions: d carries the blind-replicate rows (R6-2b) under the
first-sample rows; e has two data columns on one tree and one row set (first sample | blind replicate), the second
column without a tree as in berg 1c's shared rows (v3 a/b).

Forms: a, b = Bergeron et al. 2023 Nature Fig. 1c (time tree, upright names, dotted leaders, silhouettes, one data
column with a thin reference line, a two-line header and one double-headed arrow; b shares a's rows and carries no
tree); c = Wei et al. 2026 Nature Extended Data histogram (terr_ED8.jpg), two stacked panels on one x axis, null bars
light grey; d = Berg 1c data column only, no tree; e = same paper Fig. 1a (berg1a_top.png, berg1_full.png): time tree
with an age axis, alternating stripes one tick step wide, one dot per lineage, black-edged mean square with a capped
interval, a count box per group at the right.
Data (read-only): R6/R6_3/r6_3_stats.tsv ca2600c9 (a, b); R6/R6_2/r6_2_summary.tsv acb54400 (d, e);
R6/R6_2/r6_2_lineages.tsv 997e9842 (c top, e rows); R6/R6_2/r6_2_fits.tsv 97078a7f (c bottom, descriptive only).
No bootstrap is redone; every number equals the v2 script's.
"""
from figcommon_v4_sty import *  # noqa: F401,F403
import re, hashlib
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
from collections import defaultdict
mt_arial()

S3 = {(r["analysis"], r["metric"], r["clade"]): r for r in rtsv(W("addon", "R6", "R6_3", "r6_3_stats.tsv"))}
g3 = lambda m, c: tuple(fnum(S3[("R6-3a_paired", m, c)][k]) for k in ("estimate", "lo", "hi"))
assert abs(g3("Dpur", "ALL")[0] - .077440575607305) < 1e-12 and S3[("R6-3a_paired", "Dpur", "ALL")]["n_clades_pos"] == "6"
SM = {r["stat"]: tuple(fnum(r[k]) for k in ("est", "lo", "hi")) for r in rtsv(W("addon", "R6", "R6_2", "r6_2_summary.tsv"))}
assert abs(SM["pooled"][0] + .0512255539301645) < 1e-12 and abs(SM["share_lt0"][0] - .707692307692308) < 1e-12
LIN = rtsv(W("addon", "R6", "R6_2", "r6_2_lineages.tsv"))
assert len(LIN) == 65
lk = lambda K: np.log10(min(max(K, .001), 50))
val = defaultdict(list)
for r in rtsv(W("addon", "R6", "R6_2", "r6_2_fits.tsv")):
    if r["fit"] == "ok" and r["canon"] not in ("merged", "", "NA") and np.isfinite(fnum(r["K"])):
        val[(r["clade"], r["set"], r["canon"])].append(lk(fnum(r["K"])))
NMED = {c: np.median([x for (cc, s, _), v in val.items() if cc == c and s == "null1" for x in v]) for c in {k[0] for k in val}}
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
# blind replicate R6-2b (re-drawn orthogroups; read-only probe of r6_2_report.txt, md5 b894c4bf)
_PB = W("addon", "R6", "probe_r6_2b.txt")
assert md5(_PB)[:8] == "b894c4bf", md5(_PB)
_L = open(_PB, encoding="utf-8").read().split("\n")
_i = _L.index(next(l for l in _L if l.strip().startswith("stat") and "est" in l))
SB = {}
for l in _L[_i + 1:]:
    m = re.match(r"\s*\d+\s+(\S+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$", l)
    if not m:
        break
    SB[m.group(1)] = tuple(float(m.group(k)) for k in (2, 3, 4))
_i = _L.index(next(l for l in _L if l.startswith("per-lineage effects")))
LINB = []
for l in _L[_i + 2:]:
    m = re.match(r"\s*\d+\s+(\S+)\s+(\S+)\s+(-?[\d.]+)\s+(\d+)\s+(\d+)\s*$", l)
    if not m:
        break
    LINB.append((m.group(1), m.group(2), float(m.group(3))))
assert len(LINB) == 81 and len(SB) == 21 and abs(SB["pooled"][0] + 0.0715) < 1e-9, (len(LINB), len(SB))
BYB = {c: [e for cc, _, e in LINB if cc == c] for c in PHYLO}
for c in PHYLO:   # clade means in the report are means of the printed (5-decimal) lineage effects
    assert abs(np.mean(BYB[c]) - SB["clade_" + c][0]) < 6e-5, (c, np.mean(BYB[c]), SB["clade_" + c][0])
REALB = np.array([e for _, _, e in LINB])

Wf, Hf = 180.0, 210.0
sv = lambda v, d=2: (MINUS if v < 0 else "") + ("%." + str(d) + "f") % abs(v)
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
ax = mm_axes(fig, Wf, Hf)
f6 = lambda v: "%.6f" % v
src = [("panel", "group", "series", "estimate", "lo95", "hi95", "n", "note")]
facts = [("key", "value", "source")]


def eff_col(x0, lo, hi, unit, ytop, ybot, ticks, yaxis, xlabel, spine=True, ref=None):
    """berg 1a data column: alternating stripes one tick step wide from the axis minimum (grey first), thin left
    spine, baseline with outward ticks and true-minus labels, centred x label. Returns value -> x.
    berg 1a's black vertical line is its spine, at its own reference 0; here 0 falls inside the axis, so the same
    line is drawn at `ref` (berg 1c's form) instead."""
    u = lambda v: x0 + (v - lo) * unit
    step = ticks[1] - ticks[0]
    t, on = lo, True
    while t < hi - 1e-9:
        t2 = min(t + step, hi)
        if on:
            ax.fill_betweenx([ybot, ytop], u(t), u(t2), color=STRIPE, lw=0, zorder=0.3)
        t, on = t2, not on
    if spine:
        ax.plot([u(lo), u(lo)], [yaxis, ytop], color="k", lw=0.7, solid_capstyle="butt", zorder=1.6)
    if ref is not None:
        ax.plot([u(ref), u(ref)], [yaxis, ytop], color="k", lw=0.5, solid_capstyle="butt", zorder=1.7)
    ax.plot([u(lo), u(hi)], [yaxis, yaxis], color="k", lw=0.6, solid_capstyle="butt", zorder=1.6)
    for tv in ticks:
        ax.plot([u(tv), u(tv)], [yaxis, yaxis + 1.0], color="k", lw=0.5, solid_capstyle="butt", zorder=1.6)
        ax.text(u(tv), yaxis + 1.4, (MINUS + "%g" % -tv) if tv < 0 else "%g" % tv, ha="center", va="top",
                fontsize=5.0, color=INK, zorder=10)
    ax.text(0.5 * (u(lo) + u(hi)), yaxis + 1.4 + 5.0 * 0.46, xlabel, ha="center", va="top", fontsize=6.0, color=INK,
            zorder=10)
    return u


# ---- a, b: Berg 1c, one shared row set; the tree is drawn once, in a ------------------------------------------------
PA = 5.8
YA0 = 18.0
yb = lambda r: YA0 + r * PA
YB_ALL = yb(7) + 0.55 * PA
YAX_A = YB_ALL + 0.5 * PA
rows_ab = (YA0 - 0.5 * PA, YAX_A)
YH1, YH2 = YA0 - 4.0, YA0 - 7.3
TB = focal_tree(PHYLO)
YB = {**{c: yb(r) for r, c in enumerate(PHYLO)}, **{clab(c): yb(r) for r, c in enumerate(PHYLO)}}
draw_tree(ax, TB, 3.0, 10.5, YB)
SIL = {clab(c): min(0.9 * PA, 6.0 / aspect(PNG[clab(c)])) for c in PHYLO}
XS_A = 40.0
XVA = 86.0
COLS_AB = [ccol(c) for c in PHYLO] + [ccol(ALL_LAB)]
for r, c in enumerate(PHYLO):
    name_row_c(ax, 11.5, yb(r), "%s (%d)" % (clab(c), int(S3[("R6-3a_paired", "Dpur", c)]["n"])), XS_A, clab(c),
               PNG[clab(c)], h=SIL[clab(c)])
name_row_c(ax, 11.5, YB_ALL, "%s (%d)" % (ALL_LAB, int(S3[("R6-3a_paired", "Dpur", "ALL")]["n"])), XS_A, ALL_LAB)
UA = berg_col(ax, XS_A, -0.3, 0.3, 68.0, rows_ab, YAX_A, YH1, YH2, head1="No difference",
              head2=("Null higher", "Deep higher"), arrow_to=0.12,
              xlabel="$\\Delta_{\\mathrm{pur}}$, vs matched null", ticks=[-0.3, -0.15, 0, 0.15, 0.3])
VA = np.array([g3("Dpur", c) for c in PHYLO] + [g3("Dpur", "ALL")])
assert VA[:, 1].min() >= -0.3 and VA[:, 2].max() <= 0.3, VA[:, [1, 2]]
ys_ab = [yb(r) for r in range(8)] + [YB_ALL]
berg_bars_c(ax, UA, ys_ab, VA[:, 0], VA[:, 1], VA[:, 2], COLS_AB, pitch=PA)
for y, row in zip(ys_ab, VA):
    ax.text(XVA, y, "%+.2f" % row[0], ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "a", 1.0, 6.0)
for c, row in zip(PHYLO + ["ALL"], VA):
    src.append(("a", c, "Dpur", f6(row[0]), f6(row[1]), f6(row[2]), S3[("R6-3a_paired", "Dpur", c)]["n"],
                "r6_3_stats.tsv; p_boot %s" % S3[("R6-3a_paired", "Dpur", c)]["p_boot"]))
facts.append(("R6_3a_Dpur", "%.3f (%.3f to %.3f), p_boot %s, %s of 8 clades > 0" % (
    g3("Dpur", "ALL") + (S3[("R6-3a_paired", "Dpur", "ALL")]["p_boot"],
                         S3[("R6-3a_paired", "Dpur", "ALL")]["n_clades_pos"])), "r6_3_stats.tsv"))

XS_B = 126.0
XVB = 166.0
UB = berg_col(ax, XS_B, -0.06, 0.06, 300.0, rows_ab, YAX_A, YH1, YH2, head1="No difference",
              head2=("Null higher", "Deep higher"), arrow_to=0.024,
              xlabel="$\\Delta_{\\mathrm{pos}}$, vs matched null", ticks=[-0.06, -0.03, 0, 0.03, 0.06])
VB = np.array([g3("Dpos", c) for c in PHYLO] + [g3("Dpos", "ALL")])
assert VB[:, 1].min() >= -0.06 and VB[:, 2].max() <= 0.06, VB[:, [1, 2]]
berg_bars_c(ax, UB, ys_ab, VB[:, 0], VB[:, 1], VB[:, 2], COLS_AB, pitch=PA)
for y, row in zip(ys_ab, VB):
    ax.text(XVB, y, "%+.3f" % row[0], ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "b", 104.0, 6.0)
for c, row in zip(PHYLO + ["ALL"], VB):
    src.append(("b", c, "Dpos", f6(row[0]), f6(row[1]), f6(row[2]), S3[("R6-3a_paired", "Dpos", c)]["n"],
                "r6_3_stats.tsv; p_boot %s" % S3[("R6-3a_paired", "Dpos", c)]["p_boot"]))
facts.append(("R6_3a_Dpos", "%.3f (%.3f to %.3f), p_boot %s, %s of 8 clades > 0" % (
    g3("Dpos", "ALL") + (S3[("R6-3a_paired", "Dpos", "ALL")]["p_boot"],
                         S3[("R6-3a_paired", "Dpos", "ALL")]["n_clades_pos"])), "r6_3_stats.tsv"))
facts.append(("rows_ab", "All clades = equal-weighted mean of the 8 clade estimates, equals the paired ALL estimate",
              "r6_3_stats.tsv"))

# ---- c: Wei et al. 2026 ED histogram, two stacked panels on one x axis, null bars light grey ------------------------
T2, hh, WC = 76.0, 24.0, 60.0
bins = np.arange(-1.1, 1.11, .1)
assert REAL.min() > bins[0] and REAL.max() < bins[-1] and NULLE.min() > bins[0] and NULLE.max() < bins[-1], \
    (REAL.min(), REAL.max(), NULLE.min(), NULLE.max())
CB = ((REAL, None, SM["pooled"][0],
       ["Pooled %s" % ("%.3f" % SM["pooled"][0]).replace("-", "−"), "%.0f%% < 0" % (100 * SM["share_lt0"][0])],
       "Deep lineages (%d)" % len(REAL)),
      (NULLE, NULL, 0.0, ["Null", "%.0f%% < 0" % (100 * (NULLE < 0).mean())],
       "Null pseudo-lineages (%d)" % len(NULLE)))
for j, (v, col, at, lines, nm) in enumerate(CB):
    axh = axes_mm(fig, Wf, Hf, 12.0, T2 + j * (hh + 3.0), WC, hh)
    ed8_hist(axh, v, bins, at, lines, xlabel="Lineage effect, Δ median $\\log_{10}$ $K$" if j else None,
             ylabel="Lineages", color=col, head=1.6)
    axh.set_xlim(bins[0], bins[-1])
    if not j:
        axh.set_xticklabels([])
    cond_box(axh, nm, fs=5)
letter(fig, Wf, Hf, "c", 1.0, 72.0)
for r in LIN:
    src.append(("c", "%s|%s" % (r["clade"], r["canon"]), "deep_lineage", f6(fnum(r["eff"])), "", "", r["n_hog"],
                "r6_2_lineages.tsv"))
for (c, s, cn), (e, n) in sorted(EFF.items()):
    if s == "null1":
        src.append(("c", "%s|%s" % (c, cn), "null_pseudo_lineage", f6(e), "", "", n,
                    "r6_2_fits.tsv, descriptive"))
facts.append(("lineage_effects", "deep %d lineages, %.1f%% < 0 (%.3f to %.3f); null %d pseudo-lineages, %.1f%% < 0" % (
    len(REAL), 100 * (REAL < 0).mean(), REAL.min(), REAL.max(), len(NULLE), 100 * (NULLE < 0).mean()),
    "r6_2_lineages.tsv; r6_2_fits.tsv [null descriptive]"))

# ---- d: Berg 1c data column only, the robustness settings (no tree) --------------------------------------------------
DS = [("pooled", "Pooled"), ("lineage_weighted", "Lineage weighted"), ("without_top2_clades", "−Top 2 clades"),
      ("without_largest_lineage", "−Largest lineage")] + [("loco_" + c, "−" + clab(c)) for c in PHYLO
                                                          if "loco_" + c in SM]
DSB = [("pooled", "Blind: pooled"), ("without_top2_clades", "Blind: −Top 2 clades"),
       ("without_largest_lineage", "Blind: −Largest lineage")]
PD = 4.0
YD0 = 82.0
NDA = len(DS)
yd = lambda r: YD0 + r * PD + (1.0 if r >= NDA else 0.0)    # small gap before the blind rows
NDT = NDA + len(DSB)
YAX_D = yd(NDT - 1) + 0.5 * PD
XS_D = 132.0
XVD = 170.5
UD = berg_col(ax, XS_D, -0.185, 0.025, 155.0, (YD0 - 0.5 * PD, YAX_D), YAX_D, YD0 - 4.0, YD0 - 7.3,
              head1="No effect", head2=("Lower", "Higher"), arrow_to=0.012,
              xlabel="Pooled lineage effect\nΔ median $\\log_{10}$ $K$", ticks=[-0.15, -0.1, -0.05, 0.0])
VD = np.array([SM[k] for k, _ in DS] + [SB[k] for k, _ in DSB])
assert VD[:, 1].min() >= -0.185 and VD[:, 2].max() <= 0.025, VD[:, [1, 2]]
COLS_D = [DEEP if not k.startswith("loco_") else ccol(k[5:]) for k, _ in DS] + [DEEP2] * len(DSB)
for r, ((k, lab), col) in enumerate(zip(DS + DSB, COLS_D)):
    name_row(ax, 98.0, yd(r), lab, XS_D, color=col)
berg_bars_c(ax, UD, [yd(r) for r in range(NDT)], VD[:, 0], VD[:, 1], VD[:, 2], COLS_D, pitch=PD)
for r, row in enumerate(VD):
    ax.text(XVD, yd(r), "%+.3f" % row[0], ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "d", 92.0, 72.0)
for (k, _), row in zip(DS, VD):
    src.append(("d", k, "R6-2", f6(row[0]), f6(row[1]), f6(row[2]), "", "r6_2_summary.tsv"))
for (k, _), row in zip(DSB, VD[NDA:]):
    src.append(("d", k, "R6-2b_blind", f6(row[0]), f6(row[1]), f6(row[2]), "", "probe_r6_2b.txt"))
facts.append(("R6_2_claim", "pooled %.4f (%.4f to %.4f); share<0 %.3f (%.3f to %.3f); without top2 %.4f (%.4f to %.4f)"
              "; SUPPORTED" % (SM["pooled"] + SM["share_lt0"] + SM["without_top2_clades"]), "r6_2_report.txt"))

# ---- e: Berg 1a: age tree, stripes, one dot per lineage, mean square; two data columns on one row set --------------
TE, PE = 154.0, 5.6
EC = list(PHYLO)
ye = lambda r: TE + r * PE
YAX_E = ye(len(EC) - 1) + PE
TET = focal_tree(EC)
YEM = {**{c: ye(r) for r, c in enumerate(EC)}, **{clab(c): ye(r) for r, c in enumerate(EC)}}
XE = draw_tree(ax, TET, 4.0, 22.0, YEM)
XS_E1, XS_E2, ELO, EHI, EU = 52.0, 115.0, -0.75, 1.05, 30.0
for r, c in enumerate(EC):
    name_row_c(ax, 23.0, ye(r), clab(c), XS_E1, clab(c), PNG[clab(c)], h=min(0.9 * PE, 5.0 / aspect(PNG[clab(c)])))
ax.plot([4.0, 22.0], [YAX_E, YAX_E], color="k", lw=0.6, solid_capstyle="butt", zorder=1.6)
for a in (600, 400, 200, 0):
    xa = XE(a)
    ax.plot([xa, xa], [YAX_E, YAX_E + 1.0], color="k", lw=0.5, solid_capstyle="butt", zorder=1.6)
    ax.text(xa, YAX_E + 1.4, "%d" % a, ha="center", va="top", fontsize=5.0, color=INK, zorder=10)
ax.text(13.0, YAX_E + 1.4 + 5.0 * 0.46, "Age (million years)", ha="center", va="top", fontsize=6.0, color=INK, zorder=10)
ETK = [-0.5, 0.0, 0.5, 1.0]
BYC = {c: [fnum(r["eff"]) for r in LIN if r["clade"] == c] for c in EC}
COLSPEC = ((XS_E1, BYC, lambda c: SM.get("clade_" + c), "First sample (%d lineages)" % len(REAL), "R6-2"),
           (XS_E2, BYB, lambda c: SB.get("clade_" + c), "Blind replicate (%d lineages)" % len(REALB), "R6-2b"))
for xs, BY, MS, head, tag in COLSPEC:
    UE = eff_col(xs, ELO, EHI, EU, TE - 0.5 * PE, YAX_E, ETK, YAX_E, "Δ median $\\log_{10}$ $K$ per lineage", ref=0.0)
    ax.text(UE(0.5 * (ELO + EHI)), TE - 0.5 * PE - 1.6, head, ha="center", va="baseline", fontsize=6.0, color=INK,
            zorder=10)
    for r, c in enumerate(EC):
        y, cc, v = ye(r), ccol(c), BY[c]
        if not v:
            continue
        vin = [x for x in v if x >= ELO]
        assert max(v) <= EHI, (tag, c, max(v))
        ax.scatter([UE(x) for x in vin], [y] * len(vin), s=1.8, c=cc, lw=0, alpha=.85, zorder=3)
        for x in v:
            if x < ELO:     # off-axis lineage: open triangle at the axis edge, value beside it
                ax.scatter([UE(ELO) + 0.9], [y], s=5.0, marker="<", facecolor="none", edgecolor=cc, linewidths=0.5,
                           zorder=4)
                ax.text(UE(ELO) + 2.1, y - 1.6, sv(x, 2), ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
        m = MS(c)
        k = sum(1 for x in v if x < 0)
        if m is not None and len(v) >= 2:
            assert abs(float(np.mean(v)) - m[0]) < 6e-5, (tag, c)
            h = 0.45 * PE
            ax.plot([UE(m[1]), UE(m[2])], [y, y], color=CI, lw=0.6, solid_capstyle="butt", zorder=3)
            for xx in (UE(m[1]), UE(m[2])):
                ax.plot([xx, xx], [y - 0.3 * h, y + 0.3 * h], color=CI, lw=0.6, solid_capstyle="butt", zorder=3)
            ax.scatter([UE(m[0])], [y], s=7.0, marker="s", facecolor=cc, edgecolor="k", linewidths=0.4, zorder=5)
            src.append(("e", c, tag + "_clade_mean", f6(m[0]), f6(m[1]), f6(m[2]), "%d/%d" % (k, len(v)),
                        "lineages <0 / summarised"))
        ax.text(UE(EHI) + 0.8, y, "%d/%d" % (k, len(v)), ha="left", va="center", fontsize=5.0, color=INK,
                zorder=10)
    facts.append(("%s_lineages" % tag, "%d lineages, %d < 0; pooled %.4f" % (
        sum(len(v) for v in BY.values()), sum(1 for v in BY.values() for x in v if x < 0),
        (SM if tag == "R6-2" else SB)["pooled"][0]), "r6_2_summary.tsv" if tag == "R6-2" else "probe_r6_2b.txt"))
letter(fig, Wf, Hf, "e", 1.0, 145.0)

res = qa_v2(fig)
clean = qa_report(res, "Fig3")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig3_qa_v2.txt"))
write_facts("fig3_source_v4.tsv", src)
write_facts("fig3_facts_v4.tsv", facts)
save_v2(fig, "Fig3")
print("Fig3 written; QA clean:", clean)
