"""Fig. 3 v5 - weaker purifying selection and recurrence across deep lineages, redrawn in denser published forms
(panel map: paper/figure_benchmark_map.md, "Fig. 2-4 v5"); letters keep the v4 contents:
  a, b = stil5a (Stiller et al. 2024 Nature Fig. 5a): per clade, filled box plots of the per-orthogroup change in
         the BUSTED-PH omega classes from background to foreground on deep (blue) and matched null (grey) branches
         (400 paired orthogroups per clade), black outliers, a row of stars or exact bootstrap P of the paired
         excess over each clade, framed legend with the box key. a = Delta_pur, b = Delta_pos.
  c    = spid2b (Tong et al. 2022 Nat. Commun. Fig. 2b): framed horizontal violins with every lineage as a black
         point and a black median point on an interquartile line: deep lineages of the first and the blind sample
         and the null pseudo-lineages of the first sample.
  d    = berg 1c data column (v4 d): pooled lineage effect under the robustness settings, first and blind samples.
  e    = cell3 (Xu et al. 2025 Cell Fig. 3): tree of the 81 deep lineages (LSD2 time tree pruned to one species per
         lineage) with clade labels and node diamonds, depth-class squares, names, value-sized circles with numbers
         (species, orthogroups per sample) and bar columns (lineage effect per sample) in the depth-class colour,
         grey bands for alternate clades, pill headers and a depth legend.
Data (read-only, as v4): R6/R6_3/r6_3_fits.tsv, r6_3_stats.tsv; R6/R6_2/r6_2_summary.tsv, r6_2_lineages.tsv,
r6_2_fits.tsv, canon.tsv; R6/probe_r6_2b.txt (md5 asserted); meta/labels_v2.tsv (species, depth), chemo_habitat.tsv.
"""
from figcommon_v5_sty import *  # noqa: F401,F403
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

from scipy.stats import gaussian_kde
from matplotlib.patches import FancyBboxPatch, Circle, Polygon
sv = lambda v, d=2: (MINUS if v < 0 else "") + ("%." + str(d) + "f") % abs(v)
f6 = lambda v: "%.6f" % v
src = [("panel", "group", "series", "estimate", "lo95", "hi95", "n", "note")]
facts = [("key", "value", "source")]

# ---- per-orthogroup omega-class changes on the 400 paired orthogroups per clade (R6-3a) ---------------------------
OM = {}
for r in rtsv(W("addon", "R6", "R6_3", "r6_3_fits.tsv")):
    if r["bph_ok"] != "1":
        continue
    OM.setdefault((r["clade"], r["hog"]), {})[r["label"]] = (
        np.log10(fnum(r["wpur_T"]) + .001) - np.log10(fnum(r["wpur_B"]) + .001), fnum(r["ppos_T"]) - fnum(r["ppos_B"]))
PR = {}
for c in PHYLO:
    hs = sorted(h for (cc, h), d in OM.items() if cc == c and "real" in d and "null1" in d)
    assert len(hs) == int(S3[("R6-3a_paired", "Dpur", c)]["n"]) == 400, (c, len(hs))
    PR[c] = (np.array([OM[(c, h)]["real"] for h in hs]), np.array([OM[(c, h)]["null1"] for h in hs]))
    for j, m in enumerate(("Dpur", "Dpos")):
        ex = float((PR[c][0][:, j] - PR[c][1][:, j]).mean())
        assert abs(ex - g3(m, c)[0]) < 1e-9, (c, m, ex)      # paired excess = r6_3_stats.tsv estimate

# ---- blind-sample lineages with their orthogroup counts (same probe file, md5 asserted above) -----------------------
_i = _L.index(next(l for l in _L if l.startswith("per-lineage effects")))
LBD = {}
for l in _L[_i + 2:]:
    m = re.match(r"\s*\d+\s+(\S+)\s+(\S+)\s+(-?[\d.]+)\s+(\d+)\s+(\d+)\s*$", l)
    if not m:
        break
    LBD[m.group(2)] = (m.group(1), float(m.group(3)), int(m.group(4)), int(m.group(5)))
assert len(LBD) == 81 and [(c, k, e) for k, (c, e, _, _) in LBD.items()] == LINB
LAD = {r["canon"]: (r["clade"], fnum(r["eff"]), int(r["n_hog"]), int(r["size"])) for r in LIN}
assert set(LAD) <= set(LBD) and len(LAD) == 65

Wf, Hf = 180.0, 168.0
fig = plt.figure(figsize=(Wf * MM, Hf * MM))
ax = mm_axes(fig, Wf, Hf)
NULLF = "#d9d9d9"                                    # stil5a series colours: deep blue / light grey


def p_lab(p):
    """stil5a group label: stars when the bootstrap P is below 0.01, otherwise the exact P."""
    return ("***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else
            "$P$ = %s" % (("%.3f" % p) if p < 0.1 else ("%.2f" % p)))   # v6: conventional stars


def stil5a(x0, y0, w, h, j, lo, hi, ticks, ylab, letter_xy, key=None):
    a = axes_mm(fig, Wf, Hf, x0, y0, w, h)
    a.set_xlim(-0.55, len(PHYLO) - 0.45); a.set_ylim(lo, hi); a.set_facecolor("none")
    for s_ in ("top", "right"):
        a.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        a.spines[s_].set_linewidth(0.6)
    a.set_yticks(ticks); a.set_yticklabels([(MINUS if t < 0 else "") + ("%g" % abs(t)) for t in ticks])
    a.tick_params(axis="y", length=2, width=0.5, pad=1.0, labelsize=5.0, colors=INK)
    a.set_ylabel(ylab, fontsize=5.5, color=INK, labelpad=1.5)
    a.axhline(0.0, color=NULL, lw=0.4, ls=(0, (3, 2)), zorder=0.5)
    BW, DXB = 0.34, 0.19
    prev_text_top = False
    for g, c in enumerate(PHYLO):
        for k, (vals, col) in enumerate(((PR[c][0][:, j], DEEP2), (PR[c][1][:, j], NULLF))):
            x = g + (-DXB if k == 0 else DXB)
            q1, md, q3 = np.quantile(vals, [.25, .5, .75])
            lw_ = vals[vals >= q1 - 1.5 * (q3 - q1)].min(); hw_ = vals[vals <= q3 + 1.5 * (q3 - q1)].max()
            a.add_patch(Rectangle((x - BW / 2, q1), BW, q3 - q1, fc=col, ec="k", lw=0.5, zorder=3))
            a.plot([x - BW / 2, x + BW / 2], [md, md], color="k", lw=0.8, solid_capstyle="butt", zorder=4)
            a.plot([x, x], [q3, hw_], color="k", lw=0.5, solid_capstyle="butt", zorder=2)
            a.plot([x, x], [q1, lw_], color="k", lw=0.5, solid_capstyle="butt", zorder=2)
            out = vals[(vals < lw_) | (vals > hw_)]
            a.scatter(np.full(len(out), x), out, s=1.6, c="k", lw=0, zorder=4)
            src.append(("ab"[j], c, ("deep", "matched_null")[k], f6(md), f6(q1), f6(q3), len(vals),
                        "box: median, quartiles; whiskers %.4f to %.4f; %d outliers" % (lw_, hw_, len(out))))
        m = ("Dpur", "Dpos")[j]
        pb = fnum(S3[("R6-3a_paired", m, c)]["p_boot"])
        low = pb >= 0.01 and prev_text_top                                         # stil5a: two label rows
        prev_text_top = pb >= 0.01 and not low
        ytop = lo + (0.825 if low else 0.905) * (hi - lo)
        a.plot([g - DXB, g + DXB], [ytop, ytop], color="k", lw=0.5, solid_capstyle="butt", zorder=5)
        a.text(g, ytop + 0.012 * (hi - lo), p_lab(pb), ha="center", va="bottom",
               fontsize=6.0 if pb < 0.01 else 5.0, color=INK, zorder=10)
        e_ = g3(m, c)
        src.append(("ab"[j], c, "paired_excess", f6(e_[0]), f6(e_[1]), f6(e_[2]), 400,
                    "r6_3_stats.tsv; p_boot %s" % S3[("R6-3a_paired", m, c)]["p_boot"]))
    a.set_xticks(range(len(PHYLO)))
    a.set_xticklabels([clab(c) for c in PHYLO], rotation=40, ha="right", rotation_mode="anchor")
    a.tick_params(axis="x", length=2, width=0.5, pad=1.0, labelsize=5.0, colors=INK)
    letter(fig, Wf, Hf, "ab"[j], *letter_xy)
    return a


axa = stil5a(12.0, 8.0, 56.0, 30.0, 0, -2.9, 4.5, [-2, -1, 0, 1, 2],
             "Δpur (log10 $\\omega$, foreground $-$ background)", (1.0, 4.0))
axb = stil5a(12.0, 55.5, 56.0, 21.0, 1, -0.5, 1.0, [-0.4, -0.2, 0, 0.2, 0.4, 0.6],
             "Δpos (sites with $\\omega$ > 1)", (1.0, 52.5))

# ---- stil5a legend: framed box with the two series, the star key and the box key (shared by a and b) -------------
LX0, LY0, LW, LH = 69.6, 13.0, 18.4, 50.0
ax.add_patch(Polygon([(LX0, LY0), (LX0 + LW, LY0), (LX0 + LW, LY0 + LH), (LX0, LY0 + LH)], closed=True, fc="white",
                     ec="k", lw=0.5, zorder=7))
for k, (col, l1, l2) in enumerate(((DEEP2, "Deep", "branches"), (NULLF, "Matched", "null"))):
    yy = LY0 + 3.2 + k * 6.4
    ax.add_patch(Rectangle((LX0 + 1.2, yy - 1.0), 2.6, 2.0, fc=col, ec="k", lw=0.5, zorder=8))
    ax.plot([LX0 + 1.2, LX0 + 3.8], [yy, yy], color="k", lw=0.7, zorder=9)
    ax.plot([LX0 + 2.5] * 2, [yy - 1.0, yy - 1.7], color="k", lw=0.5, zorder=9)
    ax.plot([LX0 + 2.5] * 2, [yy + 1.0, yy + 1.7], color="k", lw=0.5, zorder=9)
    ax.text(LX0 + 4.8, yy - 0.9, l1, ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
    ax.text(LX0 + 4.8, yy + 1.1, l2, ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(LX0 + 1.2, LY0 + 15.6, "n = 400 per", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(LX0 + 1.2, LY0 + 17.8, "clade and arm", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(LX0 + 1.2, LY0 + 21.4, "***", ha="left", va="center", fontsize=6.0, color=INK, zorder=10)
ax.text(LX0 + 4.6, LY0 + 21.4, "$P$ < 0.001", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(LX0 + 1.2, LY0 + 23.6, "(paired excess)", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
# box key (stil5a): outlier dot, upper whisker, Q3, median, Q1, lower whisker, each with an arrow label
KX, KY = LX0 + 2.2, LY0 + 27.0
ax.scatter([KX], [KY - 0.6], s=3.0, c="k", lw=0, zorder=9)
ax.plot([KX] * 2, [KY + 2.2, KY + 5.0], color="k", lw=0.5, zorder=9)
ax.add_patch(Rectangle((KX - 1.3, KY + 5.0), 2.6, 7.0, fc="white", ec="k", lw=0.5, zorder=8))
ax.plot([KX - 1.3, KX + 1.3], [KY + 8.5, KY + 8.5], color="k", lw=0.8, zorder=9)
ax.plot([KX] * 2, [KY + 12.0, KY + 14.8], color="k", lw=0.5, zorder=9)
for yy, lab in ((KY - 0.6, "Outlier"), (KY + 2.2, "Q3 + 1.5 IQR"), (KY + 5.0, "Q3"), (KY + 8.5, "Median"),
                (KY + 12.0, "Q1"), (KY + 14.8, "Q1 $-$ 1.5 IQR")):
    ax.annotate("", xy=(KX + 1.5, yy), xytext=(KX + 3.0, yy), annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>,head_length=0.25,head_width=0.12", lw=0.4, color="k", shrinkA=0,
                                shrinkB=0, mutation_scale=4), zorder=10)
    ax.text(KX + 3.25, yy, lab, ha="left", va="center", fontsize=5.0, color=INK, zorder=10)

# ---- c: spid2b framed horizontal violins of the lineage effects -----------------------------------------------------
CX0, CY0, CW, CH = 28.0, 89.0, 58.0, 21.0
axc = axes_mm(fig, Wf, Hf, CX0, CY0, CW, CH)
XL, XR = -0.65, 1.12
axc.set_xlim(XL, XR); axc.set_ylim(-0.62, 2.62); axc.set_facecolor("none")
for s_ in axc.spines.values():
    s_.set_visible(True); s_.set_linewidth(0.5); s_.set_color("k")
rng_c = np.random.default_rng(20261107)
CROWS = (("Deep lineages,", "first sample (%d)" % len(REAL), REAL, DEEP2, 2),
         ("Deep lineages,", "second sample (%d)" % len(REALB), REALB, DEEP3, 1),
         ("Null pseudo-", "lineages (%d)" % len(NULLE), NULLE, "#bdbdbd", 0))
GRC = np.linspace(XL, XR, 500)
for l1, l2, v, col, yy in CROWS:
    dd = gaussian_kde(v)(GRC)
    hw = dd / dd.max() * 0.40
    axc.fill_between(GRC, yy - hw, yy + hw, facecolor=col, edgecolor="k", lw=0.4, zorder=2)
    vin = v[(v >= XL) & (v <= XR)]
    axc.scatter(vin, yy + rng_c.uniform(-0.09, 0.09, len(vin)), s=1.6, c="k", lw=0, zorder=4)
    q1, md, q3 = np.percentile(v, [25, 50, 75])
    axc.plot([q1, q3], [yy, yy], color="k", lw=0.7, solid_capstyle="butt", zorder=5)
    axc.scatter([md], [yy], s=7.0, c="k", lw=0, zorder=6)
    for x_ in v[v < XL]:          # off-axis lineage: open triangle at the frame edge, value beside it
        axc.scatter([XL + 0.03], [yy], s=6.0, marker="<", facecolor="none", edgecolor="k", linewidths=0.5, zorder=6)
        axc.text(XL + 0.07, yy + 0.30, sv(x_, 2), ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
    axc.text(XR - 0.03, yy + 0.30, "%.0f%% < 0" % (100 * (v < 0).mean()), ha="right", va="center", fontsize=5.0,
             color=INK, zorder=10)
    src.append(("c", "lineage_effects", "%s %s" % (l1, l2), f6(md), f6(q1), f6(q3), len(v),
                "median, quartiles; share < 0 %.3f" % (v < 0).mean()))
axc.set_yticks([2, 1, 0]); axc.set_yticklabels(["%s\n%s" % (r_[0], r_[1]) for r_ in CROWS], linespacing=1.0)
axc.tick_params(axis="y", length=2, width=0.5, pad=1.2, labelsize=5.0, colors=INK)
axc.set_xticks([-0.5, 0.0, 0.5, 1.0]); axc.set_xticklabels([MINUS + "0.5", "0", "0.5", "1.0"])
axc.tick_params(axis="x", length=2, width=0.5, pad=1.2, labelsize=5.0, colors=INK)
axc.set_xlabel("Lineage effect, $\\Delta$ median log10 $K$", fontsize=5.5, color=INK, labelpad=1.5)
axc.axvline(0.0, color=NULL, lw=0.4, ls=(0, (3, 2)), zorder=1)
letter(fig, Wf, Hf, "c", 1.0, 87.0)
for r in LIN:                                     # every plotted lineage value (as v4 c), plus the blind sample
    src.append(("c", "%s|%s" % (r["clade"], r["canon"]), "deep_lineage_first", f6(fnum(r["eff"])), "", "", r["n_hog"],
                "r6_2_lineages.tsv"))
for k, (cl, eb, nb, _) in LBD.items():
    src.append(("c", "%s|%s" % (cl, k), "deep_lineage_second", f6(eb), "", "", nb, "probe_r6_2b.txt"))
for (c, s_, cn), (e, n) in sorted(EFF.items()):
    if s_ == "null1":
        src.append(("c", "%s|%s" % (c, cn), "null_pseudo_lineage", f6(e), "", "", n, "r6_2_fits.tsv, descriptive"))

# ---- d: berg 1c data column, pooled lineage effect under the robustness settings ----------------------------------
DS = [("pooled", "Pooled"), ("lineage_weighted", "Lineage weighted"), ("without_top2_clades", "−Top 2 clades"),
      ("without_largest_lineage", "−Largest lineage")] + [("loco_" + c, "−" + clab(c)) for c in PHYLO
                                                          if "loco_" + c in SM]
DSB = [("pooled", "Second: pooled"), ("without_top2_clades", "Second: −Top 2 clades"),
       ("without_largest_lineage", "Second: −Largest lineage")]
PD, YD0 = 2.5, 125.5
NDA = len(DS)
yd = lambda r: YD0 + r * PD + (1.0 if r >= NDA else 0.0)
NDT = NDA + len(DSB)
YAX_D = yd(NDT - 1) + 0.6 * PD
XS_D, XVD = 31.0, 86.0
UD = berg_col(ax, XS_D, -0.185, 0.025, 205.0, (YD0 - 0.5 * PD, YAX_D), YAX_D, YD0 - 3.4, YD0 - 6.6,
              head1="No effect", head2=("Lower", "Higher"), arrow_to=0.012,
              xlabel="Pooled lineage effect, $\\Delta$ median log10 $K$", ticks=[-0.15, -0.1, -0.05, 0.0])
VD = np.array([SM[k] for k, _ in DS] + [SB[k] for k, _ in DSB])
assert VD[:, 1].min() >= -0.185 and VD[:, 2].max() <= 0.025, VD[:, [1, 2]]
COLS_D = [DEEP if not k.startswith("loco_") else ccol(k[5:]) for k, _ in DS] + [DEEP2] * len(DSB)
for r, ((k, lab), col) in enumerate(zip(DS + DSB, COLS_D)):
    name_row(ax, 3.0, yd(r), lab, XS_D, color=col, fs=5.5)
berg_bars_c(ax, UD, [yd(r) for r in range(NDT)], VD[:, 0], VD[:, 1], VD[:, 2], COLS_D, pitch=PD)
for r, row in enumerate(VD):
    ax.text(XVD, yd(r), sv(row[0], 3), ha="right", va="center", fontsize=5.0, color=INK, zorder=10)
letter(fig, Wf, Hf, "d", 1.0, YD0 - 8.6)
for (k, _), row in zip(DS, VD):
    src.append(("d", k, "R6-2", f6(row[0]), f6(row[1]), f6(row[2]), "", "r6_2_summary.tsv"))
for (k, _), row in zip(DSB, VD[NDA:]):
    src.append(("d", k, "second_sample", f6(row[0]), f6(row[1]), f6(row[2]), "", "probe_r6_2b.txt"))
facts.append(("R6_2_claim", "pooled %.4f (%.4f to %.4f); share<0 %.3f (%.3f to %.3f); without top2 %.4f (%.4f to %.4f)"
              "; SUPPORTED" % (SM["pooled"] + SM["share_lt0"] + SM["without_top2_clades"]), "r6_2_report.txt"))
facts.append(("R6_3a_Dpur", "%.3f (%.3f to %.3f), p_boot %s, %s of 8 clades > 0" % (
    g3("Dpur", "ALL") + (S3[("R6-3a_paired", "Dpur", "ALL")]["p_boot"],
                         S3[("R6-3a_paired", "Dpur", "ALL")]["n_clades_pos"])), "r6_3_stats.tsv"))
facts.append(("R6_3a_Dpos", "%.3f (%.3f to %.3f), p_boot %s, %s of 8 clades > 0" % (
    g3("Dpos", "ALL") + (S3[("R6-3a_paired", "Dpos", "ALL")]["p_boot"],
                         S3[("R6-3a_paired", "Dpos", "ALL")]["n_clades_pos"])), "r6_3_stats.tsv"))

# ---- e: cell3 lineage panel ------------------------------------------------------------------------------------------
TREE = parse_nex_tree(W("final_v3", "date", "lsd.timetree.nex"))
for n_ in postorder(TREE):
    if not n_.children:
        n_.age = 0.0 if n_.age is None else n_.age
    assert n_.age is not None
LABS = {r["code"]: r for r in rtsv(os.path.join(DATA, "meta", "labels_v2.tsv"))}
CHEM = {r["code"]: r for r in rtsv(os.path.join(DATA, "meta", "chemo_habitat.tsv"))}
CAN = {r["canon"]: r for r in rtsv(W("addon", "R6", "R6_2", "canon.tsv"))}
INFO, BADDEP, SPECD = {}, [], []
DINT = {r["species"]: r for r in rtsv(os.path.join(DATA, "meta", "depth_integrated.tsv"))}
for k, (cl, eb, nb, sz) in LBD.items():
    sp = CAN[k]["species"].split(",")
    assert CAN[k]["clade"] == cl and int(CAN[k]["n_species"]) == len(sp) == sz, k
    assert all(LABS[s_]["lab"] == "D" for s_ in sp), k            # every member is labelled deep
    nm = sorted((LABS[s_]["species"], s_) for s_ in sp)
    dep = []
    for s_ in sp:                                 # species median occurrence depth, else specimen depth
        d_ = fnum(LABS[s_]["depth_median"])
        if np.isfinite(d_) and d_ >= 200:
            dep.append(d_)
            continue
        smx = fnum(DINT.get(LABS[s_]["species"], {}).get("sample_max_m"))
        if np.isfinite(smx) and smx >= 200:
            dep.append(smx); SPECD.append(LABS[s_]["species"])
        else:                                     # only a genus-level fallback (or nothing) is on record
            assert not np.isfinite(d_) or abs(d_ - fnum(DINT.get(LABS[s_]["species"], {}).get("genus_median_m"))) < 1e-9                 or DINT.get(LABS[s_]["species"], {}).get("occ_src") == "range", (s_, d_)
            BADDEP.append(LABS[s_]["species"])
    INFO[k] = dict(rep=nm[0][1], name=nm[0][0], n=len(sp), clade=cl, sp=sp,
                   depth=float(10 ** np.mean(np.log10(dep))) if dep else float("nan"),
                   chem=any(CHEM.get(s_, {}).get("chemo") == "1" for s_ in sp))
REP = {v["rep"]: k for k, v in INFO.items()}
for k, v in INFO.items():                         # a lineage is a clade: its MRCA holds no other representative
    tips = mrca(TREE, v["sp"]).tips
    assert not (set(REP) - {v["rep"]}) & tips, k


class LN:
    def __init__(self, age, kids=(), lin=None):
        self.age, self.kids, self.lin, self.y = age, list(kids), lin, None


def prune_l(n):
    if not n.children:
        return LN(0.0, (), REP[n.name]) if n.name in REP else None
    ks = [q for q in (prune_l(c_) for c_ in n.children) if q is not None]
    return None if not ks else ks[0] if len(ks) == 1 else LN(n.age, ks)


def lv(p):
    return [p.lin] if not p.kids else [x for q in p.kids for x in lv(q)]


def rot_l(p):
    for q in p.kids:
        rot_l(q)
    p.kids.sort(key=lambda q: (min(PHYLO.index(INFO[x]["clade"]) for x in lv(q)),
                               max(PHYLO.index(INFO[x]["clade"]) for x in lv(q))))


LT = prune_l(TREE)
rot_l(LT)
ORD = lv(LT)
assert len(ORD) == 81 and sorted(ORD) == sorted(LBD)
CL_ORD = [INFO[k]["clade"] for k in ORD]
assert CL_ORD == sorted(CL_ORD, key=PHYLO.index), "clades not contiguous"
YE0, PE = 13.6, 1.80
ROWY = {k: YE0 + i * PE for i, k in enumerate(ORD)}
XST, XT0, XT1, XSQ, XNM = 90.6, 103.2, 113.0, 113.5, 115.7
XC = (139.8, 144.1, 148.4)
XB = ((151.0, 163.8), (166.4, 178.6))
BL, BR = -0.5, 0.5
DEPTH_BINS = ((1000.0, "<1,000 m", "#9ecae1"), (3000.0, "1,000–3,000 m", "#4a90c4"),
              (5500.0, "3,000–5,500 m", "#1f5f99"), (1e9, ">5,500 m", "#0b2e59"))
dcol = lambda d_: "#ffffff" if not np.isfinite(d_) else next(c_ for t_, _, c_ in DEPTH_BINS if d_ < t_)
YEND = YE0 + 80 * PE
# grey bands for alternate clades (cell3 grey bands), clade labels with silhouettes in the strip
blocks = []
for c in PHYLO:
    ks = [k for k in ORD if INFO[k]["clade"] == c]
    blocks.append((c, ROWY[ks[0]] - PE / 2, ROWY[ks[-1]] + PE / 2, ks))
for i, (c, y0_, y1_, ks) in enumerate(blocks):
    if i % 2 == 0:
        ax.add_patch(Polygon([(XST, y0_), (179.6, y0_), (179.6, y1_), (XST, y1_)], closed=True, fc="#eeeeee",
                             ec="none", zorder=0.2))
    ym = 0.5 * (y0_ + y1_)
    ax.text(XST + 0.3, ym - (1.4 if y1_ - y0_ > 8 else 0.9), clab(c), ha="left", va="center", fontsize=5.0,
            color=ccol(c), zorder=10)
    sh = min(2.8, (y1_ - y0_) - 4.2) if y1_ - y0_ > 6.5 else 2.0
    silhouette_c(ax, PNG[clab(c)], XST + 0.4, ym + 1.4 + (0.0 if y1_ - y0_ > 8 else -0.3), ccol(c),
                 h=min(max(sh, 1.8), 9.5 / aspect(PNG[clab(c)])), z=6)
# lineage tree (LSD2 ages, pruned to one species per lineage), red node diamonds at the clade crowns of the sample
AR = LT.age
XL_ = lambda a_: XT0 + (AR - a_) / AR * (XT1 - XT0)


def sety_l(p):
    if not p.kids:
        p.y = ROWY[p.lin]
    else:
        for q in p.kids:
            sety_l(q)
        p.y = 0.5 * (p.kids[0].y + p.kids[-1].y)


def draw_l(p):
    xp = XL_(p.age)
    ax.plot([xp, xp], [p.kids[0].y, p.kids[-1].y], color="k", lw=0.35, solid_capstyle="projecting", zorder=4)
    for q in p.kids:
        ax.plot([xp, XT1 if not q.kids else XL_(q.age)], [q.y, q.y], color="k", lw=0.35, solid_capstyle="butt",
                zorder=4)
        if q.kids:
            draw_l(q)


sety_l(LT)
draw_l(LT)


def lnodes(p):
    return [p] + [x for q in p.kids for x in lnodes(q)]


for c, _, _, ks in blocks:
    rank = "class" if clab(c) in CLASS_FOCAL else "phylum"
    ks = [k for k in ks if LABS[INFO[k]["rep"]][rank] == clab(c)]
    cand = [p for p in lnodes(LT) if set(ks) <= set(lv(p))]
    nd = min(cand, key=lambda p: len(lv(p)))
    ax.scatter([XL_(nd.age)], [nd.y], s=5.0, marker="D", c=ED8_RED, lw=0, zorder=6)
    facts.append(("e_crown_%s" % c, "%d lineages, sampled crown %.1f Ma" % (len(ks), nd.age), "lsd.timetree.nex"))
# depth squares, names, circles, bars
fmax = {0: 9.0, 1: 200.0, 2: 200.0}
ABBR = []
for k in ORD:
    v = INFO[k]; y = ROWY[k]; col = dcol(v["depth"])
    ax.add_patch(Rectangle((XSQ, y - 0.65), 1.5, 1.3, fc=col, ec="#9a9a9a" if col == "#ffffff" else "none",
                           lw=0.3, zorder=5))
    nm_ = v["name"]
    t_ = ax.text(XNM, y, nm_, ha="left", va="center", fontsize=5.0, style="italic",
                 fontweight="bold" if v["chem"] else "normal", color=SHAL_TXT if v["chem"] else INK, zorder=10)
    if XNM + text_mm(ax, t_)[0] > XC[0] - 2.1 and " " in nm_:      # too long: genus to its initial
        nm_ = nm_[0] + ". " + nm_.split(" ", 1)[1]
        t_.set_text(nm_)
        ABBR.append(v["name"])
    na = LAD.get(k)
    for j, (xc, val) in enumerate(zip(XC, (v["n"], na[2] if na else None, LBD[k][2]))):
        if val is None:
            continue
        rr = 0.25 + 0.62 * np.sqrt(val / fmax[j])
        fc_ = col if col != "#ffffff" else "#e0e0e0"
        ax.add_patch(Circle((xc, y), rr, fc=fc_, ec="none", alpha=0.45, zorder=5))
        ax.text(xc, y, "%d" % val, ha="center", va="center", fontsize=5.0, color=INK, zorder=10)._qa_skip = True
    for j, ((x0_, x1_), eff) in enumerate(zip(XB, (na[1] if na else None, LBD[k][1]))):
        if eff is None:
            continue
        ub = lambda e_: x0_ + (min(max(e_, BL), BR) - BL) / (BR - BL) * (x1_ - x0_)
        ax.add_patch(Rectangle((min(ub(0), ub(eff)), y - 0.36 * PE), abs(ub(eff) - ub(0)), 0.72 * PE,
                               fc=col if col != "#ffffff" else "#bdbdbd", ec="none", zorder=5))
        if eff < BL or eff > BR:                  # off-scale: white break near the edge, value on the other side
            xe = ub(eff)
            xb_ = xe + (1.4 if eff < BL else -1.4)
            ax.plot([xb_ - 0.25, xb_ + 0.25], [y - 0.5, y + 0.5], color="white", lw=0.9, zorder=6)
            ax.text(ub(0) + (0.6 if eff < BL else -0.6), y, sv(eff, 2), ha="left" if eff < BL else "right",
                    va="center", fontsize=5.0, color=INK, zorder=10)
    src.append(("e", INFO[k]["clade"], k, "%s" % (f6(na[1]) if na else ""), f6(LBD[k][1]), "",
                "%d" % v["n"], "rep %s; species %s; depth gm %.0f m; chemo %d; OGs first %s second %d" % (
                    v["name"], ",".join(v["sp"]), v["depth"] if np.isfinite(v["depth"]) else -1, v["chem"],
                    na[2] if na else "NA", LBD[k][2])))
# zero lines, clade means and axes of the bar columns
for j, (x0_, x1_) in enumerate(XB):
    ub = lambda e_: x0_ + (min(max(e_, BL), BR) - BL) / (BR - BL) * (x1_ - x0_)
    ax.plot([ub(0)] * 2, [YE0 - 0.6 * PE, YEND + 0.6 * PE], color="#8c8c8c", lw=0.4, zorder=4.5)
    for c, y0_, y1_, ks in blocks:
        vals = [(LAD[k][1] if k in LAD else None) if j == 0 else LBD[k][1] for k in ks]
        vals = [x for x in vals if x is not None]
        if len(vals) >= 2:
            mval = float(np.mean(vals))
            ax.plot([ub(mval)] * 2, [y0_ + 0.15, y1_ - 0.15], color="k", lw=0.6, zorder=7)
            lab_kn = "%d/%d" % (sum(1 for x in vals if x < 0), len(vals))
            t_kn = ax.text(0, 0, lab_kn, fontsize=5.0); wkn = text_mm(ax, t_kn)[0] + 0.3; t_kn.remove()
            effs = {k: ((LAD[k][1] if k in LAD else None) if j == 0 else LBD[k][1]) for k in ks}
            spots = [(kk, sd) for sd in ("right", "left") for kk in reversed(ks)]
            xm_ = ub(mval)
            for kk, side in spots:                # bottom-most row of the block, right side first, with room for k/n
                e_ = effs[kk]
                if e_ is not None and not (BL <= e_ <= BR):   # that row carries the off-scale value label
                    continue
                lo_x, hi_x = (sorted((ub(0), ub(e_))) if e_ is not None else (ub(0), ub(0)))
                if side == "right" and max(hi_x, xm_) < x1_ - 0.15 - wkn:
                    break
                if side == "left" and min(lo_x, xm_) > x0_ + 0.15 + wkn:
                    break
            else:
                raise SystemExit("no room for k/n in %s" % c)
            ax.text(x1_ - 0.15 if side == "right" else x0_ + 0.15, ROWY[kk], lab_kn, ha=side, va="center",
                    fontsize=5.0, color=INK, zorder=10)
            ci_ = (SM if j == 0 else SB).get("clade_" + c)
            assert ci_ is None or abs(ci_[0] - mval) < 6e-5, (c, j, mval, ci_)
            src.append(("e", c, ("first", "second")[j] + "_clade_mean", f6(mval), f6(ci_[1]) if ci_ else "",
                        f6(ci_[2]) if ci_ else "",
                        "%d/%d" % (sum(1 for x in vals if x < 0), len(vals)), "lineages <0 / summarised"))
    ya = YEND + 0.9 * PE
    ax.plot([ub(BL), ub(BR)], [ya, ya], color="k", lw=0.5, solid_capstyle="butt", zorder=4)
    for t in (-0.4, 0.0, 0.4):
        ax.plot([ub(t)] * 2, [ya, ya + 0.8], color="k", lw=0.5, zorder=4)
        ax.text(ub(t), ya + 1.1, sv(t, 1) if t else "0", ha="center", va="top", fontsize=5.0, color=INK, zorder=10)
ax.text(0.5 * (XB[0][0] + XB[1][1]), YEND + 0.9 * PE + 3.6, "Lineage effect, $\\Delta$ median log10 $K$",
        ha="center", va="top", fontsize=5.5, color=INK, zorder=10)
# pill headers and column labels (cell3)
for (x0_, x1_), lab in (((XC[0] - 2.1, XC[2] + 2.1), "Counts"), ((XB[0][0], XB[1][1]), "Lineage effect")):
    ax.add_patch(FancyBboxPatch((x0_, 3.3), x1_ - x0_, 2.9, boxstyle="round,pad=0,rounding_size=1.4",
                                fc="#dbe6f1", ec="none", zorder=6))
    ax.text(0.5 * (x0_ + x1_), 4.75, lab, ha="center", va="center", fontsize=5.5, color=INK, zorder=10)._qa_skip = True
for xc, l1, l2 in zip(XC, ("Spp.", "OGs", "OGs"), ("", "1st", "2nd")):
    ax.text(xc, 8.6, l1, ha="center", va="center", fontsize=5.0, color=INK, zorder=10)
    if l2:
        ax.text(xc, 10.6, l2, ha="center", va="center", fontsize=5.0, color=INK, zorder=10)
for (x0_, x1_), l1, l2 in zip(XB, ("First sample", "Second sample"), ("(%d lineages)" % len(LAD),
                                                                    "(%d lineages)" % len(LBD))):
    ax.text(0.5 * (x0_ + x1_), 8.6, l1, ha="center", va="center", fontsize=5.0, color=INK, zorder=10)
    ax.text(0.5 * (x0_ + x1_), 10.6, l2, ha="center", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(XNM, 10.6, "Deep lineage", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
# depth legend (cell3 lower left) and the chemosynthetic-habitat key
LGY = YEND + 3.0
for i, (_, lab, col) in enumerate(DEPTH_BINS):
    xx, yy = XST + (i % 2) * 15.5, LGY + (i // 2) * 2.6
    ax.add_patch(Rectangle((xx, yy - 0.65), 1.5, 1.3, fc=col, ec="none", zorder=6))
    ax.text(xx + 2.1, yy, lab, ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.add_patch(Rectangle((XST, LGY + 2 * 2.6 - 0.65), 1.5, 1.3, fc="white", ec="#9a9a9a", lw=0.3, zorder=6))
ax.text(XST + 2.1, LGY + 2 * 2.6, "No species-level depth value", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(XST + 31.5, LGY, "Bold", ha="left", va="center", fontsize=5.0, style="italic", fontweight="bold",
        color=SHAL_TXT, zorder=10)
ax.text(XST + 36.6, LGY, "vent, seep or fall", ha="left", va="center", fontsize=5.0, color=INK, zorder=10)
ax.text(XST + 31.5, LGY + 2.6, "taxon, or chemosymbiotic host", ha="left", va="center", fontsize=5.0, color=INK,
        zorder=10)
letter(fig, Wf, Hf, "e", 89.0, 4.0)

res = qa_v2(fig)
clean = qa_report(res, "Fig3")
facts.append(("e_depth_none", "; ".join(BADDEP) or "none", "no species median >= 200 m and no specimen depth"))
facts.append(("e_depth_specimen", "; ".join(SPECD) or "none", "depth_integrated.tsv sample_max_m used"))
facts.append(("e_abbreviated_names", "; ".join(ABBR), "names shortened to fit"))
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig3_qa_v2.txt"))
write_facts("fig3_source_v5.tsv", src)
write_facts("fig3_facts_v5.tsv", facts)
save_v2(fig, "Fig3")
print("Fig3 written; QA clean:", clean)
