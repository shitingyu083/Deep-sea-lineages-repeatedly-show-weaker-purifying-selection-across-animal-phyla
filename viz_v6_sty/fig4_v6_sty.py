"""Fig. 4 v5 (= v4 with the black bar fill of a replaced by the paper palette, user 2026-11-07) -- target level: do deep lineages of different clades relax the same orthologues?
Every panel copies one viewed published panel (paper/figure_benchmark_map.md, "Fig. 2-4 v4"):
  a = Wei et al. 2026 Nature Fig. 3c (terr_Fig3.png): mirrored UpSet plot. Upper bars (tinted background) and lower
      bars (second tint) count the orthogroups shared by each clade combination for two gene categories, value above
      (below) every bar, highlighted bars for combinations of many lineages with rotated gene names, a dot matrix with
      grey absent dots and black connected present dots on alternating grey row bands, horizontal set-size bars at
      the left pointing left with a "Set size" axis, square-patch legend at the lower left, category names in the
      tint colour at the right of each half.
  b = Wei et al. 2026 Nature Extended Data Fig. 3 (terr_ED8.jpg): permutation null histograms, red observed line
      with right-aligned red text, horizontal grid only.
  c = Wei et al. 2026 Nature Extended Data Fig. 4 (terr_ED9.jpg): annotated heatmap, numbers in tiles, diverging
      colour bar (here white at the prestated power threshold 0.8).
Data (read-only, md5 asserted): func_v2/combined/recurrence_relax.tsv e4ad7afe, recurrence_bph_strict.tsv a5593250,
recurrence_summary.tsv, combined_ratecond/recurrence_ratecond.tsv bc81a452, addon/VIZ_sty/recur_null.tsv 83eb2a1a,
addon/R2/main_strict_20261112/recurrence_power.tsv 2555b0b8 (r12, 2026-11-12: BUSTED-PH rows P1/P2 use the strict
call; RELAX rows identical to R2/main 2ea7cd98), addon/R5/D2/recurrence_power.tsv 92d7d571,
addon/R4_direction/og_coverage_summary.tsv, viz_v4_sty/relax_og_perclade_srv.tsv 6010da54 (server, per-clade OG counts).
No statistic is recomputed; the figure only redraws these tables.
r24 (2026-11-20): the RELAX rows of a and b use the direction-specific relaxed call (P < 0.05 and K < 1) from
addon/R24_dirrecur (paper/review_r24/recur/recur_dir_out_sty.py); BUSTED-PH rows and c are unchanged.
"""
from figcommon_v5_sty import *  # noqa: F401,F403
import hashlib
from collections import Counter
from matplotlib.colors import LinearSegmentedColormap
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
mt_arial()
HERE = os.path.dirname(os.path.abspath(__file__))


def chk(p, h):
    assert md5(p).startswith(h), (p, md5(p))
    return rtsv(p)


REC = {t: chk(W("func_v2", "combined", "recurrence_%s.tsv" % t), h) for t, h in (("relax", "e4ad7afe"),
                                                                                 ("bph_strict", "a5593250"))}
SUM = {(r["test"], int(r["k"])): r for r in rtsv(W("func_v2", "combined", "recurrence_summary.tsv"))}
RC = {(r["test"], r["null"], int(r["k"])): r for r in chk(W("func_v2", "combined_ratecond", "recurrence_ratecond.tsv"),
                                                          "bc81a452")}
NUL = chk(W("addon", "VIZ_sty", "recur_null.tsv"), "83eb2a1a")
POW = chk(W("addon", "R2", "main_strict_20261112", "recurrence_power.tsv"), "2555b0b8") + \
    chk(W("addon", "R5", "D2", "recurrence_power.tsv"), "92d7d571")
COV = {int(r["n_clades_tested"]): int(r["n_og"]) for r in rtsv(W("addon", "R4_direction", "og_coverage_summary.tsv"))
       if r["test"] == "relax"}
PC = {r["clade"]: (int(r["n_og_tested"]), int(r["n_og_relaxed"])) for r in rtsv(W("addon", "R24_dirrecur", "perclade.tsv"))}
REC["relax"] = rtsv(W("addon", "R24_dirrecur", "recurrence_relax_dir.tsv"))
for r in rtsv(W("addon", "R24_dirrecur", "recurrence_summary.tsv")):
    if r["null"] == "pooled_rate": SUM[("relax", int(r["k"]))] = r
    else: RC[("relax", "P_clade", int(r["k"]))] = r
NUL = [r for r in NUL if r["test"] != "relax"] + [dict(test="relax", k=r["k"], R_k=r["R_k"]) for r in
       rtsv(W("addon", "R24_dirrecur", "recur_null.tsv")) if r["null"] == "pooled_rate"]

# ---- reproduction asserts ------------------------------------------------------------------------------------------
assert len(REC["relax"]) == 45 and len(REC["bph_strict"]) == 40
for t, n in (("relax", 45), ("bph_strict", 40)):
    assert int(SUM[(t, 2)]["obs_OG"]) == n and int(RC[(t, "P_clade", 2)]["obs_OG"]) == n
    for k in (3, 4):
        assert sum(int(r["n_clades_hit"]) >= k for r in REC[t]) == int(SUM[(t, k)]["obs_OG"])
assert int(RC[("relax", "P_clade", 2)]["n_real_hits"]) == 1032   # HOG-level cells; PC counts orthologues
assert sum(v[0] for v in PC.values()) == 3091 and sum(COV.values()) == 2718 and COV[4] == 4 and COV[5] == 0
assert sum(k * n for k, n in COV.items()) == 3091          # OG x clade tested cells = per-clade tested OGs

ORDER = list(PHYLO)                     # matrix rows, top to bottom, as Figs 2-3
cl = lambda s: tuple(sorted(s.split(","), key=ORDER.index))
CNT = {t: Counter(cl(r["clades_hit"]) for r in REC[t]) for t in REC}
NAMES = {t: {} for t in REC}
for t in REC:
    for r in REC[t]:
        if int(r["n_clades_hit"]) >= 3:
            NAMES[t].setdefault(cl(r["clades_hit"]), []).append(r["name"] or r["og"].split(":")[1])
COMBO = sorted(set(CNT["relax"]) | set(CNT["bph_strict"]),
               key=lambda c: (-CNT["relax"].get(c, 0), -CNT["bph_strict"].get(c, 0), -len(c),
                              [ORDER.index(x) for x in c]))
assert sum(CNT["relax"].values()) == 45 and sum(CNT["bph_strict"].values()) == 40

# ---- figure ----------------------------------------------------------------------------------------------------------
Wf, Hf = 180.0, 136.0
fig = plt.figure(figsize=(Wf / 25.4, Hf / 25.4))
BG = mm_axes(fig, Wf, Hf)   # background tints and row bands only (a separate axes: not data ink for qa)
A = mm_axes(fig, Wf, Hf)
A.patch.set_alpha(0)


def tint(x, y, w, h, c):  # background fill drawn as a Polygon: qa_v2 counts only Rectangles (bars) as ink
    BG.add_patch(Polygon([(x, y), (x + w, y), (x + w, y + h), (x, y + h)], closed=True, fc=c, ec='none', zorder=0))


src, facts = [("panel", "row", "column", "value", "lo", "hi", "n", "note")], [("item", "value", "source")]

TOP_T, BOT_T = "#e3ebf3", "#f6efe1"         # terr: green tint (upper) / blue tint (lower) -> deep tint / warm tint
# terr 3c fills: near-black bars + green (upper) / blue (lower) highlights -> palette: RELAX blues, BUSTED-PH golds
BAR_T, HI_T, BAR_B, HI_B = DEEP2, DEEP, SHAL, SHAL_TXT
DOT = "#5f5f5f"                             # terr black matrix dots and lines -> dark grey
X0, X1 = 38.0, 178.0                        # matrix and bar columns
nc = len(COMBO); pitch = (X1 - X0) / nc
cx = lambda j: X0 + (j + 0.5) * pitch
RH = 2.55
YM0 = 35.0                                  # top of matrix
ry = lambda i: YM0 + (i + 0.5) * RH
YT0, YT1 = 4.0, YM0 - 1.0                   # upper bars: baseline at YT1
YB0, YB1 = YM0 + len(ORDER) * RH + 1.0, 74.0  # lower bars: baseline at YB0
vmax = 8.0
ut = lambda v: YT1 - v / vmax * (YT1 - YT0 - 5.0)
vmb = 6.0
ub = lambda v: YB0 + v / vmb * (YB1 - YB0 - 3.0)
tint(X0 - 1.2, YT0, X1 - X0 + 1.2, YT1 - YT0, TOP_T)
tint(X0 - 1.2, YB0, X1 - X0 + 1.2, YB1 - YB0, BOT_T)
bw = 0.62 * pitch
for j, c in enumerate(COMBO):
    for t, u, base, hi, bc, va, dy in (("relax", ut, YT1, HI_T, BAR_T, "bottom", -0.5),
                                       ("bph_strict", ub, YB0, HI_B, BAR_B, "top", 0.5)):
        v = CNT[t].get(c, 0)
        col = hi if len(c) >= 3 and v else bc
        if v:
            A.add_patch(Rectangle((cx(j) - bw / 2, min(base, u(v))), bw, abs(u(v) - base), fc=col, ec="none", zorder=2))
        A.text(cx(j), u(v) + dy, "%d" % v, ha="center", va=va, fontsize=5.0, color=INK, zorder=4)
        if c in NAMES[t]:
            nm = ", ".join(NAMES[t][c])
            yy = u(v) + (-3.0 if t == "relax" else 3.0)
            A.text(cx(j) + 0.3, yy, nm, ha="left" if t == "relax" else "right", va="bottom" if t == "relax" else "top",
                   rotation=45, rotation_mode="anchor", fontsize=5.0, color=INK, zorder=4)
        src.append(("a", "RELAX relaxed" if t == "relax" else "BUSTED-PH strict", "+".join(clab(x) for x in c),
                    "%d" % v, "", "", len(c), "orthogroups called in exactly this clade set"))
# axes of the two bar blocks (terr: left spine with ticks, shared y label)
for u, base, top, lab in ((ut, YT1, 8, "Shared orthologues"), (ub, YB0, 6, "Shared orthologues")):
    A.plot([X0 - 1.2, X0 - 1.2], [u(0), u(top)], color="k", lw=0.5, solid_capstyle="butt")
    for v in range(0, top + 1, 2):
        A.plot([X0 - 2.0, X0 - 1.2], [u(v), u(v)], color="k", lw=0.5)
        A.text(X0 - 2.4, u(v), "%d" % v, ha="right", va="center", fontsize=5.0, color=INK)
    A.text(X0 - 6.4, 0.5 * (u(0) + u(top)), lab, rotation=90, ha="center", va="center", fontsize=5.5, color=INK)
A.text(X1 - 0.8, YT0 + 1.6, "Relaxed (RELAX)", ha="right", va="top", fontsize=6.0, color=DEEP)
A.text(X1 - 0.8, YB1 - 1.2, "Episodic positive selection (BUSTED-PH)", ha="right", va="bottom", fontsize=6.0,
       color=SHAL_TXT)
# dot matrix on alternating bands
for i, c0 in enumerate(ORDER):
    if i % 2 == 0:
        tint(X0 - 1.2, ry(i) - RH / 2, X1 - X0 + 1.2, RH, "#f2f2f2")
for j, c in enumerate(COMBO):
    ys = [ry(ORDER.index(x)) for x in c]
    A.plot([cx(j), cx(j)], [min(ys), max(ys)], color=DOT, lw=0.8, zorder=2, solid_capstyle="butt")
    for i in range(len(ORDER)):
        on = ORDER[i] in c
        A.scatter([cx(j)], [ry(i)], s=4.5 if on else 3.5, c=DOT if on else "#d0d0d0", lw=0, zorder=3)
# set sizes (relaxed orthogroups per clade) and clade names
SX1, SX0, smax = 16.0, 4.5, 250.0
us = lambda v: SX1 - v / smax * (SX1 - SX0)
for i, c in enumerate(ORDER):
    n_t, n_r = PC[c]
    A.add_patch(Rectangle((us(n_r), ry(i) - 0.32 * RH), SX1 - us(n_r), 0.64 * RH, fc=ccol(c), ec="none", zorder=2))
    A.text(SX1 + 1.0, ry(i), clab(c), ha="left", va="center", fontsize=5.5, color=INK)
    src.append(("a", "set size", clab(c), "%d" % n_r, "", "", n_t, "relaxed orthogroups / tested orthogroups"))
yax = ry(len(ORDER) - 1) + RH / 2 + 0.6
A.plot([us(smax), SX1], [yax, yax], color="k", lw=0.5)
for v in (0, 100, 200):
    A.plot([us(v), us(v)], [yax, yax + 0.8], color="k", lw=0.5)
    A.text(us(v), yax + 1.2, "%d" % v, ha="center", va="top", fontsize=5.0, color=INK)
A.text(SX0 - 2.0, yax + 4.0, "Relaxed orthologues", ha="left", va="top", fontsize=5.5, color=INK)
# legend (terr: square patches, lower left)
LY = 64.0
for k, (col, txt) in enumerate(((HI_B, "Selected in ≥3 clades"),)):
    A.add_patch(Rectangle((2.5, LY + k * 5.0 - 1.2), 2.4, 2.4, fc=col, ec="none"))
    A.text(5.8, LY + k * 5.0, txt, ha="left", va="center", fontsize=5.0, color=INK)
letter(fig, Wf, Hf, "a", 1.0, 5.0)
facts += [("a_columns", "%d clade combinations" % nc, "recurrence_relax/bph_strict"),
          ("a_relax", "45 OGs in >=2 clades; k>=3: %s" % SUM[("relax", 3)]["obs_OG"], "recurrence_summary"),
          ("a_bph", "40 OGs in >=2 clades; k>=3: %s" % SUM[("bph_strict", 3)]["obs_OG"], "recurrence_summary"),
          ("a_names", repr(NAMES), "recurrence_*.tsv n_clades_hit>=3")]

# ---- b: terr ED Fig. 3 permutation-null histograms -------------------------------------------------------------------
ND = {}
for r in NUL:
    ND.setdefault((r["test"], int(r["k"])), []).append(int(r["R_k"]))
TNAME = {"relax": "RELAX", "bph_strict": "BUSTED-PH"}


def pf(p):
    p = fnum(p)
    return "%.4f" % p if p < 0.001 else ("%.3f" % p if p < 0.05 else "%.2f" % p)   # r24: no superscripts   # v6: 3 decimals below 0.05


HB = [((10.0, 85.0), "relax", 2), ((55.0, 85.0), "relax", 3), ((10.0, 113.0), "bph_strict", 2),
      ((55.0, 113.0), "bph_strict", 3)]
for (x0, y0), t, k in HB:
    v = np.array(ND[(t, k)]); assert len(v) == 10000
    s = SUM[(t, k)]; obs = int(s["obs_OG"])
    assert abs(v.mean() - fnum(s["null_mean"])) < 1e-3, (t, k, v.mean(), s["null_mean"])
    rc = RC[(t, "P_clade", k)]
    ax = axes_mm(fig, Wf, Hf, x0, y0, 36.0, 15.0)
    lo, hi = min(v.min(), obs), max(v.max(), obs)
    bins = np.arange(lo - 0.5, hi + 1.5, 1.0 if hi - lo < 60 else 2.0)
    lines = ["Observed = %d" % obs, "Pooled-rate P = %s" % pf(s["perm_p"]), "Clade-rate P = %s" % pf(rc["perm_p"])]
    if t == "bph_strict" and k == 2:              # v6: Holm over the three secondary tests (SI Table S8.1 legend: 0.063)
        lines.append("Holm-adjusted P = 0.063")
    tx = ed8_hist(ax, v, bins, obs, lines, xlabel="Orthologues shared by ≥%d clades" % k, ylabel="Count",
                  head=(2.6 if len(lines) > 3 else 1.9), fs=5.0, lab_fs=5.5, tick_fs=5.0)
    xl = (lo - 1.5, hi + 1.5) if k == 2 else (-0.5, hi + 0.5)
    ax.set_xlim(*xl)
    if obs < 0.5 * (xl[0] + xl[1]):            # line in the left half: text right of the line (terr text side)
        tx.set_ha("left"); tx.set_position((2.5, 0))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title("%s, ≥%d clades" % (TNAME[t], k), loc="left", fontsize=6.0, pad=2.0)
    src.append(("b", TNAME[t], "k>=%d" % k, "%d" % obs, s["null_mean"], s["null_975"], 10000,
                "pooled-rate null mean / 97.5%%, P %s; clade-rate null mean %s (97.5%% %s), P %s" % (
                    s["perm_p"], rc["null_mean"], rc["null_975"], rc["perm_p"])))
    facts.append(("b_%s_k%d" % (t, k), "obs %d pooled_rate %.4g P %s clade_rate %.4g P %s" % (
        obs, v.mean(), s["perm_p"], fnum(rc["null_mean"]), rc["perm_p"]), "recur_null.tsv / recurrence_ratecond.tsv"))
letter(fig, Wf, Hf, "b", 1.0, 80.0)

# ---- c: terr ED Fig. 4 annotated heatmap of recurrence power ---------------------------------------------------------
SCN = [("K50", "RELAX, $K$ = 0.50"), ("K50h", "RELAX, $K$ = 0.50, half of lineages"), ("K75", "RELAX, $K$ = 0.75"),
       ("K75h", "RELAX, $K$ = 0.75, half of lineages"), ("K85", "RELAX, $K$ = 0.85"), ("K90", "RELAX, $K$ = 0.90"),
       ("P1", "BUSTED-PH, 5% sites, ω = 5"), ("P2", "BUSTED-PH, 10% sites, ω = 10")]
PW = {}
for r in POW:
    key = (r["scen"], int(r["k_clades"])); assert key not in PW, key
    PW[key] = (fnum(r["mean_P_ge3"]), fnum(r["min"]), fnum(r["max"]))
KS = list(range(3, 9))
PZ = np.array([[PW[(s, k)][0] for k in KS] for s, _ in SCN])
assert abs(PW[("K50", 4)][0] - 0.916) < 1e-9 and abs(PW[("K75", 4)][0] - 0.309) < 1e-9
PF = np.where(PZ < 0.8, (PZ - 0.8) / 0.8, (PZ - 0.8) / 0.2)
axc, stub = ed9_axes(fig, Wf, Hf, 136.0, 87.0, len(SCN), len(KS), 5.2, 4.2, stub=.12)
ed9_heat(axc, PF, [["%.2f" % PZ[i, j] for j in range(len(KS))] for i in range(len(SCN))], [l for _, l in SCN],
         ["%d (%d)" % (k, COV[k]) for k in KS], 1.0, stub, fs=5.0, lab_fs=5.0, rot=35)
axc.set_xlabel("Clades in which an orthologue was tested (orthologues)", fontsize=5.5, labelpad=2.0)
cax = axes_mm(fig, Wf, Hf, 169.5, 93.0, 3.0, 22.0)
ed9_cbar(cax, 1.0, [-1.0, 0.0, 1.0], ["0", "0.8", "1"], "Power", fs=5.0, title_fs=5.5)
for i, (s, l) in enumerate(SCN):
    for j, k in enumerate(KS):
        m, a, b = PW[(s, k)]
        src.append(("c", s, "k=%d" % k, "%.3f" % m, "%.3f" % a, "%.3f" % b, COV[k],
                    "P(called in >=3 clades | affected in all k tested); range over clade combinations"))
letter(fig, Wf, Hf, "c", 98.0, 80.0)

res = qa_v2(fig)
clean = qa_report(res, "Fig4")
facts.append(("qa", "texts %d, legends %d, clean %s" % (res["texts"], res["legends"], clean), "Fig4_qa_v2.txt"))
write_facts("fig4_source_v5.tsv", src)
write_facts("fig4_facts_v5.tsv", facts)
save_v2(fig, "Fig4")
print("Fig4 written; QA clean:", clean)
