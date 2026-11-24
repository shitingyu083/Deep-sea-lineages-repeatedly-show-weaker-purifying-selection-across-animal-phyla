"""Figure 3: selection regimes in deep lineages versus matched shallow nulls.
a, per-clade fraction of loci with a test call in the deep foreground (real) versus the matched shallow null
(mean of null1-5) for RELAX, BUSTED-PH strict and CSUBST perm; b, RELAX calls split into relaxation (K<1) and
intensification (K>1) with the null 95% interval; c, weighted Stouffer z across eight clades with leave-one-clade-out
range; d, odds ratios real versus null from the A4 confound models."""
import math
import numpy as np
import matplotlib.pyplot as plt
from figcommon_sty import *

ORDER = ["Mollusca", "Cnidaria", "Echinodermata", "Annelida", "Arthropoda", "Actinopteri_m", "Chondrichthyes", "Porifera"]
LAB = lambda c: CLADE_LABEL.get(c, c)
T3 = [("RELAX p<0.05", "RELAX"), ("BUSTED-PH strict (FG p<0.05 and diff p<0.05)", "BUSTED-PH strict"),
      ("CSUBST perm p<0.05", "CSUBST perm")]
T5 = [("RELAX p<0.05", "RELAX"), ("BUSTED-PH FG p<0.05", "BUSTED-PH FG"),
      ("BUSTED-PH strict (FG p<0.05 and diff p<0.05)", "BUSTED-PH strict"), ("CSUBST perm p<0.05", "CSUBST perm"),
      ("CSUBST qualified pairs", "CSUBST qualified"), ("CSUBST nochemo qualified pairs", "CSUBST qual., no chemo.")]
facts = []
b25 = {(r["clade"], r["test"]): r for r in rtsv(os.path.join(BUILD, "src_25b_perclade.tsv")) if r["set"] == "rescued"}
c25 = {(r["clade"], r["test"]): r for r in rtsv(os.path.join(BUILD, "src_25c_perclade.tsv"))}

fig = plt.figure(figsize=(180 * MM, 128 * MM))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.02], hspace=0.5, left=0.115, right=0.985, top=0.95, bottom=0.08)
ga = gs[0].subgridspec(1, 5, wspace=0.45, width_ratios=[1, 1, 1, 0.12, 1.25])
gb = gs[1].subgridspec(1, 3, wspace=0.75, width_ratios=[1.1, 1, 1.25])


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


# a: real vs null per clade
for k, (tk, tl) in enumerate(T3):
    ax = fig.add_subplot(ga[0, k])
    if k == 0:
        panel(ax, "a", -0.5)
    for i, c in enumerate(ORDER):
        r = b25[(c, tk)]; rv, nv = fnum(r["real_mean"]), fnum(r["null_mean"])
        if tk.startswith("CSUBST"):
            rv, nv = rv * 100, nv * 100
        p = fnum(c25[(c, tk)]["p_holm_within_test"])
        ax.plot([nv, rv], [i, i], color=LIGHT, lw=1.2, zorder=1)
        ax.scatter([nv], [i], s=9, color=NULL, zorder=2, lw=0)
        ax.scatter([rv], [i], s=11, color=DEEP if rv > nv else SHAL, zorder=3, lw=0)
        ax.text(1.0, i, stars(p), transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=6)
        facts.append(("Fig3a", tl, c, "real %.4f null %.4f Holm %.3g" % (rv, nv, p)))
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([LAB(c) for c in ORDER] if k == 0 else [])
    ax.set_ylim(len(ORDER) - 0.5, -0.5); ax.set_title(tl, fontsize=6.5, pad=3)
    ax.set_xlabel("Loci called (%)" if tk.startswith("CSUBST") else "Loci called (fraction)")
    if tk.startswith("CSUBST"):
        ax.set_xscale("symlog", linthresh=1); ax.set_xlim(-0.3, 25); ax.set_xticks([0, 1, 5, 20]); ax.set_xticklabels(["0", "1", "5", "20"]); ax.minorticks_off()
    elif k == 0:
        ax.set_xlim(0.15, 0.7)
ax0 = fig.axes[0]
ax0.scatter([], [], s=11, color=DEEP, label="Deep (real)"); ax0.scatter([], [], s=9, color=NULL, label="Matched shallow null")
ax0.legend(loc="lower right", bbox_to_anchor=(1.06, -0.01), fontsize=5, handletextpad=0.1, borderaxespad=0.1)

# c (top right): combined z with LOCO range
loco = {}
for r in rtsv(os.path.join(BUILD, "src_25c_derived_sensitivity.tsv")):
    if r["variant"] == "LOCO":
        loco.setdefault(r["test"], []).append(fnum(r["z"]))
comb = {r["test"]: r for r in rtsv(W("combine25c", "final_main.combined.tsv"))}
ax = fig.add_subplot(ga[0, 4]); panel(ax, "b", -0.62)
for i, (tk, tl) in enumerate(T5):
    z = fnum(comb[tk]["stouffer_z"]); pr = comb[tk]["set"] == "primary"
    col = DEEP if pr else NULL
    if tk in loco:
        ax.plot([min(loco[tk]), max(loco[tk])], [i, i], color=col, lw=2.2, alpha=0.35, solid_capstyle="butt")
    ax.plot(z, i, "o" if pr else "s", color=col, ms=3.2)
    ax.text(z + 0.25, i - 0.32, "%s/%s" % (comb[tk]["n_clades_p<0.05"], comb[tk]["k_clades"]), fontsize=4.8, va="center")
    facts.append(("Fig3b", tl, "z %.3f LOCO %s Holm %s" % (z, ("%.3f-%.3f" % (min(loco[tk]), max(loco[tk]))) if tk in loco else "NA", comb[tk]["p_holm"]), ""))
ax.axvline(norm_z := 1.6449, color=NULL, lw=0.5, ls=":")
ax.set_yticks(range(len(T5))); ax.set_yticklabels([t[1] for t in T5]); ax.set_ylim(len(T5) - 0.5, -0.6)
ax.set_xlim(0, 8.3); ax.set_xlabel("Combined Stouffer z (8 clades)")
ax.set_title("Bar, leave-one-clade-out range", fontsize=6, pad=3)

# b (bottom left): relaxed vs intensified
rd = {(r["clade"], r["ind"]): r for r in rtsv(os.path.join(BUILD, "src_relax_direction.tsv"))}
ax = fig.add_subplot(gb[0, 0]); panel(ax, "c", -0.45)
rows = ORDER + ["all"]
for i, c in enumerate(rows):
    for j, (ind, col) in enumerate((("relax_relaxed", DEEP), ("relax_intens", SHAL))):
        r = rd.get((c, ind))
        if r is None:
            continue
        y = i + (j - 0.5) * 0.34
        ax.plot([fnum(r["lo"]), fnum(r["hi"])], [y, y], color=NULL, lw=2, alpha=0.5, solid_capstyle="butt")
        ax.plot(fnum(r["real_rate"]), y, "o", color=col, ms=2.8 if c != "all" else 3.6)
        facts.append(("Fig3c", c, ind, "real %.4f null %.4f [%.4f, %.4f]" % tuple(fnum(r[k]) for k in ("real_rate", "rate", "lo", "hi"))))
ax.set_yticks(range(len(rows))); ax.set_yticklabels([LAB(c) if c != "all" else "Pooled" for c in rows])
ax.set_ylim(len(rows) - 0.5, -0.5); ax.axhline(len(rows) - 1.5, color=LIGHT, lw=0.4)
ax.set_xlabel("Loci with RELAX p<0.05 (fraction)")
ax.plot([], [], "o", color=DEEP, ms=3, label="Relaxed (K<1)"); ax.plot([], [], "o", color=SHAL, ms=3, label="Intensified (K>1)")
ax.plot([], [], color=NULL, lw=2, alpha=0.5, label="Null 95% interval")
ax.legend(loc="center right", fontsize=5, handlelength=1.2, borderaxespad=0.1)

# d: A4 odds ratios
mods = {}
for r in rtsv(W("addon", "A4", "models.tsv")):
    if r["term"] == "is_real":
        mods[(r["ind"], r["model"])] = r
IND = [("relax", "RELAX (any)"), ("relax_relaxed", "RELAX relaxed"), ("relax_intens", "RELAX intensified"),
       ("bph", "BUSTED-PH FG"), ("bph_strict", "BUSTED-PH strict"), ("cs_perm", "CSUBST perm"), ("cs_qual", "CSUBST qualified")]
MS = [("M0", "Clade only", LIGHT, "o"), ("M1", "Primary (+ annotation route, tips, gap)", DEEP, "o"),
      ("M2", "+ background route, TSA, length", "#5b84ad", "s"), ("C1", "Within-HOG conditional", INK, "D")]
ax = fig.add_subplot(gb[0, 1:]); panel(ax, "d", -0.3)
for i, (ind, il) in enumerate(IND):
    for j, (m, ml, col, mk) in enumerate(MS):
        r = mods[(ind, m)]; e = fnum(r["est"]); se = fnum(r["se_cl"]) if m != "C1" else fnum(r["se"])
        y = i + (j - 1.5) * 0.18
        ax.plot([math.exp(e - 1.96 * se), math.exp(e + 1.96 * se)], [y, y], color=col if m != "M0" else NULL, lw=0.6)
        ax.plot(math.exp(e), y, mk, color=col, mec=col if m != "M0" else NULL, mew=0.4, ms=2.6, label=ml if i == 0 else None)
        facts.append(("Fig3d", ind, m, "OR %.3f [%.3f, %.3f] p %s" % (math.exp(e), math.exp(e - 1.96 * se), math.exp(e + 1.96 * se), r["p_cl"] if m != "C1" else r["p"])))
ax.axvline(1, color=NULL, lw=0.5, ls=":"); ax.set_xscale("log")
ax.set_xticks([0.25, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8"]); ax.minorticks_off()
ax.set_yticks(range(len(IND))); ax.set_yticklabels([t[1] for t in IND]); ax.set_ylim(len(IND) - 0.5, -0.5)
ax.set_xlabel("Odds ratio, deep foreground vs matched null (95% CI, HOG-clustered)")
ax.legend(loc="lower right", fontsize=5, handlelength=0.8, borderaxespad=0.1, handletextpad=0.3)
for i in range(1, len(IND)):
    ax.axhline(i - 0.5, color=LIGHT, lw=0.3)

save(fig, "figures", "Fig3")
with open(os.path.join(BUILD, "fig3_facts.tsv"), "w", encoding="utf-8") as f:
    for r in facts:
        f.write("\t".join(r) + "\n")
