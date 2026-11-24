"""Figure 4: tests of shared genes and functions, and what they can exclude.
a, orthogroups called in >= k clades versus a clade-label permutation null; b, competitive GO tests (QQ);
c, RERconverge accelerated loci versus permutation null across label sets and trees; d, orthologue test coverage by
number of clades and recurrence power (2026-10-30; the CAFE5 panel moved to ED Fig. 7c)."""
import math, re
import numpy as np
import matplotlib.pyplot as plt
from figcommon_sty import *

facts = []
fig = plt.figure(figsize=(180 * MM, 110 * MM))
gs = fig.add_gridspec(2, 2, wspace=0.42, hspace=0.6, left=0.135, right=0.935, top=0.94, bottom=0.1, width_ratios=[1, 1.15])

# a: recurrence
rec = rtsv(os.path.join(BUILD, "src_recurrence.tsv"))
TESTS = [("relax", "RELAX"), ("bph_strict", "BUSTED-PH\nstrict"), ("hyphy_any", "Any HyPhy"), ("cs_qual", "CSUBST\nqualified")]
ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", -0.2)
w = 0.26
for j, k in enumerate(("2", "3", "4")):
    for i, (t, tl) in enumerate(TESTS):
        r = [x for x in rec if x["test"] == t and x["k"] == k][0]
        x = i + (j - 1) * w; o, m, hi = fnum(r["obs_OG"]), fnum(r["null_mean"]), fnum(r["null_975"])
        ax.bar(x, m, w * 0.9, color=LIGHT, lw=0); ax.plot([x - w * 0.4, x + w * 0.4], [hi, hi], color=NULL, lw=0.6)
        ax.plot(x, o, "o", ms=2.8, color=DEEP if fnum(r["perm_p"]) < 0.05 else INK)
        if fnum(r["perm_p"]) < 0.05:
            ax.text(x, o + 4, "P=%s" % ("%.4f" % fnum(r["perm_p"]) if fnum(r["perm_p"]) < 0.001 else "%.3f" % fnum(r["perm_p"])),
                    ha="center", fontsize=4.6, rotation=90, va="bottom")
        if i == 0:
            ax.text(x, -9, "k=%s" % k, ha="center", fontsize=4.6, va="top")
        facts.append(("Fig4a", t, "k" + k, "obs %s null %.2f (97.5%% %s) p %s" % (r["obs_OG"], m, r["null_975"], r["perm_p"])))
ax.set_xticks(range(len(TESTS))); ax.set_xticklabels([t[1] for t in TESTS]); ax.tick_params(axis="x", pad=8)
ax.set_ylim(0, 200); ax.set_ylabel("Orthogroups called in ≥k clades")
ax.bar([0], [0], color=LIGHT, label="Permutation null (mean; tick 97.5%)"); ax.plot([], [], "o", color=DEEP, ms=3, label="Observed")
ax.legend(loc="upper left", fontsize=5, handlelength=1, frameon=False)

# b: competitive GO QQ
ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", -0.14)
cols = {"relax": DEEP, "bph_strict": "#5b84ad", "hyphy_any": SHAL, "cs_qual": INK}
mx = 0
for t, tl in TESTS:
    rr = rtsv(W("addon", "D3_func_competitive", "competitive_%s.tsv" % t))
    p = np.sort(np.array([fnum(r["p"]) for r in rr if not math.isnan(fnum(r["p"]))]))
    n = len(p); e = -np.log10((np.arange(1, n + 1) - 0.5) / n); o = -np.log10(np.clip(p, 1e-300, 1))
    nq = sum(fnum(r["q"]) < 0.05 for r in rr)
    ax.plot(e, o, ".", ms=1.6, color=cols[t], label="%s (%d/%d at q<0.05)" % (" ".join(tl.split()), nq, n))
    mx = max(mx, e.max()); facts.append(("Fig4b", t, "terms %d q<0.05 %d min p %.3g" % (n, nq, p[0]), ""))
ax.plot([0, mx], [0, mx], color=NULL, lw=0.5, ls=":")
ax.set_xlabel("Expected −log$_{10}$ P"); ax.set_ylabel("Observed −log$_{10}$ P (GO term)")
ax.legend(loc="upper left", fontsize=5, markerscale=3, handletextpad=0.1)

# c: RER
rer = [r for r in rtsv(os.path.join(BUILD, "src_rer.tsv")) if r["test"] == "binary_deep" and r["statistic"] == "n_pos_p05"]
RL = {("lab", "fixed"): "All labels", ("lab_nochemo", "fixed"): "No chemosynthetic", ("lab_strict", "fixed"): "Strict labels",
      ("lab_1000", "fixed"): "≥1,000 m labels", ("lab", "rer_v3_ml"): "ML tree", ("lab", "rer_v3_unc"): "Unconstrained tree"}
ax = fig.add_subplot(gs[1, 0]); panel(ax, "c", -0.3)
order = [("lab", "fixed"), ("lab", "rer_v3_ml"), ("lab", "rer_v3_unc"), ("lab_strict", "fixed"), ("lab_nochemo", "fixed"), ("lab_1000", "fixed")]
for i, key in enumerate(order):
    r = [x for x in rer if (x["labels"], x["tree"]) == key][0]
    ax.plot([fnum(r["null_q025"]), fnum(r["null_q975"])], [i, i], color=NULL, lw=2.4, alpha=0.5, solid_capstyle="butt")
    ax.plot(fnum(r["null_mean"]), i, "|", color=NULL, ms=5)
    p = fnum(r["p_upper"]); ax.plot(fnum(r["observed"]), i, "o", ms=3, color=DEEP if p < 0.05 else INK)
    ax.text(124, i, "P=%.3f" % p, va="center", fontsize=5)
    facts.append(("Fig4c", RL[key], "obs %s null %.1f [%s, %s] p %.4f" % (r["observed"], fnum(r["null_mean"]), r["null_q025"], r["null_q975"], p), ""))
ax.set_yticks(range(len(order))); ax.set_yticklabels([RL[k] for k in order]); ax.set_ylim(len(order) - 0.5, -0.5)
ax.set_xlim(15, 140); ax.set_xlabel("Loci with accelerated deep rates (P<0.05)")

# d: what the recurrence test can exclude (2026-10-30 revision; CAFE5 panel moved to ED Fig. 7c).
# Bars: eggNOG orthologues in the RELAX recurrence frame by number of clades in which they were tested
# (work/addon/R4_direction/og_coverage_summary.tsv, script 68). Lines: probability that an orthologue affected in all
# of its k tested clades is significant in >= 3 of them (work/addon/R2/main/recurrence_power.tsv; S10c).
ax = fig.add_subplot(gs[1, 1]); panel(ax, "d", -0.14)
cov = {int(r["n_clades_tested"]): int(r["n_og"]) for r in rtsv(W("addon", "R4_direction", "og_coverage_summary.tsv")) if r["test"] == "relax"}
ks = list(range(1, 9))
kp = [k for k in ks if cov.get(k, 0) > 0]
ax.bar(kp, [cov[k] for k in kp], 0.62, color=NULL, alpha=0.45, lw=0)
for k in ks:
    if cov.get(k, 0) > 0:
        ax.text(k, cov[k] * 1.25, "%d" % cov[k], ha="center", va="bottom", fontsize=4.8, color=INK)
    facts.append(("Fig4d", "coverage", "k%d" % k, "n_og %d" % cov.get(k, 0)))
ax.set_yscale("log"); ax.set_ylim(0.8, 3e5); ax.set_xlim(0.4, 8.6); ax.set_xticks(ks)
ax.set_xlabel("Clades in which an orthologue was tested"); ax.set_ylabel("Orthologues (RELAX frame)")
ax2 = ax.twinx()
RP = rtsv(W("addon", "R2", "main", "recurrence_power.tsv"))
for sc, col, ls, lb in (("K50", DEEP, "-", "K = 0.5, all lineages"), ("K50h", "#5b84ad", "--", "K = 0.5, half of lineages"), ("K75", SHAL, "-", "K = 0.75")):
    rr = sorted([r for r in RP if r["scen"] == sc], key=lambda r: int(r["k_clades"]))
    x = [int(r["k_clades"]) for r in rr]; y = [fnum(r["mean_P_ge3"]) for r in rr]
    ax2.plot(x, y, ls, marker="o", ms=2, lw=0.8, color=col, label=lb)
    for r in rr:
        facts.append(("Fig4d", sc, "k%s" % r["k_clades"], "P_ge3 %s [%s, %s]" % (r["mean_P_ge3"], r["min"], r["max"])))
ax2.axhline(0.8, color=NULL, lw=0.5, ls=":"); ax2.set_ylim(0, 1.05); ax2.set_ylabel("Power to call in ≥3 tested clades")
ax.axvspan(4.5, 8.6, color=LIGHT, alpha=0.25, lw=0)
ax.text(6.55, 1.6, "no orthologue tested\nin more than 4 clades", ha="center", va="bottom", fontsize=4.8, color=INK)
ax2.legend(loc="upper left", fontsize=4.6, frameon=False, handlelength=2, borderaxespad=0.2)

save(fig, "figures", "Fig4")
with open(os.path.join(BUILD, "fig4_facts.tsv"), "w", encoding="utf-8") as f:
    for r in facts:
        f.write("\t".join(r) + "\n")
