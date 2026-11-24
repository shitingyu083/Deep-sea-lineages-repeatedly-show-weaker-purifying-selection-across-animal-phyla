"""Figure 2: proteome composition. a, PGLS deep effects of six traits across specifications; b, contrasts at the
65 disjoint deep origins; c, clade-centred f_T against depth; d, gene-level sign fraction versus the calibrated null."""
import math, re
import numpy as np
import matplotlib.pyplot as plt
from figcommon_sty import *

P6 = ["busco_f_T", "whole_pI", "whole_R_RK", "whole_acid_base", "whole_net_charge7_per_res", "busco_f_F"]
GMAP = {"busco_f_T": "f_T", "whole_pI": "pI", "whole_R_RK": "R_RK", "whole_acid_base": "acid_base",
        "whole_net_charge7_per_res": "net_charge7_per_res", "busco_f_F": "f_F"}
SHORT = {"busco_f_T": "$f_T$", "whole_pI": "pI", "whole_R_RK": "R/(R+K)", "whole_acid_base": "Acid/base",
         "whole_net_charge7_per_res": "Net charge", "busco_f_F": "$f_F$"}
facts = []

ctl = {k: {r["trait"]: r for r in rtsv(W("final_v3", d, "controls_summary.tsv"))}
       for k, d in (("lab", "pgls_controls"), ("strict", "pgls_controls_strict"), ("1000", "pgls_controls_1000"))}
reann = {}
for l in open(W("addon", "A1_reannot", "refit_report.txt"), encoding="utf-8"):
    f = l.rstrip("\n").split("\t")
    if len(f) == 8 and f[1] == "S_route3":
        est = float(f[5].split()[0]); lo, hi = [float(x) for x in re.findall(r"[-0-9.e]+", f[6])[:2]]
        reann[f[0]] = (est, lo, hi)
SPECS = [("All labels", "lab", "all", DEEP, "o"), ("+ habitat temperature", "lab", "deepT", "#5b84ad", "s"),
         ("No chemosynthetic taxa", "lab", "nochemo", "#8fb3d4", "D"), ("Strict labels", "strict", "all", SHAL, "^"),
         ("Labels ≥ 1,000 m", "1000", "all", "#a07020", "v"), ("Re-annotated routes", None, None, INK, "P")]

fig = plt.figure(figsize=(180 * MM, 118 * MM))
gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], height_ratios=[1, 1], wspace=0.55, hspace=0.6,
                      left=0.115, right=0.985, top=0.93, bottom=0.17)

# a: forest
ax = fig.add_subplot(gs[:, 0]); panel(ax, "a", -0.36, 1.0)
for i, t in enumerate(P6):
    for j, (lab, ls, sp, col, mk) in enumerate(SPECS):
        y = i + (j - 2.5) * 0.13
        if ls is None:
            e, lo, hi = reann[t]
        else:
            r = ctl[ls][t]; e, se = fnum(r[sp + "_est"]), fnum(r[sp + "_se"]); lo, hi = e - 1.96 * se, e + 1.96 * se
        ax.plot([lo, hi], [y, y], color=col, lw=0.6); ax.plot(e, y, mk, ms=2.6, color=col, mec=col, label=lab if i == 0 else None)
        facts.append(("Fig2a", t, lab, "%.4f [%.4f, %.4f]" % (e, lo, hi)))
ax.axvline(0, color=NULL, lw=0.5, ls=":")
ax.set_yticks(range(len(P6))); ax.set_yticklabels([SHORT[t] for t in P6]); ax.invert_yaxis()
ax.set_xlabel("Deep effect (s.d., 95% CI)")
ax.legend(loc="upper center", bbox_to_anchor=(0.4, -0.07), ncol=2, handlelength=1, columnspacing=0.8, fontsize=5.5)
for i in range(1, len(P6)):
    ax.axhline(i - 0.5, color=LIGHT, lw=0.3)

# b: origin contrasts
con = [r for r in rtsv(W("addon", "A2", "contrasts.tsv")) if r["set"] == "disjoint" and r["adj"] == "raw"]
summ = {(r["trait"], r["set"], r["adj"]): r for r in rtsv(W("addon", "A2", "summary.tsv"))}
ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", -0.3)
rng = np.random.default_rng(7)
for i, t in enumerate(P6):
    v = np.array([fnum(r["contrast"]) for r in con if r["trait"] == t])
    s = summ[(t, "disjoint", "raw")]
    col = DEEP if t == "busco_f_T" else NULL
    ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), np.clip(v, -2.5, 2.5), s=2.5, color=col, alpha=0.6, lw=0)
    ax.plot([i - 0.3, i + 0.3], [fnum(s["median"])] * 2, color=INK, lw=1)
    hp = fnum(s["wilcox_holm"])
    ax.text(i, 2.7 + 0.55 * (i % 2), "%d/%d\nP=%s" % (int(s["n_pos"]), len(v), ("%.3f" % hp) if hp >= 0.001 else "%.1e" % hp),
            ha="center", va="bottom", fontsize=4.6, linespacing=1.0)
    facts.append(("Fig2b", t, "n %d pos %s median %.3f Holm %.3g" % (len(v), s["n_pos"], fnum(s["median"]), hp), ""))
ax.axhline(0, color=NULL, lw=0.5, ls=":")
ax.set_xticks(range(len(P6))); ax.set_xticklabels([SHORT[t] for t in P6], rotation=45, ha="right")
ax.set_ylim(-2.7, 4.1); ax.set_yticks([-2, -1, 0, 1, 2]); ax.set_ylabel("Deep − shallow contrast (s.d.)")
ax.set_title("65 disjoint deep origins", fontsize=6.5, pad=3)

# c: clade-centred f_T vs depth
met = {r["code"]: r for r in rtsv(W("physchem", "species_metrics.tsv"))}
spT = rtsv(W("final_v3", "pgls_controls", "species_T.tsv"))
rows = [(r["clade"], fnum(r["ldepth"]), fnum(met[r["code"]]["busco_f_T"]), r["L"]) for r in spT if r["code"] in met]
rows = [r for r in rows if not (math.isnan(r[1]) or math.isnan(r[2]))]
byc = {}
for c, d, f, L in rows:
    byc.setdefault(c, []).append(f)
mu = {c: (np.mean(v), np.std(v) or 1.0) for c, v in byc.items()}
x = np.array([r[1] for r in rows]); yv = np.array([(r[2] - mu[r[0]][0]) / mu[r[0]][1] for r in rows])
isd = np.array([r[3] == "D" for r in rows])
ax = fig.add_subplot(gs[0, 2]); panel(ax, "c", -0.3)
ax.scatter(x[~isd], yv[~isd], s=2.5, color=SHAL, alpha=0.7, lw=0, label="Shallow")
ax.scatter(x[isd], yv[isd], s=2.5, color=DEEP, alpha=0.7, lw=0, label="Deep")
edges = np.arange(0, 4.01, 0.5); cx, cy = [], []
for a, b in zip(edges[:-1], edges[1:]):
    m = (x >= a) & (x < b)
    if m.sum() >= 8:
        cx.append((a + b) / 2); cy.append(np.mean(yv[m]))
ax.plot(cx, cy, color=INK, lw=0.9, marker="o", ms=2)
ax.axvline(math.log10(200), color=NULL, lw=0.5, ls=":")
ax.set_xticks([0, 1, 2, 3, 4]); ax.set_xticklabels(["1", "10", "100", "1,000", "10,000"])
ax.set_xlabel("Median depth (m)"); ax.set_ylabel("$f_T$, clade-centred (s.d.)")
ax.legend(loc="upper left", handletextpad=0.1, borderaxespad=0.2, markerscale=2)
facts.append(("Fig2c", "n species", str(len(rows)), "binned means " + ", ".join("%.2f:%.3f" % p for p in zip(cx, cy))))

# d: gene-level sign fraction vs calibrated null
gn = {r["trait"]: r for r in rtsv(W("addon", "A1", "gene_null_summary.tsv"))}
ax = fig.add_subplot(gs[1, 1:]); panel(ax, "d", -0.13)
for i, t in enumerate(P6):
    r = gn[GMAP[t]]
    ax.plot([i, i], [fnum(r["f_sim_lo"]), fnum(r["f_sim_hi"])], color=NULL, lw=3, solid_capstyle="butt", alpha=0.6)
    ax.plot(i, fnum(r["f_sim_mean"]), "_", color=NULL, ms=6)
    col = DEEP if fnum(r["p_cal"]) < 0.05 else INK
    ax.plot(i, fnum(r["f_obs"]), "o", color=col, ms=3.2)
    ax.text(i + 0.12, fnum(r["f_obs"]), ("P=%.2f" % fnum(r["p_cal"])) if fnum(r["p_cal"]) > 0.0101 else "P=0.0099 (min.)", fontsize=5, va="center")
    facts.append(("Fig2d", t, "f_obs %.3f null %.3f [%.3f, %.3f] p_cal %.3f" % tuple(fnum(r[k]) for k in ("f_obs", "f_sim_mean", "f_sim_lo", "f_sim_hi", "p_cal")), ""))
ax.axhline(0.5, color=LIGHT, lw=0.5)
ax.set_xticks(range(len(P6))); ax.set_xticklabels([SHORT[t] for t in P6]); ax.set_xlim(-0.5, len(P6) - 0.2)
ax.set_ylabel("Loci with positive deep effect")
ax.set_title("666 loci; grey, 95% range of 100 simulated label sets", fontsize=6, pad=4)

save(fig, "figures", "Fig2")
with open(os.path.join(BUILD, "fig2_facts.tsv"), "w", encoding="utf-8") as f:
    for r in facts:
        f.write("\t".join(r) + "\n")
