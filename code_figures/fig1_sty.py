"""Figure 1: sampling design and repeated deep-sea colonisation (dated tree, labels, transitions, origin ages)."""
import math, statistics as st
from collections import Counter
import numpy as np
from figcommon_sty import *

tree = parse_newick(open(W("final_v3", "date", "lsd.timetree.nwk")).read())
tips, internal = ape_numbering(tree)
lab = {r["code"]: r for r in rtsv(os.path.join(DATA, "meta", "labels_v2.tsv"))}
# depth from root; scale to Ma with the LSD2 root age (833.0 Ma, lsd.iqtree / node_ages.tsv)
for n in postorder(tree)[::-1]:
    n.x = 0.0 if n.parent is None else n.parent.x + n.length
H = max(t.x for t in tips); SC = 833.0 / H
for k, t in enumerate(tips):
    t.y = k
for n in postorder(tree):
    if n.children:
        n.y = (n.children[0].y + n.children[-1].y) / 2
        n.tips = set().union(*[c.tips for c in n.children])
    else:
        n.tips = {n.name}
N = len(tips)
ang = lambda y: math.pi / 2 - 2 * math.pi * (y + 0.5) / N * 0.97
rad = lambda x: x * SC

fig = plt.figure(figsize=(180 * MM, 118 * MM))
ax = fig.add_axes([0.0, 0.02, 0.60, 0.94]); ax.set_aspect("equal"); ax.axis("off")
R = 833.0
for ma in (200, 400, 600):
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot((R - ma) * np.cos(th), (R - ma) * np.sin(th), color=LIGHT, lw=0.4, ls=(0, (2, 2)), zorder=0)
    ax.text(0, R - ma + 6, "%d Ma" % ma, ha="center", va="bottom", fontsize=5, color=NULL)


def state(n):
    s = {lab[t]["lab"] for t in n.tips}
    return DEEP if s == {"D"} else SHAL if s == {"S"} else NULL


for n in postorder(tree):
    if n.parent is None:
        continue
    a, r0, r1 = ang(n.y), rad(n.parent.x), rad(n.x)
    ax.plot([r0 * math.cos(a), r1 * math.cos(a)], [r0 * math.sin(a), r1 * math.sin(a)], color=state(n), lw=0.35,
            solid_capstyle="butt")
    if n.children:
        a0, a1 = ang(n.children[0].y), ang(n.children[-1].y)
        th = np.linspace(a0, a1, max(2, int(abs(a1 - a0) * 60)))
        ax.plot(r1 * np.cos(th), r1 * np.sin(th), color=NULL if state(n) == NULL else state(n), lw=0.35)
# tip label ring
for t in tips:
    L = lab[t.name]["lab"]; a = ang(t.y)
    c = DEEP if L == "D" else SHAL if L == "S" else "#ffffff"
    ax.plot([(R + 8) * math.cos(a), (R + 48) * math.cos(a)], [(R + 8) * math.sin(a), (R + 48) * math.sin(a)],
            color=c, lw=1.0, solid_capstyle="butt")
# deep-clade crown nodes (deep_origins.tsv, 97 maximal deep clades)
byset = {frozenset(n.tips): n for n in postorder(tree)}
ncrown = 0
for r in rtsv(W("final_v3", "asr", "deep_origins.tsv")):
    s = frozenset(r["tips"].split(","))
    n = byset.get(s)
    if n is None:
        continue
    ncrown += 1
    a, rr = ang(n.y), rad(n.x)
    ax.plot(rr * math.cos(a), rr * math.sin(a), "o", ms=1.6, mfc="white", mec=DEEP, mew=0.4, zorder=5)
# clade arcs
FOCAL = ["Actinopteri", "Chondrichthyes", "Mollusca", "Annelida", "Arthropoda", "Echinodermata", "Cnidaria", "Porifera"]
for cl in FOCAL:
    ys = [t.y for t in tips if lab[t.name]["clade"] == cl]
    a0, a1 = ang(min(ys)), ang(max(ys)); th = np.linspace(a0, a1, 100); rr = R + 70
    ax.plot(rr * np.cos(th), rr * np.sin(th), color=INK, lw=1.2, solid_capstyle="butt")
    am = (a0 + a1) / 2; rt = R + 95
    ha = "left" if math.cos(am) > 0.15 else "right" if math.cos(am) < -0.15 else "center"
    ax.text(rt * math.cos(am), rt * math.sin(am), cl, ha=ha, va="center", fontsize=6.5, fontstyle="italic"
            if False else "normal")
ax.set_xlim(-R - 330, R + 330); ax.set_ylim(-R - 180, R + 180)
ax.text(-R - 320, R + 165, "a", fontsize=8, fontweight="bold")
ax.plot([], [], color=DEEP, lw=2, label="Deep (>200 m), n = %d" % sum(lab[t.name]["lab"] == "D" for t in tips))
ax.plot([], [], color=SHAL, lw=2, label="Shallow, n = %d" % sum(lab[t.name]["lab"] == "S" for t in tips))
ax.plot([], [], "o", ms=3, mfc="white", mec=DEEP, mew=0.6, ls="none", label="Crown of a deep clade (%d)" % ncrown)
ax.legend(loc="lower left", bbox_to_anchor=(-0.02, -0.04), handlelength=1.2, fontsize=5.5)
ax.set_ylim(-R - 330, R + 180)

# b sampling per focal clade
axb = fig.add_axes([0.70, 0.62, 0.27, 0.32])
cnt = {cl: Counter(lab[t.name]["lab"] for t in tips if lab[t.name]["clade"] == cl) for cl in FOCAL}
yy = np.arange(len(FOCAL))[::-1]
axb.barh(yy, [cnt[c]["D"] for c in FOCAL], color=DEEP, height=0.7, label="Deep")
axb.barh(yy, [-cnt[c]["S"] for c in FOCAL], color=SHAL, height=0.7, label="Shallow")
axb.set_yticks(yy); axb.set_yticklabels(FOCAL); axb.axvline(0, color=INK, lw=0.5)
axb.set_xticks([-50, 0, 50, 100]); axb.set_xticklabels(["50", "0", "50", "100"]); axb.set_xlabel("Species (shallow | deep)")
for y, c in zip(yy, FOCAL):
    axb.text(cnt[c]["D"] + 2, y, str(cnt[c]["D"]), va="center", fontsize=5.5)
    axb.text(-cnt[c]["S"] - 2, y, str(cnt[c]["S"]), va="center", ha="right", fontsize=5.5)
axb.set_xlim(-75, 120); panel(axb, "b", x=-0.45)

# c stochastic-map transitions
tr = rtsv(W("final_v3", "asr", "transitions.tsv"))
sd = np.array([int(r["S,D"]) for r in tr]); ds = np.array([int(r["D,S"]) for r in tr])
axc = fig.add_axes([0.70, 0.12, 0.12, 0.33])
bins = np.arange(min(sd.min(), ds.min()) - 2, max(sd.max(), ds.max()) + 4, 4)
axc.hist(sd, bins=bins, color=DEEP, alpha=0.85, label="Shallow→deep")
axc.hist(ds, bins=bins, color=SHAL, alpha=0.7, label="Deep→shallow")
axc.set_xlabel("Transitions per map"); axc.set_ylabel("Stochastic maps (of %d)" % len(sd))
axc.legend(loc="upper right", bbox_to_anchor=(1.12, 1.0), handlelength=0.6, fontsize=5)
panel(axc, "c", x=-0.55, y=1.14)

# d stem ages of the disjoint deep origins (A2)
org = [r for r in rtsv(W("addon", "A2", "origins.tsv")) if r["set"] == "disjoint"]
ages = np.array([fnum(r["stem_age"]) for r in org])
axd = fig.add_axes([0.88, 0.12, 0.10, 0.33])
yv = np.sort(ages)
axd.plot(np.ones_like(yv) + np.random.default_rng(1).uniform(-0.25, 0.25, len(yv)), yv, "o", ms=1.8, mfc=DEEP, mec="none", alpha=0.8)
axd.set_yscale("log"); axd.set_xlim(0.4, 1.6); axd.set_xticks([1]); axd.set_xticklabels(["%d origins" % len(yv)])
axd.set_ylabel("Stem age (Ma)"); axd.axhline(np.median(yv), color=INK, lw=0.6)
panel(axd, "d", x=-0.9, y=1.14)
save(fig, "figures", "Fig1")
with open(os.path.join(BUILD, "fig1_facts.tsv"), "w") as f:
    f.write("tips\t%d\ncrowns_drawn\t%d\nS_D_mean\t%.1f\nD_S_mean\t%.1f\norigins\t%d\norigin_age_median\t%.1f\norigin_age_range\t%.2f-%.1f\n" % (
        N, ncrown, sd.mean(), ds.mean(), len(ages), np.median(ages), ages.min(), ages.max()))
print(open(os.path.join(BUILD, "fig1_facts.tsv")).read())
