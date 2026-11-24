"""r10c (2026-11-07): text for the v5 Figs 2-4 (viz_v5_sty) in submission_ready_20261107/manuscript/main.md.
Patches the r10b main.md in place (md5 asserted, so it runs once): legends of Figs 2-4 and three Results sentences
whose panel references or panel contents changed. Every number is from viz_v5_sty/out/fig2_facts_v5.tsv,
fig3_source_v5.tsv or the r10 text; no statistic is new except the Fig. 2a share of significant deep fits at K < 1
and the Fig. 2b paired-difference medians (fig2_facts_v5.tsv).
Usage (from paper/): python scripts/r10c_fig_text_sty.py submission_ready_20261107"""
import os, sys, hashlib

PKG = sys.argv[1]
MD = os.path.join(PKG, "manuscript", "main.md")
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5(MD).startswith("125a749b"), md5(MD)
T = open(MD, encoding="utf-8").read()

REP = [
    ("Median K on deep branches fell to 0.80–0.83 in echinoderms, cnidarians and molluscs and to 0.94 in annelids, "
     "whereas null medians stayed at 0.98–1.06 in every clade (Fig. 2a).",
     "In cnidarians, molluscs and echinoderms, 97–98% of significant fits on deep branches had K < 1 (Fig. 2a), and "
     "within orthogroups K was lower on the deep than on the matched null branch, with median paired differences of "
     "−0.11, −0.09 and −0.08 log10 units in echinoderms, molluscs and cnidarians and −0.04 in annelids, against −0.02 "
     "to 0.01 in the other four clades (Fig. 2b)."),
    ("and that of intensification −0.38 (−0.51 to −0.26; Fig. 2b,c).",
     "and that of intensification −0.38 (−0.51 to −0.26; Fig. 2c)."),
    ("the pooled effect was −0.071 (−0.091 to −0.049; Fig. 3d).",
     "the pooled effect was −0.071 (−0.091 to −0.049; Fig. 3c,d)."),
    ("Left bars, orthologues relaxed in each clade (1,344 calls over 3,091 tested clade-by-orthologue cells).",
     "Left bars, orthologues relaxed in each clade, in the clade colours of Figs 2 and 3 (1,344 calls over 3,091 "
     "tested clade-by-orthologue cells)."),
]

FIG2 = ("**Fig. 2 | Relaxed selection is more frequent on deep than on matched shallow branches.** a, RELAX fits of 1,000 "
        "orthogroups per clade on the deep foreground (top) and of 400 of them on a matched shallow null foreground of "
        "the same trees (bottom): selection-intensity parameter K against −log10 P. Fits with P < 0.05 are coloured by "
        "clade and others are grey; K ≤ 0.001 and P ≤ 10^−10^ are drawn at the axis limits. Numbers give the "
        "percentage of orthogroups called relaxed (P < 0.05 and K < 1; left) and intensified (P < 0.05 and K > 1; "
        "right). b, Paired difference in log10 K between the deep and the matched null fit of each orthogroup (398–400 "
        "orthogroups per clade) beside the LSD2 time tree of the eight clades; densities are drawn between −1.5 and "
        "1.5, bars mark clade medians and the dashed line marks no difference. c, Log odds ratio of a relaxed (left) or "
        "intensified (right) call on deep versus null branches; All clades, clade-adjusted estimate from logistic "
        "regression with clade as a covariate. Error bars, 95% confidence intervals from 2,000 resamples of "
        "orthogroups within clade, keeping the deep and null fits of an orthogroup together; numbers, estimates with "
        "two-sided bootstrap P (***P < 0.01, **P < 0.05, *P < 0.1). Fits without P or K (eight deep, two null) are "
        "counted as neither relaxed nor intensified. Clade colours are the same in Figs 2–4. The preregistered "
        "permutation combination over five matched nulls is in Extended Data Fig. 5.")
FIG3 = ("**Fig. 3 | Relaxation reflects weaker purifying selection and recurs across most deep lineages.** a, b, Change "
        "in BUSTED-PH ω classes from background to foreground branches in the 400 orthogroups per clade with deep "
        "(blue) and matched null (grey) fits: a, log10 proportion-weighted ω of purifying classes (Δ_pur); b, "
        "proportion of sites in positive classes (Δ_pos). Boxes, median and interquartile range; whiskers, 1.5 × "
        "interquartile range; points, outlying orthogroups; labels, two-sided bootstrap P of the paired excess of deep "
        "over null (***P < 0.01). c, Lineage effects (median log10 K of a deep lineage minus the median of the null "
        "pseudo-lineages of its clade, from group-mode RELAX) of the 65 deep lineages of the first sample, the 81 of a "
        "blind sample of new orthogroups and the 83 null pseudo-lineages of the first sample, each present in at "
        "least five orthogroups. Points, lineages; black point and line, median and interquartile range; triangles, "
        "lineages beyond the axis; numbers, percentage below zero. d, Pooled lineage effect (mean of clade means, "
        "clades with at least two lineages) under the robustness settings in the first and blind samples; error bars, "
        "95% intervals from 2,000 resamples of orthogroups within clade. e, The 81 blind-sample lineages on the LSD2 "
        "tree, one species shown per lineage (diamonds, crowns of the sampled lineages of each clade), with depth class "
        "(geometric mean over the lineage's species of median occurrence depth, or of specimen depth where no "
        "occurrence depth exists), numbers of species and of orthogroups in each sample (circle area scales with the "
        "number) and lineage effects in each sample (bars; black lines, clade means; numbers, lineages below zero of "
        "lineages summarised). The Chondrichthyes set includes one hagfish lineage. Bold names, lineages with a vent, "
        "seep or fall taxon or a chemosymbiotic host.")

for a, b in REP:
    assert T.count(a) == 1, (T.count(a), a[:70])
    T = T.replace(a, b)
B = T.split("\n\n")
for k, new in ((2, FIG2), (3, FIG3)):
    j = [i for i, b in enumerate(B) if b.startswith("**Fig. %d |" % k)]
    assert len(j) == 1, k
    B[j[0]] = new
T = "\n\n".join(B)
assert "Fig. 2b,c" not in T.split("## Figure legends")[0]
open(MD, "w", encoding="utf-8", newline="\n").write(T)
wc = lambda s: len(s.split())
print("main.md", md5(MD)[:8], "| legend words Fig2 %d Fig3 %d" % (wc(FIG2), wc(FIG3)))
