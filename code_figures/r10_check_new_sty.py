# r10 new-number check: every number new in r10 (and the r9 new numbers) must appear in main.md or the SI AND its
# source token in the source file. Usage (from paper/): python scripts/r10_check_new_sty.py submission_ready_20261107
import sys, os
PKG = sys.argv[1]
M = open(os.path.join(PKG, "manuscript", "main.md"), encoding="utf-8").read()
SI = open(os.path.join(PKG, "supplementary", "Supplementary_Information.md"), encoding="utf-8").read()
S = "data_srv/"
FA = S + "work/addon/R6/r6_1_failure_audit_20261105.txt"
QV = S + "work/addon/G3run/g3r6_quality_v2_20261105.txt"
RC = S + "work/func_v2/combined_ratecond/recurrence_ratecond.tsv"
RS = S + "work/func_v2/combined/recurrence_summary.tsv"
R7 = S + "addon_R6/R6_7/r6_7_report.txt"
R6 = S + "addon_R6/R6_6/r6_6_report.txt"
DV = S + "addon_R6/R6_7/DEVIATIONS_r6.md"
CV = S + "work/addon/R4_direction/og_coverage_summary.tsv"
PW2 = S + "work/addon/R2/main/recurrence_power.tsv"
PW5 = S + "work/addon/R5/D2/recurrence_power.tsv"
F1 = S + "work/addon/F1_modules/out/F1_report.txt"
F4 = "../viz_v4_sty/out/fig4_facts_v4.tsv"
F4S = os.path.join(PKG, "qa_v2", "fig4_source_v5.tsv")
F3S = "../viz_v4_sty/out/fig3_source_v4.tsv"
F2F = os.path.join(PKG, "qa_v2", "fig2_facts_v5.tsv")
C = [  # r9 items
    ("137 failed HOGs", "137 (14.0%) had at least one failed fit", FA, ">=1 arm failed 137; total 981", "SI"),
    ("fail per arm", "5.1–7.5% per arm", FA, "50 / 981 (5.1%)", "M"),
    ("McNemar 62/43", "62 versus 43 orthogroups", FA, "real-only failure 62, null1-only failure 43", "M"),
    ("McNemar P", "P = 0.079", FA, "exact McNemar p 0.0785", "M"),
    ("grid min lower", "+0.72", FA, "+0.94/+0.72", "M"),
    ("repl 758/776", "758 of 776 fits (97.7%)", DV, "A = 758/776 = 0.977", "M"),
    ("repl logOR", "1.135 to 1.215", DV, "is 1.135 original and 1.215 replicate", "M"),
    ("sat excl top", "1.30, 1.16–1.43", R7, "excluding top 1.298 (1.160-1.434)", "M"),
    ("sat all", "against 1.10 in all", R7, "real log OR all 1.101", "M"),
    ("sat adjusted", "gave 1.27 (1.13–1.40)", R7, "1.265 (1.132-1.402)", "M"),
    ("sat top quartile", "(0.42, 0.19–0.65)", R7, "top quartile 0.421 (0.190-0.649)", "M"),
    ("sat paired ratio", "median paired ratio 0.65", R7, "median ratio 0.645", "M"),
    ("ratecond null mean", "71.1 expected", RC, "relax\tP_clade\t2\t73\t71.1326", "M"),
    ("ratecond P", "P = 0.42", RC, "0.415658", "M"),
    ("chrom-level blind", "−0.125, −0.198 to −0.073", QV, "est -0.1249 [-0.1977, -0.0733] clades 7 lin 33", "M"),
    ("R6-6 decision", "did not predict", R6, "DECISION (ADDENDUM_R6_6): NOT SUPPORTED", "SI"),
    ("R6-7 fits", "7,992 real and 3,198 null", R7, "real 7992  null1 3198", "SI"),
    # r10 items
    # r10c (v5 Fig. 2a,b): share of significant deep fits at K < 1, paired orthogroup differences
    ("sig K<1 Cnidaria", "97–98% of significant fits on deep branches had K < 1", F2F, "530 of 548 significant deep fits at K < 1 (96.7%)", "M"),
    ("sig K<1 Mollusca", "97–98% of significant fits", F2F, "641 of 652 significant deep fits at K < 1 (98.3%)", "M"),
    ("sig K<1 Echinodermata", "97–98% of significant fits", F2F, "430 of 443 significant deep fits at K < 1 (97.1%)", "M"),
    ("paired Echinodermata", "median paired differences of −0.11, −0.09 and −0.08 log10 units in echinoderms, molluscs and cnidarians", F2F, "paired_dlog10K_Echinodermata\tn 400 median -0.109", "M"),
    ("paired Mollusca", "−0.09 and −0.08 log10 units", F2F, "paired_dlog10K_Mollusca\tn 400 median -0.093", "M"),
    ("paired Cnidaria", "−0.08 log10 units in echinoderms", F2F, "paired_dlog10K_Cnidaria\tn 400 median -0.078", "M"),
    ("paired Annelida", "and −0.04 in annelids", F2F, "paired_dlog10K_Annelida\tn 399 median -0.039", "M"),
    ("paired others low", "against −0.02 to 0.01 in the other four clades", F2F, "paired_dlog10K_Arthropoda\tn 399 median -0.024", "M"),
    ("paired others high", "against −0.02 to 0.01", F2F, "paired_dlog10K_Chondrichthyes\tn 398 median 0.013", "M"),
    ("paired others Porifera", "against −0.02 to 0.01", F2F, "paired_dlog10K_Porifera\tn 400 median 0.010", "M"),
    ("paired others Actinopteri", "against −0.02 to 0.01", F2F, "paired_dlog10K_Actinopteri_m\tn 400 median -0.018", "M"),
    ("paired n", "398–400 orthogroups per clade", F2F, "paired_dlog10K_Chondrichthyes\tn 398", "M"),
    ("coverage >=2", "333 were tested in two or more clades", CV, "relax\t2\t297\t333", "M"),
    ("coverage >=3", "36 in three or more", CV, "relax\t3\t32\t36", "M"),
    ("pairs", "22 clade pairs and two larger combinations", F4S, "Porifera+Annelida+Mollusca+Chondrichthyes", "M"),
    ("prereg k2", "52.6 expected, P = 1 × 10^−4^", RS, "relax\t2\t73\t52.5525\t63\t9.999e-05", "M"),
    ("relax k3", "(P = 0.91 and 0.89", RC, "relax\tP_clade\t3\t2\t3.6041\t7\t0.894811", "M"),
    ("relax k3 prereg", "P = 0.91", RS, "relax\t3\t2\t2.8802\t5\t0.906709", "M"),
    ("bph k3", "2 in three or more (P = 0.71", RS, "bph_strict\t3\t2\t2.0086\t4\t0.714829", "M"),
    ("power K50 k4", "with probability 0.92", PW2, "K50\t4\t0.916", "M"),
    ("power P1 k4", "with probability 0.97", PW2, "P1\t4\t0.975", "M"),
    ("power K75 k4", "the probability was 0.31", PW2, "K75\t4\t0.309", "M"),
    ("power K50h k4", "below 0.6", PW2, "K50h\t4\t0.567", "M"),
    ("power K85 k4", "below 0.6", PW5, "K85\t4\t0.060", "M"),
    ("cells", "1,344 calls over 3,091", "../viz_v4_sty/relax_og_perclade_srv.tsv", "Mollusca\t377\t242", "M"),
    ("F1 M1 s1", "θ was 0.006 in the first sample (P = 0.57", F1, "M1_membrane_ion testable clades 8, theta 0.0057, P_lower 0.5709", "SI"),
    ("F1 M1 s2", "−0.046 in the blind sample (P = 0.10)", F1, "M1_membrane_ion testable clades 8, theta -0.0456, P_lower 0.1009", "SI"),
    ("F1 M2 s1", "θ was 0.044 in the first sample (P = 0.86", F1, "theta 0.0442, P_lower 0.8644", "SI"),
    ("F1 M2 s2", "−0.105 in the blind sample (P = 0.012; five of six", F1, "testable clades 6, theta -0.1047, P_lower 0.0121, P_two 0.0242, clades < 0: 5", "SI"),
    ("F1 light", "θ = 0.065 (P = 0.92)", F1, "C_light testable clades 7, theta 0.0652, P_lower 0.9222", "SI"),
    ("F1 M1 OGs", "of 111 metazoan orthologous groups with a relaxed M1 member, 107", F1, "OGs 111 in 1 clade 107", "SI"),
    ("F1 M2 OGs", "of 71 for M2, 66", F1, "OGs 71 in 1 clade 66", "SI"),
    ("F1 decision", "No module met the rule", F1, "M2_proteostasis: NOT SUPPORTED", "SI"),
]
_r = [l.split("	") for l in open(F4S, encoding="utf-8") if l.startswith("a	RELAX relaxed	")]
_pairs = sum(1 for x in _r if int(x[3]) > 0 and int(x[6]) == 2); _big = sum(1 for x in _r if int(x[3]) > 0 and int(x[6]) > 2)
assert (_pairs, _big) == (22, 2) and sum(int(x[3]) for x in _r) == 73, (_pairs, _big)
src, bad = {}, 0
out = ["| Item | Text | Where | Source | Status |", "|---|---|---|---|---|"]
for lab, mt, f, tok, where in C:
    if f not in src:
        src[f] = open(f, encoding="utf-8").read()
    txt = M if where == "M" else M + SI
    ok = (mt in txt) and (tok in src[f])
    bad += not ok
    out.append(f"| {lab} | {mt} | {where} | {f} | {'OK' if ok else 'MISSING'} |")
open(os.path.join(PKG, "QC_number_check_r10_new.md"), "w", encoding="utf-8").write("\n".join(out) + f"\n\n{len(C)} checks, {bad} missing.\n")
print(len(C), "checks,", bad, "missing")
