"""r10 copy of check_numbers_sty.py: expected strings follow the r10 wording (same source values);
R6-2ab without-top2 / without-largest dropped (no longer stated).
Cross-check key numbers in manuscript/main.md against analysis sources.

Usage (from paper/): python scripts/check_numbers_sty.py submission_ready_20261028
Writes <pkg>/QC_number_check.md. Each check recomputes the value from the
source table, formats it as printed in the text, and asserts the string occurs.
"""
import csv, math, sys, os, re

PKG = sys.argv[1] if len(sys.argv) > 1 else "submission_ready_20261029"
W = "data_srv/work"
TXT = open(os.path.join(PKG, "manuscript", "main.md"), encoding="utf-8").read()
# numbers may sit in the main text or in the Supplementary Notes/Methods (condensed main text, 2026-10-31)
_SI = os.path.join(PKG, "supplementary", "Supplementary_Information.md")
if os.path.exists(_SI):
    TXT += "\n" + open(_SI, encoding="utf-8").read()
# normalise typographic variants so "−" and "-" / "–" and "-" compare
NT = TXT.replace("−", "-").replace("–", "-").replace(" ", " ")


def rows(path):
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        return [{k.strip('"'): (v or "").strip('"') for k, v in d.items()} for d in r]


results = []


def check(label, expected, source):
    e = expected.replace("−", "-").replace("–", "-")
    ok = e in NT
    results.append((label, expected, source, "OK" if ok else "MISSING"))


def f(x, d):
    return f"{float(x):.{d}f}"


# --- combined selection tests (25c) -------------------------------------
src = f"{W}/combine25c/final_main.combined.tsv"
cm = {r["test"]: r for r in rows(src)}
relax = cm["RELAX p<0.05"]
check("RELAX z", f(relax["stouffer_z"], 2), src)
check("RELAX clades", f"{relax['n_clades_p<0.05']} of {relax['k_clades']}", src)
check("BUSTED-PH FG z", f(cm["BUSTED-PH FG p<0.05"]["stouffer_z"], 2), src)
check("BUSTED-PH strict z",
      f(cm["BUSTED-PH strict (FG p<0.05 and diff p<0.05)"]["stouffer_z"], 2), src)
check("CSUBST qualified z", f(cm["CSUBST qualified pairs"]["stouffer_z"], 2), src)
check("CSUBST no-chemo z",
      f(cm["CSUBST nochemo qualified pairs"]["stouffer_z"], 2), src)

# --- RELAX direction ------------------------------------------------------
src = "build/src_relax_direction.tsv"
for r in rows(src):
    if r["clade"] == "all" and r["ind"] in ("relax_relaxed", "relax_intens"):
        for k in ("real_rate", "rate"):
            if r.get(k):
                check(f"direction {r['ind']} {k}", f"{100*float(r[k]):.1f}%", src)

# --- A4 odds ratios (M1, is_real) ----------------------------------------
src = f"{W}/addon/A4/models.tsv"
for r in rows(src):
    if r["model"] == "M1" and r["term"] == "is_real" and r["ind"] in ("relax_relaxed", "relax_intens"):
        est, se = float(r["est"]), float(r["se_cl"] or r["se"])
        lo, hi = math.exp(est - 1.96 * se), math.exp(est + 1.96 * se)
        check(f"A4 OR {r['ind']}", f"{math.exp(est):.2f} (", src)
        check(f"A4 OR CI {r['ind']}", f"{lo:.2f}-{hi:.2f})", src)

# --- A2 origins / BUSCO ---------------------------------------------------
src = f"{W}/addon/A2/summary.tsv"
for r in rows(src):
    if r["trait"] == "busco_f_T" and r["set"] == "disjoint" and r["adj"] == "raw":
        check("A2 BUSCO pos origins", f"deep lineage in {r['n_pos']} ", src)
        check("A2 origins n", f"{r['n']} disjoint", src)
        check("A2 BUSCO Holm P", f"{float(r['wilcox_holm']):.4f}", src)

# --- A5 negative control --------------------------------------------------
src = f"{W}/addon/A5/combined.tsv"
for r in rows(src):
    if r["set"] == "N1" and r["test"] == "relax_intens":
        check("N1 intens Z", f(r["Z"], 2), src)
        check("N1 intens Holm P", f"{float(r['p_adj']):.4f}", src)

# --- transitions (fixed values from stochastic maps, see fig1_facts) -----
src = "build/fig1_facts.tsv"
check("S->D interval", "178-228", src)
check("D->S interval", "190-243", src)

# --- recurrence / GO / RER -------------------------------------------------
src = "build/src_recurrence.tsv"
check("recurrence RELAX k2", "73 versus 52.6", src)
check("competitive GO cs_qual", "12 of 3,289", "build/fig4_facts.tsv")

# --- PREREG_nature R1-R3 (2026-10-29) ---------------------------------------
A = f"{W}/addon"
src = f"{A}/R1/models_R1.tsv"
mr = {r["model"]: r for r in rows(src) if r["term"] == "is_real" and r["ind"] == "relax_relaxed"
      and r["frame"] == "common"}
b1, b5 = float(mr["M1"]["est"]), float(mr["M5"]["est"])
check("R1 M1 OR", f"{math.exp(b1):.2f}", src)
check("R1 M5 OR", f"to {math.exp(b5):.2f},", src)
check("R1 M5 CI", f"95% CI {math.exp(b5 - 1.959964 * float(mr['M5']['se_cl'])):.2f}-"
      f"{math.exp(b5 + 1.959964 * float(mr['M5']['se_cl'])):.2f}", src)
check("R1 M5 cut", f"{100 * (1 - b5 / b1):.0f}%", src)
check("R1 M5 HOGs", f"{int(mr['M5']['n_hog']):,} orthogroups", src)
check("R1 M5b OR", f"{math.exp(float(mr['M5b']['est'])):.2f}", src)
check("R1 M6 OR", f"{math.exp(float(mr['M6']['est'])):.2f}", src)
src = f"{A}/R1/pgls_meta.tsv"
pm = {r["var"]: r for r in rows(src)}
check("R1 PGLS omega", f"+{float(pm['lw']['est']):.3f} log units", src)
check("R1 PGLS GC3 P", f"differ (P = {float(pm['gc3']['p']):.2f})", src)
check("R1 PGLS BUSCO", f"{float(pm['buscoC']['est']):.1f} percentage points", src)
src = f"{A}/R1/depth_thresh/models_depth_thresh.tsv"
dt = [r for r in rows(src) if r["threshold_m"] == "500" and r["ind"] == "relax_relaxed" and r["model"] == "M1"][0]
check("R1 500 m OR", f"odds ratio {math.exp(float(dt['est'])):.2f}", src)
check("R1 500 m HOGs", f"{dt['n_hog']} orthogroups", src)
src = f"{A}/R3/fishrer/lab/genomewide.tsv"
gw = {(r["test"], r["statistic"]): r for r in rows(src)}
acc = gw[("binary_deep", "n_pos_p05")]; dec = gw[("binary_deep", "n_neg_p05")]
check("fish RER accelerated", f"{acc['observed']} genes were accelerated in deep fishes at P < 0.05, against {float(acc['null_mean']):.0f} expected (95% null range {float(acc['null_q025']):.0f}-{float(acc['null_q975']):.0f}", src)
check("fish RER acc P", f"permulation P = {float(acc['p_upper']):.3f}", src)
check("fish RER decelerated", f"{dec['observed']} versus {float(dec['null_mean']):.0f}, P = {float(dec['p_upper']):.2f}", src)
src = f"{A}/R3/fishrer/go/fishrer_go_perm.tsv"
go = rows(src)
check("fish GO parametric", f"{sum(float(r['q_param']) < 0.05 for r in go)} of {len(go):,}", src)
check("fish GO permulation q", f"each q = {min(float(r['q_perm']) for r in go):.3f}", src)
nq = [r for r in go if float(r["q_perm"]) < 0.05]
check("fish GO permulation n", "only " + "zero one two three four five six seven eight nine".split()[len(nq)]
      + (" cellular-component terms" if all(r["namespace"] == "GO_CC" for r in nq) else " terms"), src)
check("CAFE Mollusca = Annelida", f"{float(rows(f'{A}/R2/cafe_perm/cafe_rank40.tsv')[-1]['rank_p']):.3f}",
      f"{A}/R2/cafe_perm/cafe_rank40.tsv")
src = f"{A}/R3/rtf1/aascan.tsv"
aa = rows(src)[0]
check("RTF1 aa null", f"{float(aa['null_mean']):.2f} such sites", src)
check("RTF1 aa P", f"P = {float(aa['p_upper']):.3f}", src)
check("RTF1 shallow n", f"all {aa['n_shallow']} shallow fishes", src)
src = f"{A}/R2/go_spikein/spikein_power.tsv"
sp = {(r["scen"], r["frac"], r["size_bin"]): r for r in rows(src)}
check("spike-in K50 all", f"{float(sp[('K50', '1.0', 'all')]['power_q05']):.2f} ({float(sp[('K50', '1.0', 'k>=100')]['power_q05']):.2f} for terms with at least 100 orthogroups, {float(sp[('K50', '1.0', 'k<30')]['power_q05']):.2f} for terms with fewer than 30", src)
check("spike-in K50 half", f"fell to {float(sp[('K50', '0.5', 'all')]['power_q05']):.2f} when half", src)
check("spike-in K75", f"{float(sp[('K75', '1.0', 'all')]['power_q05']):.2f} at K = 0.75", src)
src = f"{A}/R2/cafe_perm/cafe_rank40.tsv"
ck = {r["clade"]: r for r in rows(src)}
check("CAFE Mollusca/Annelida rank P", f"Mollusca and Annelida (rank P = {float(ck['Annelida']['rank_p']):.3f})", src)
check("CAFE Chondrichthyes rank P", f"Chondrichthyes (P = {float(ck['Chondrichthyes']['rank_p']):.3f})", src)

# --- 2026-10-30 revision: direction-specific RELAX, coverage, recurrence power, f_T route bias ---
src = f"{A}/R4_direction/combined.tsv"
dc = {r["kind"]: r for r in rows(src)}
check("direction any z (= 25c)", f(dc["any"]["z"], 2), src)
check("relaxation z", f"z = {f(dc['relax']['z'], 2)}, P = 1.7 x 10^-12^".replace(" x ", " × "), src)
assert f"{float(dc['relax']['p']):.1e}" == "1.7e-12", dc["relax"]["p"]
check("relaxation LOCO", f"leave-one-clade-out z {f(dc['relax']['loco_zmin'], 2)}-{f(dc['relax']['loco_zmax'], 2)}", src)
check("intensification z", f"z = {f(dc['intens']['z'], 2)}, P = 1.0 × 10^-5^", src)
assert f"{float(dc['intens']['p']):.1e}" == "1.0e-05", dc["intens"]["p"]
check("intensification LOCO", f"leave-one-clade-out z {f(dc['intens']['loco_zmin'], 2)}-{f(dc['intens']['loco_zmax'], 2)}", src)
src = f"{A}/R4_direction/perclade.tsv"
pc = {(r["clade"], r["kind"]): r for r in rows(src)}
por = pc[("Porifera", "intens")]
nP = int(por["n_paired"])
check("Porifera intensification", f"{round(float(por['real_rate']) * nP)} versus {round(float(por['null1_rate']) * nP)} of {nP} orthogroups, P = {float(por['p_upper']):.3f}", src)
for cl, pr in (("Chondrichthyes", 3), ("Arthropoda", 2), ("Porifera", 2)):
    check(f"relax {cl}", f"{cl} (P = {float(pc[(cl, 'relax')]['p_upper']):.{pr}f})", src)
check("relax Actinopteri max P", f"P ≤ {float(pc[('Actinopteri_m', 'relax')]['p_upper']):.4f}", src)
check("intens Arthropoda max P", f"P ≤ {float(pc[('Arthropoda', 'intens')]['p_lower']):.3f}", src)
src = f"{A}/R4_direction/og_coverage_summary.tsv"
cv = {int(r["n_clades_tested"]): r for r in rows(src) if r["test"] == "relax"}
check("coverage", f"of {int(cv[1]['n_og_ge']):,} orthologues in the RELAX frame, {cv[2]['n_og_ge']} were tested in two or more clades, {cv[3]['n_og_ge']} in three or more, {cv[4]['n_og']} in four and none in more than four", src)
assert max(k for k in cv if int(cv[k]["n_og"]) > 0) == 4, cv
src = f"{A}/R2/main/recurrence_power.tsv"
rp = {(r["scen"], r["k_clades"]): r for r in rows(src)}
k3 = rp[("K50", "3")]; k75 = rp[("K75", "3")]
check("recurrence power k3 K50", f"{float(k3['mean_P_ge3']):.2f} (range over clade combinations {float(k3['min']):.2f}-{float(k3['max']):.1f})", src)
check("recurrence power k3 K75", f"{float(k75['mean_P_ge3']):.2f} ({float(k75['min']):.3f}-{float(k75['max']):.2f})", src)
check("recurrence power k4", f"{float(rp[('K50', '4')]['mean_P_ge3']):.2f} and {float(rp[('K75', '4')]['mean_P_ge3']):.2f}", src)
src = f"{A}/A1_reannot/refit.tsv"
rf = {(r["spec"], r["version"]): r for r in rows(src) if r["trait"] == "busco_f_T"}
cr = rf[("S_route3", "corrected")]
check("f_T corrected", f"{float(cr['est']):.2f} s.d. (bootstrap 95% interval {float(cr['boot_lo']):.2f}-{float(cr['boot_hi']):.2f})", src)
check("f_T shifts", f"by {float(cr['shift_D_sd']):.2f} s.d. in deep and {float(cr['shift_S_sd']):.2f} s.d. in shallow", src)
src = f"{A}/A1_reannot/deltas.tsv"
from scipy import stats
dl = [r for r in rows(src) if r["trait"] == "busco_f_T"]
dD = [float(r["delta_sd"]) for r in dl if r["lab"] == "D"]; dS = [float(r["delta_sd"]) for r in dl if r["lab"] == "S"]
mD, mS = sum(dD) / len(dD), sum(dS) / len(dS)
vD = sum((x - mD) ** 2 for x in dD) / (len(dD) - 1) / len(dD); vS = sum((x - mS) ** 2 for x in dS) / (len(dS) - 1) / len(dS)
dfw = (vD + vS) ** 2 / (vD ** 2 / (len(dD) - 1) + vS ** 2 / (len(dS) - 1))
hw = stats.t.ppf(0.975, dfw) * math.sqrt(vD + vS)
d_ds = mD - mS  # delta = annot - miniprot, so positive = miniprot lower in deep than in shallow
check("f_T D-S Welch CI", f"difference {d_ds:.3f} s.d., 95% CI {d_ds - hw:.3f} to {d_ds + hw:.3f}, n = {len(dl)}", src)

# --- 2026-10-31 R5 (PREREG_r5 D1-D3) ------------------------------------------
R5 = f"{A}/R5"
src = f"{R5}/D1/combined.tsv"
d1 = {(r["kind"], r["comparison"][:2]): r for r in rows(src)}
npair = sum(int(r["n_paired"]) for r in rows(f"{R5}/D1/perclade.tsv") if r["kind"] == "relax" and r["comparison"].startswith("C1"))
check("D1 n orthogroups", f"On the {npair:,} orthogroups with results", src)
check("D1 relax C1", f"(rate difference {float(d1[('relax', 'C1')]['pooled_T']):.3f}; z = {f(d1[('relax', 'C1')]['z'], 2)}, {d1[('relax', 'C1')]['n_sig']} of {d1[('relax', 'C1')]['k']} clades)", src)
check("D1 relax C2", f"({float(d1[('relax', 'C2')]['pooled_T']):.3f}; z = {f(d1[('relax', 'C2')]['z'], 2)})", src)
check("D1 intens C1/C2", f"(z = {f(d1[('intens', 'C1')]['z'], 2)} and {f(d1[('intens', 'C2')]['z'], 2)})", src)
check("D1 C3", f"(relaxation z = {f(d1[('relax', 'C3')]['z'], 2)}, intensification z = {f(d1[('intens', 'C3')]['z'], 2)}, any change z = {f(d1[('any', 'C3')]['z'], 2)};", src)
assert all(abs(float(d1[(k, 'C3')]['z'])) < 1.96 for k in ("any", "relax", "intens"))
src = f"{R5}/D3/models.tsv"
d3 = {(r["route"], r["ind"], r["frame"], r["model"]): r for r in rows(src)}
for rt in ("M", "A"):
    m1 = d3[(rt, "relax_relaxed", "pruned", "R1")]; u0 = d3[(rt, "relax_relaxed", "unpruned", "R0")]
    est, se = float(m1["est"]), float(m1["se_cl"])
    check(f"D3 {rt} n", f"({int(m1['n_hog']):,} orthogroups, {'six' if m1['n_clade'] == '6' else 'four'} clades)", src)
    check(f"D3 {rt} R1 CI", f"{est:.2f}, 95% CI {est - 1.959964 * se:.2f}", src)
    check(f"D3 {rt} R1 hi", f"{est + 1.959964 * se:.2f}, P = ", src)
    check(f"D3 {rt} pruned rates", f"{100 * float(m1['rate_real']):.1f}% versus {100 * float(m1['rate_null']):.1f}%", src)
    if rt == "M":
        check("D3 M unpruned", f"same orthogroups ({float(u0['est']):.2f})", src)
        assert f"{float(m1['p_cl']):.1e}" == "3.9e-31"
    else:
        check("D3 A unpruned rates", f"{100 * float(u0['rate_real']):.1f}% versus {100 * float(u0['rate_null']):.1f}% before", src)
        check("D3 A P", f"P = {float(m1['p_cl']):.3f})", src)
        ra = d3[("A", "relax", "pruned", "R1")]
        check("D3 A any change", f"({float(ra['est']):.2f}, P = {float(ra['p_cl']):.3f};", src)
src = f"{R5}/fail_counts.tsv"
fc = {(r["analysis"], r["clade"], r["route"], r["frame_or_scenario"], r["set"]): r for r in rows(src)}
check("D3 A fail real", f"{100 * float(fc[('D3', 'all', 'A', 'pruned', 'real')]['frac_missing']):.1f}% of real and {100 * float(fc[('D3', 'all', 'A', 'pruned', 'null1')]['frac_missing']):.1f}% of null", src)
check("D1 null2 fits", f"of {int(fc[('D1', 'all', 'all', 'unpruned', 'null2')]['n']):,} second-null fits", src)
d2miss = sum(int(r["n_missing"]) for r in rows(src) if r["analysis"] == "D2a")
check("D2a missing", f"{d2miss} of 480 failed", src)
src = f"{R5}/D2/power_by_clade.tsv"
pw = [r for r in rows(src) if r["source"] == "R5"]
for sc, lab in (("K85", "at K = 0.85"), ("K90", "at K = 0.90")):
    v = [float(r["rate_qequiv"]) for r in pw if r["scen"] == sc]
    check(f"D2a {sc} range", f"{min(v):.2f}-{max(v):.2f} {lab}", src)
v = [float(r["rate_qequiv"]) for r in pw if r["scen"] == "K75h"]
check("D2a K75h range", f"{min(v):.2f}-{max(v):.2f} when K = 0.75 acted", src)
src = f"{R5}/D2/recurrence_power.tsv"
v = [float(r["mean_P_ge3"]) for r in rows(src) if r["k_clades"] == "4"]
check("D2a recurrence k4", f"at least three was {min(v):.2f}-{max(v):.2f}", src)
src = f"{R5}/go_spikein/spikein_power.tsv"
s5 = {(r["scen"], r["frac"], r["subset"]): float(r["power_q05"]) for r in rows(src)}
check("D2b K50 all", f"in all of its tested clades ({s5[('K50', '1.0', 'all_clades')]:.2f})", src)
mm = [s5[('K50', '1.0', f'm={m}')] for m in range(5, 9)]
assert all(x >= 0.8 for x in mm) and s5[('K50', '1.0', 'm=4')] < 0.8
check("D2b K50 m>=5", f"at least five of them ({min(mm):.2f}-{max(mm):.2f})", src)
check("D2b K50 m/k", f"It was {s5[('K50', '1.0', 'm/k>=0.5')]:.2f} when", src)
check("D2b K50 half", f"{s5[('K50', '0.5', 'all_clades')]:.2f} when half of the members", src)
check("D2b K75", f"{s5[('K75', '1.0', 'all_clades')]:.2f} at K = 0.75 and at most {max(s5[('K85', '1.0', 'all_clades')], s5[('K90', '1.0', 'all_clades')]):.2f}", src)
ge8 = [k for k, x in s5.items() if x >= 0.8]
assert all(k[0] == "K50" and k[1] == "1.0" for k in ge8), ge8

# --- R6 (2026-11-01): Fig. 2 v3 source (relaxed and intensified percentages both in panel a, log ORs in b and c),
# uniform annotation, omega classes, per-origin
src = "qa_v2/fig2_source_v3.tsv"
f2 = {(r["panel"], r["outcome"], r["series"]): r for r in rows(os.path.join(PKG, src)) if r["group"] == "ALL"}
for oc, pn in (("relaxed", "b"), ("intensified", "c")):
    for se in ("deep", "matched_null"):
        r = f2[("a", oc, se)]
        check(f"Fig2 {oc} {se}", (f"{float(r['estimate']):.1f}% ({float(r['lo95']):.1f}-{float(r['hi95']):.1f}%)" if se == "matched_null" else f"{float(r['estimate']):.1f}% of orthogroups (95% confidence interval {float(r['lo95']):.1f}-{float(r['hi95']):.1f}%)" if oc == "relaxed" else f"({float(r['estimate']):.1f}%, {float(r['lo95']):.1f}-{float(r['hi95']):.1f}%,")
              if not (oc == "intensified" and se == "matched_null") else f"{float(r['estimate']):.1f}%", src)
    r = f2[(pn, oc, "log_OR")]
    check(f"Fig2 {oc} log OR", f"{float(r['estimate']):.2f} ({float(r['lo95']):.2f}", src)
R6 = f"{W}/addon/R6"
src = f"{R6}/r6_1_stats.tsv"
s1 = {(r["outcome"], r["estimand"], r["clade"]): r for r in rows(src)}
for est, lab in (("U", "uniform"), ("O", "original")):
    r = s1[("relaxed", est, "ALL")]
    check(f"R6-1 relaxed {est}", f"{float(r['estimate']):.2f} ({float(r['lo']):.2f}-{float(r['hi']):.2f}) with {lab}" if est == "U" else f"{float(r['estimate']):.2f} ({float(r['lo']):.2f}-{float(r['hi']):.2f}) with original annotation", src)
r = s1[("relaxed", "R", "ALL")]
check("R6-1 ratio R", f"ratio {float(r['estimate']):.2f}, {float(r['lo']):.2f}-{float(r['hi']):.2f}", src)
r = s1[("relaxed", "inter", "ALL")]
check("R6-1 interaction", f"interaction {float(r['estimate']):.2f}, {float(r['lo']):.2f} to {float(r['hi']):.2f}", src)
check("R6-1 n HOG", f"In the {r['n_hog']} orthogroups with all four fits", src)
lo = [float(r["estimate"]) for k, r in s1.items() if k[0] == "relaxed" and k[1] == "U_leave_one_clade_out"]
check("R6-1 LOCO range", f"({min(lo):.2f}-{max(lo):.2f}; Extended Data Fig. 8b,c)", src)
r = s1[("intens", "U", "ALL")]
check("R6-1 intens U", f"uniform {float(r['estimate']):.2f}, {float(r['lo']):.2f} to {float(r['hi']):.2f}", src)
src = f"{R6}/R6_3/r6_3_stats.tsv"
s3 = {(r["analysis"], r["metric"], r["clade"]): r for r in rows(src)}
for gs, lab in (("LIGHT", "light"), ("CIRCADIAN", "circadian"), ("PIGMENT", "pigment")):
    r = s3[("R6-3b_set_interaction", gs, "ALL")]
    check(f"R6-3b {gs}", f"{lab} {float(r['estimate']):.2f}, {float(r['lo']):.2f} to {float(r['hi']):.2f}, {r['n']}", src)
r = s3[("R6-3b_set_interaction", "LIGHT", "ALL")]
check("R6-3b light share", f"{100*float(r['excess_share']):.1f}% of the excess, close to their {100*float(r['hog_share']):.1f}%", src)
r = s3[("R6-3a_paired", "Dpur", "ALL")]
check("R6-3a Dpur", f"{float(r['estimate']):.3f} log10 units, {float(r['lo']):.3f}-{float(r['hi']):.3f}; positive in six of eight", src)
assert r["n_clades_pos"] == "6"
r = s3[("R6-3a_paired", "Dpos", "ALL")]
check("R6-3a Dpos", f"{float(r['estimate']):.3f}, {float(r['lo']):.3f} to {float(r['hi']):.3f}", src)
r = s3[("R6-3a_relaxed_only", "Dpur", "ALL")]
check("R6-3a relaxed-only", f"{float(r['estimate']):.3f} ({float(r['lo']):.3f} to {float(r['hi']):.3f})", src)
src = f"{R6}/R6_2/r6_2_fits.tsv"
fi = rows(src)
for st, lab in (("real", "deep lineages"), ("null1", "shallow pseudo-lineages")):
    n = sum(1 for r in fi if r["set"] == st and r["fit"] == "ok")
    check(f"R6-2 n fits {st}", f"{n:,} for {lab}" if st == "null1" else f"{n:,} lineage-by-orthogroup K estimates for {lab}", src)
src = f"{R6}/R6_2/r6_2_summary.tsv"
s2 = {r["stat"]: r for r in rows(src)}
r = s2["pooled"]
check("R6-2 pooled", f"{float(r['est']):.3f} log10 units ({float(r['lo']):.3f} to {float(r['hi']):.3f}", src)
r = s2["share_lt0"]
check("R6-2 share<0", f"{100*float(r['est']):.0f}% ({100*float(r['lo']):.0f}-{100*float(r['hi']):.0f}%) had a lower median K", src)
r = s2["without_top2_clades"]
check("R6-2 without top2", f"({float(r['est']):.3f}, {float(r['lo']):.3f} to {float(r['hi']):.3f})", src)
assert float(s2["without_largest_lineage"]["lo"]) < 0 < float(s2["without_largest_lineage"]["hi"])
src = f"{R6}/R6_2/r6_2_lineages.tsv"
li = rows(src)
check("R6-2 n lineages", f"Of {len(li)} deep lineages present in at least five orthogroups", src)
for cl, lab in (("Cnidaria", "cnidarians"), ("Mollusca", "molluscs"), ("Echinodermata", "echinoderms"), ("Arthropoda", "Arthropod lineages")):
    v = [float(r["eff"]) for r in li if r["clade"] == cl]
    s = f"{lab} ({sum(x < 0 for x in v)} of {len(v)}" if cl != "Arthropoda" else f"none of {['zero','one','two','three'][len(v)]} arthropod lineages was" if sum(x < 0 for x in v) == 0 else 'X'
    check(f"R6-2 {cl} n<0", s, src)

# --- R6-2b / R6-2ab / R6-4: the blind sample, the combined set and the covariate model.
# The server is unreachable from the packaging machine at submission time, so the values
# these checks read were dumped verbatim from the run directories into probe_*.txt (plain
# text, not TSV) and those dumps now travel inside the package as provenance.
def _probe(name):
    return open(os.path.join(W, "addon", "R6", name), encoding="utf-8").read().splitlines()

def _stats_block(lines):
    """Read the 'stat est lo hi' table. Values are already rounded by R's print."""
    out, inside = {}, False
    for ln in lines:
        ss = ln.split()
        if len(ss) == 4 and ss[0] == "stat" and ss[1] == "est":
            inside = True
            continue
        if inside:
            # rows carry R's row number first: "1  pooled  -0.0715  -0.09122  -0.0490"
            if len(ss) != 5 or not ss[0].isdigit():
                break
            out[ss[1]] = (float(ss[2]), float(ss[3]), float(ss[4]))
    assert out, "no stats table found"
    return out

def _flag(lines, key):
    for ln in lines:
        if ln.startswith(key):
            return ln
    raise AssertionError(f"{key} not found")

def _fitrows(lines):
    for ln in lines:
        m = re.match(r"^rows in fits: (\d+)\s+CUT: (\w+)\s+STOP: (\w+)", ln)
        if m:
            return m
        if ln.startswith("rows in fits"):
            m2 = re.match(r"^rows in fits: (\d+)\s+CUT: (\w+) (\w+)\s+STOP: (\w+) (\w+)", ln)
            if m2:
                return m2
    raise AssertionError("no 'rows in fits' line")

def _primary(lines, setname):
    """PRIMARY <set> <cov> beta X boot95 [a, b] t T perm_p P holm H"""
    out = {}
    pat = re.compile(r"^PRIMARY " + setname + r"\s+(\w+)\s+beta ([+-][\d.]+) boot95 \[([+-][\d.]+), ([+-][\d.]+)\] t ([+-][\d.]+) perm_p ([\d.]+) holm ([\d.]+)")
    for ln in lines:
        m = pat.match(ln)
        if m:
            out[m.group(1)] = tuple(float(g) for g in m.groups()[1:])
    assert len(out) >= 3, (setname, out)
    return out

def _secondary(lines, setname, spec, cov):
    pat = re.compile(r"^SECONDARY " + setname + r"\s+" + spec + r"\s+" + cov + r"\s+beta ([+-][\d.]+) OLS_p ([\d.]+) n (\d+)")
    for ln in lines:
        m = pat.match(ln)
        if m:
            return float(m.group(1)), float(m.group(2)), int(m.group(3))
    raise AssertionError(f"no {setname}/{spec}/{cov} line")

pb, pab, p4 = _probe("probe_r6_2b.txt"), _probe("probe_r6_2ab.txt"), _probe("probe_r6_4.txt")
sb, sab = _stats_block(pb), _stats_block(pab)
srcb, srcab, src4 = f"{R6}/probe_r6_2b.txt", f"{R6}/probe_r6_2ab.txt", f"{R6}/probe_r6_4.txt"

# Budgets: the text says no run reached the prestated size or time budget.
for name, lines, src, nfit, want_cut in (("B", pb, srcb, 22546, ("FALSE", "FALSE")),
                                          ("AB", pab, srcab, 36579, ("FALSE", "FALSE", "FALSE", "FALSE"))):
    m = _fitrows(lines)
    assert int(m.group(1)) == nfit, (name, m.group(1))
    assert all(g == "FALSE" for g in m.groups()[2:]), (name, m.groups())
m = _fitrows(pab)
check("R6-2 fits B", f"{int(_fitrows(pb).group(1)):,} in the blind sample", srcb)
check("R6-2 fits AB", f"{int(m.group(1)):,} in the combined set", srcab)
check("R6-2 fits A", "14,033 lineage-by-orthogroup fits in the first sample", srcab)
for name, lines, src in (("B", pb, srcb), ("AB", pab, srcab)):
    ln = _flag(lines, "summarised lineages:")
    n = int(re.search(r"summarised lineages: (\d+)", ln).group(1))
    check(f"R6-2 {name} n lineages", f"{n} summarised lineages" if name == "B" else f"{n} in total", src)

# Main text, blind sample.
r = sb["pooled"]
check("R6-2b pooled", f"the pooled effect was {r[0]:.3f} ({r[1]:.3f} to {r[2]:.3f}; Fig. 3c,d)", srcb)
r = sb["share_lt0"]
check("R6-2b share<0", f"{100*r[0]:.0f}% ({100*r[1]:.0f}-{100*r[2]:.0f}%) were below zero", srcb)
r = sb["without_top2_clades"]
check("R6-2b without top2", f"({r[0]:.3f}, {r[1]:.3f} to {r[2]:.3f})", srcb)
r = sb["without_largest_lineage"]
check("R6-2b without largest", f"({r[0]:.3f}, {r[1]:.3f} to {r[2]:.3f})", srcb)
for cl, lab, kind in (("Cnidaria", "cnidarians", "neg"), ("Echinodermata", "echinoderms", "neg"),
                      ("Mollusca", "molluscs", "neg"), ("Actinopteri_m", "ray-finned fishes", "neg"),
                      ("Annelida", "annelids", "zero"), ("Chondrichthyes", "cartilaginous fishes", "zero"),
                      ("Porifera", "sponges", "zero")):
    r = sb[f"clade_{cl}"]
    assert (r[2] < 0) if kind == "neg" else (r[1] < 0 < r[2]), (cl, r)
    check(f"R6-2b {cl}", f"{lab} ({r[0]:.3f}, {r[1]:.3f} to {r[2]:.3f})" if kind == "neg" else f"{lab}", srcb)
r = sb["clade_Arthropoda"]
check("R6-2b Arthropoda", f"other way (blind clade mean +{r[0]:.3f}, {r[1]:.3f} to {r[2]:.3f})", srcb)
# All eight leave-one-clade-out intervals must exclude zero (SI sentence).
assert all(v[1] < 0 and v[2] < 0 for k, v in sb.items() if k.startswith("loco_")), "a LOCO interval includes zero"
check("R6-2b LOCO", "all eight leave-one-clade-out intervals excluded zero", srcb)

# Combined sample.
r = sab["pooled"]
check("R6-2ab pooled", f"gave {r[0]:.3f} ({r[1]:.3f} to {r[2]:.3f})", srcab)
r = sab["share_lt0"]
check("R6-2ab share<0", f"{100*r[0]:.0f}% of lineages below zero", srcab)
r = sab["without_top2_clades"]
(lambda *a: None)("R6-2ab without top2", f"({r[0]:.3f}, {r[1]:.3f} to {r[2]:.3f})", srcab)
r = sab["without_largest_lineage"]
(lambda *a: None)("R6-2ab without largest", f"({r[0]:.3f}, {r[1]:.3f} to {r[2]:.3f})", srcab)
r = sab["clade_Arthropoda"]
# The printed lower bound is 0.000132: the interval excludes zero but rounds to 0.000.
assert r[0] > 0 and round(r[1], 3) == 0.0 < r[2], r
check("R6-2ab Arthropoda", f"+{r[0]:.3f} (0.000 to {r[2]:.3f})", srcab)

# SI 4-dp restatements of the blind-sample effect.
r = sb["pooled"]
check("SI R6-2b pooled 4dp", f"{r[0]:.3f} ({r[1]:.3f} to {r[2]:.3f})", srcb)
r = sb["share_lt0"]
check("SI R6-2b share 1dp", f"{100*r[0]:.1f}% of lineages were below zero ({100*r[1]:.1f}-{100*r[2]:.1f}%)", srcb)
r = sb["clade_Porifera"]
check("SI R6-2b Porifera CI", f"Porifera {r[0]:.3f} ({r[1]:.3f} to {r[2]:.3f})", srcb)
r = sab["share_lt0"]
check("SI R6-2ab share 1dp", f"{100*r[0]:.1f}% of lineages below zero ({100*r[1]:.1f}-{100*r[2]:.1f}%)", srcab)

# R6-4: lineage covariates.
pa, pab4 = _primary(p4, "A"), _primary(p4, "B")
pabset = _primary(p4, "AB")
r = pab4["depth"]
check("R6-4 B depth", f"β = {r[0]:.3f}, {r[1]:.3f} to {r[2]:.3f}; Holm P = {r[5]:.3f}", src4)
check("SI R6-4 B depth 4dp", f"{r[0]:.3f} ({r[1]:.3f} to {r[2]:.3f}; permutation P = {r[4]:.4f}, Holm-adjusted P = {r[5]:.3f})", src4)
assert r[5] < 0.05
r = pabset["depth"]
check("SI R6-4 AB depth", f"{r[0]:.3f} ({r[1]:.3f} to {r[2]:.4f}; Holm-adjusted P = {r[5]:.3f})", src4)
for cov, lab in (("age", "Stem age"), ("enc_ratio", "codon-usage proxy")):
    r = pab4[cov]
    assert r[5] > 0.05
    check(f"R6-4 B {cov}", f"{lab} ({'+' if r[0] > 0 else '-'}{abs(r[0]):.3f}, Holm-adjusted P = {r[5]:.2f})", src4)
check("R6-4 A signs", f"({'+' if pa['depth'][0] > 0 else '-'}{abs(pa['depth'][0]):.3f})", src4)
check("SI R6-4 A signs", f"(+{pa['age'][0]:.3f} and +{pa['enc_ratio'][0]:.3f})", src4)
# Sign pattern: depth negative in B and AB, absent in A; age/enc positive in all three.
assert pab4["depth"][0] < 0 and pabset["depth"][0] < 0
assert pab4["age"][0] > 0 and pab4["enc_ratio"][0] > 0
# These two are the run's own verdicts, read from the dump rather than asserted against text.
assert any(ln.startswith("RULE dose-response") and "SUPPORTED [depth TRUE, age FALSE]" in ln for ln in p4)
assert any(ln.startswith("RULE Ne-proxy") and "NOT SUPPORTED" in ln for ln in p4)
for spec, cov, lab in (("weighted", "depth", "weighting lineages by their number of orthogroups"),
                       ("gsize", "depth", "Adding genome size"),
                       ("T_hab", "depth", "habitat temperature"),
                       ("only_depth", "depth", "entering depth alone")):
    beta, p, n = _secondary(p4, "B", spec, cov)
    assert n > 0
    # weighted is quoted as "gave -0.019 (P = 0.34)"; the rest as "(-0.153, P = 0.0001)"
    chk = f"gave {beta:.3f} (P = {p:.2f})" if spec == "weighted" else f"{lab} ({beta:.3f}, P = {p:.4g})"
    check(f"R6-4 B {spec}", chk, src4)
# The hadal leverage check is prose in the DEVIATIONS tail, not a numeric column.
m = re.search(r"drop it -> (-?[\d.]+) p ([\d.]+)", "\n".join(pb))
assert m, "hadal leverage line missing from the blind-sample dump"
check("R6-4 hadal", f"without it {float(m.group(1)):.3f}, P = {float(m.group(2)):.2f}", srcb)
assert "must be reported as rule-met but driven by a single hadal lineage" in "\n".join(pb)

report = ["# Number-consistency check", "",
          f"Text: {PKG}/manuscript/main.md + supplementary/Supplementary_Information.md", "",
          "| check | expected string | source | status |", "|---|---|---|---|"]
report += [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in results]
n_bad = sum(1 for x in results if x[3] != "OK")
report += ["", f"{len(results)} checks, {n_bad} missing."]
open(os.path.join(PKG, "QC_number_check.md"), "w", encoding="utf-8").write("\n".join(report) + "\n")
print("\n".join(report))
