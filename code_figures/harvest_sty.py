"""Harvest every number quoted in the manuscript from the server bundle (data_srv).

Writes build/numbers_final.tsv (key, value, source) and build/src_*.tsv tables used by the figure script.
Only reads files; derived statistics (unweighted / leave-one-clade-out Stouffer, transition summaries) are
computed here from the frozen per-clade outputs and labelled as such.
"""
import csv, math, os, re, statistics as st, sys
from collections import Counter, defaultdict
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data_srv")
B = os.path.join(ROOT, "build")
os.makedirs(B, exist_ok=True)
W = lambda *p: os.path.join(D, "work", *p)
CLADES = ["Actinopteri_m", "Mollusca", "Cnidaria", "Annelida", "Echinodermata", "Arthropoda",
          "Chondrichthyes", "Porifera"]
NUM = []


def put(key, value, source):
    NUM.append((key, value, source))


def rtsv(path):
    with open(path, encoding="utf-8") as f:
        rows = [r for r in csv.reader(f, delimiter="\t") if r and not r[0].startswith("#")]
    h = [c.strip('"') for c in rows[0]]
    return [dict(zip(h, [c.strip('"') for c in r])) for r in rows[1:]]


def wtsv(name, rows, cols):
    with open(os.path.join(B, name), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


g = lambda x: float("nan") if x in ("", "NA", "nan", None) else float(x)

# ---- sampling and labels
lab = rtsv(os.path.join(D, "meta", "labels_v2.tsv"))
put("n_species_label_table", len(lab), "meta/labels_v2.tsv rows (all screened species)")
tips = set(re.findall(r"[(,]([A-Za-z0-9_]+):", open(W("final_v3", "species_tree.nwk")).read()))
put("n_tree_tips", len(tips), "work/final_v3/species_tree.nwk")
lab = [r for r in lab if r["code"] in tips]
put("n_species_in_tree_labelled", len(lab), "labels_v2 rows whose code is a tip of species_tree.nwk")
c = Counter(r["lab"] for r in lab)
put("lab_counts", "D %d / S %d / U %d" % (c["D"], c["S"], c["U"]), "meta/labels_v2.tsv lab")
for k in ("lab_strict", "lab_1000"):
    cc = Counter(r[k] for r in lab)
    put(k + "_counts", "D %d / S %d / other %d" % (cc["D"], cc["S"], len(lab) - cc["D"] - cc["S"]),
        "meta/labels_v2.tsv " + k)
put("n_phyla", len(set(r["phylum"] for r in lab)), "meta/labels_v2.tsv phylum")
put("source_routes", "; ".join("%s %d" % kv for kv in Counter(r["source"] for r in lab).most_common()),
    "meta/labels_v2.tsv source")
per = defaultdict(Counter)
for r in lab:
    per[r["clade"]][r["lab"]] += 1
rows = [{"clade": k, "D": v["D"], "S": v["S"], "U": v["U"], "n": sum(v.values())} for k, v in sorted(per.items())]
wtsv("src_labels_by_clade.tsv", rows, ["clade", "D", "S", "U", "n"])
evid = Counter(r["evidence_level"] for r in lab if r["lab"] == "D")
put("deep_evidence_levels", "; ".join("%s %d" % kv for kv in evid.most_common()), "meta/labels_v2.tsv evidence_level (lab D)")

loci = rtsv(W("tree_v3", "loci_table.tsv"))
put("n_loci_table", "%d rows, kept=1: %d" % (len(loci), sum(r["kept"] == "1" for r in loci)), "work/tree_v3/loci_table.tsv")

# ---- ASR / transitions
tr = rtsv(W("final_v3", "asr", "transitions.tsv"))
sd = [int(r["S,D"]) for r in tr]; ds = [int(r["D,S"]) for r in tr]
q = lambda v, p: sorted(v)[int(p * (len(v) - 1))]
put("smap_S_to_D", "mean %.1f (95%% %d-%d), n maps %d" % (st.mean(sd), q(sd, .025), q(sd, .975), len(sd)),
    "work/final_v3/asr/transitions.tsv S,D")
put("smap_D_to_S", "mean %.1f (95%% %d-%d)" % (st.mean(ds), q(ds, .025), q(ds, .975)), "work/final_v3/asr/transitions.tsv D,S")
do = rtsv(W("final_v3", "asr", "deep_origins.tsv"))
nt = [int(r["n_tips"]) for r in do]
put("asr_deep_clades", "%d maximal deep clades (single-tip %d, largest %d tips)" % (len(do), sum(x == 1 for x in nt), max(nt)),
    "work/final_v3/asr/deep_origins.tsv")
fit = {r["model"]: r for r in rtsv(W("final_v3", "asr", "asr_fit.tsv"))}
put("mk_best", "ARD (dAIC ER %.2f; irreversible models dAIC > 126)" % g(fit["ER"]["dAIC"]), "work/final_v3/asr/asr_fit.tsv")

# ---- dating
rep = open(W("mcmctree_bb", "interp", "report.txt"), encoding="utf-8").read()
m = re.search(r"r\(clock-2 posterior mean, LSD2\) ([0-9.]+)", rep)
put("mcmctree_r_lsd2", m.group(1), "work/mcmctree_bb/interp/report.txt")
m = re.search(r"root: point ([0-9.]+) Ma, LSD2 ([0-9.]+) Ma, draws ([0-9.]+) \(([0-9.]+)-([0-9.]+)\)", rep)
put("root_ages", "MCMCTree-backbone point %s Ma (draws %s, %s-%s); LSD2 %s Ma" % (m.group(1), m.group(3), m.group(4), m.group(5), m.group(2)),
    "work/mcmctree_bb/interp/report.txt")
m = re.search(r"backbone nodes (\d+) of (\d+)", rep)
put("mcmctree_backbone_nodes", "%s of %s" % m.groups(), "work/mcmctree_bb/interp/report.txt")

# ---- topology tests
au = open(W("final_v3", "au.iqtree"), encoding="utf-8").read()
blk = au[au.index("USER TREES"):]
aus = re.findall(r"^\s+(\d)\s+(-[0-9.]+)\s+([0-9.]+).*?([0-9.e-]+) [+-] $", blk, re.M)
put("au_test", "; ".join("tree %s deltaL %s p-AU %s" % (a[0], a[2], a[3]) for a in aus) + " (1 fixed, 2 wASTRAL, 3 concatenated ML)",
    "work/final_v3/au.iqtree USER TREES")
gcf = [r for r in rtsv(W("final_v3", "gcf.cf.stat"))]
v = [g(r["gCF"]) for r in gcf if r["gCF"] not in ("", "NA")]
put("gcf_median", "%.1f%% (n branches %d)" % (st.median(v), len(v)), "work/final_v3/gcf.cf.stat")

# ---- 25b per clade real vs null
TESTS = ["RELAX p<0.05", "BUSTED-PH FG p<0.05", "BUSTED-PH strict (FG p<0.05 and diff p<0.05)",
         "CSUBST perm p<0.05", "CSUBST qualified pairs"]
r25 = []
for cl in CLADES + ["Actinopteri"]:
    for sub, tag in (("", "rescued"), ("norescue", "norescue")):
        f = W("codon", cl, sub, "25b_rescued.log") if sub == "" else W("codon", cl, "norescue", "25b.log")
        if not os.path.exists(f):
            alt = [x for x in os.listdir(W("codon", cl, sub)) if x.endswith(".log")] if os.path.isdir(W("codon", cl, sub)) else []
            if not alt:
                continue
            f = W("codon", cl, sub, alt[0])
        for r in rtsv(f):
            if r.get("test") in TESTS:
                r25.append({"clade": cl, "set": tag, **r})
wtsv("src_25b_perclade.tsv", r25, ["clade", "set", "test", "real_frame", "n_real", "real_mean", "n_null", "null_mean",
                                   "n_paired", "paired_T", "paired_perm_p"])

# ---- 25c combined + derived sensitivities
for s in ("final_main", "si_norescue", "si_fishold"):
    for r in rtsv(W("combine25c", s + ".combined.tsv")):
        put("25c_%s|%s" % (s, r["test"]), "z %s, p %s, Holm %s, clades p<0.05 %s/%s" % (
            r["stouffer_z"], r["p_stouffer"], r["p_holm"], r["n_clades_p<0.05"], r["k_clades"]), "work/combine25c/%s.combined.tsv" % s)
pc = rtsv(W("combine25c", "final_main.perclade.tsv"))
wtsv("src_25c_perclade.tsv", pc, list(pc[0].keys()))


def stouffer(ps, ws):
    zs = [norm.isf(min(max(p, 1e-4), 1 - 1e-4)) for p in ps]
    z = sum(w * x for w, x in zip(ws, zs)) / math.sqrt(sum(w * w for w in ws))
    return z, norm.sf(z)


der = []
for t in TESTS:
    rr = [r for r in pc if r["test"] == t and r["paired_perm_p"] not in ("", "NA", "nan")]
    ps = [g(r["paired_perm_p"]) for r in rr]; ws = [math.sqrt(g(r["n_paired"])) for r in rr]
    z, p = stouffer(ps, [1] * len(ps))
    der.append({"test": t, "variant": "unweighted", "left_out": "", "z": "%.3f" % z, "p": "%.3g" % p})
    for i, r in enumerate(rr):
        z, p = stouffer(ps[:i] + ps[i + 1:], ws[:i] + ws[i + 1:])
        der.append({"test": t, "variant": "LOCO", "left_out": r["clade"], "z": "%.3f" % z, "p": "%.3g" % p})
wtsv("src_25c_derived_sensitivity.tsv", der, ["test", "variant", "left_out", "z", "p"])
for t in TESTS:
    lo = [d for d in der if d["test"] == t and d["variant"] == "LOCO"]
    un = [d for d in der if d["test"] == t and d["variant"] == "unweighted"][0]
    worst = max(lo, key=lambda d: float(d["p"]))
    put("25c_derived|" + t, "unweighted z %s p %s; LOCO max p %s (without %s)" % (un["z"], un["p"], worst["p"], worst["left_out"]),
        "build/src_25c_derived_sensitivity.tsv (computed by harvest_sty.py from final_main.perclade.tsv)")

# ---- A4 confounds and direction of K
for r in rtsv(W("addon", "A4", "decision.tsv")):
    put("A4|" + r["ind"], "M1 OR %.2f, clustered p %.2g; C1 p %.2g; %s" % (g(r["M1_OR"]), g(r["M1_p_cl"]), g(r["C1_p"]), r["decision"]),
        "work/addon/A4/decision.tsv")
fpr = rtsv(W("addon", "A4", "null_fpr.tsv"))
seen = set(); kd = []
for r in fpr:
    key = (r["clade"], r["ind"])
    if key in seen:
        continue
    seen.add(key); kd.append(r)
wtsv("src_relax_direction.tsv", kd, ["clade", "ind", "n", "k", "real_rate", "rate", "lo", "hi"])
for r in kd:
    if r["clade"] == "all":
        put("relax_rates_all|" + r["ind"], "real %.3f vs matched null %.3f" % (g(r["real_rate"]), g(r["rate"])), "work/addon/A4/null_fpr.tsv")
pcl = rtsv(W("addon", "A4", "perclade.tsv"))
for ind in ("relax_relaxed", "relax_intens"):
    rr = [r for r in pcl if r["ind"] == ind and r["term"] == "is_real"]
    up = sum(g(r["est"]) > 0 and g(r["p_cl"]) < 0.05 for r in rr); dn = sum(g(r["est"]) < 0 and g(r["p_cl"]) < 0.05 for r in rr)
    put("A4_perclade|" + ind, "real>null p<0.05 in %d, real<null p<0.05 in %d of %d clades" % (up, dn, len(rr)), "work/addon/A4/perclade.tsv (M1 per clade)")

# ---- composition
ctl = {r["trait"]: r for r in rtsv(W("final_v3", "pgls_controls", "controls_summary.tsv"))}
P6 = ["busco_f_T", "whole_pI", "whole_R_RK", "whole_acid_base", "whole_net_charge7_per_res", "busco_f_F"]
comp = []
for t in P6:
    r = ctl[t]
    row = {"trait": t}
    for s in ("all", "deepT", "deepTGC", "nochemo", "deepdose", "ldepthT", "cs"):
        row[s + "_est"] = r.get(s + "_est", ""); row[s + "_p"] = r.get(s + "_p", "")
    comp.append(row)
a1 = {r["trait"]: r for r in rtsv(W("addon", "A1", "decision.tsv"))}
a1r = {r["trait"]: r for r in [dict(zip(["trait", "spec", "n", "shift", "orig", "corr", "boot", "h"], l.rstrip("\n").split("\t")))
                               for l in open(W("addon", "A1_reannot", "refit_report.txt"), encoding="utf-8") if "\tS_route3\t" in l]}
a2 = {(r["trait"], r["set"], r["adj"]): r for r in rtsv(W("addon", "A2", "summary.tsv"))}
gn = {r["trait"]: r for r in rtsv(W("addon", "A1", "gene_null_summary.tsv"))}
gmap = {"busco_f_T": "f_T", "whole_pI": "pI", "whole_R_RK": "R_RK", "whole_acid_base": "acid_base",
        "whole_net_charge7_per_res": "net_charge7_per_res", "busco_f_F": "f_F"}
for row in comp:
    t = row["trait"]
    row.update({"holm": a1[t]["holm"], "mp_est": a1[t]["mp_est"], "mp_p": a1[t]["mp_p"], "an_est": a1[t]["an_est"], "an_p": a1[t]["an_p"],
                "p_sim": a1[t]["p_sim"], "a1_call": a1[t]["call"] or "route- or quality-sensitive",
                "reannot_corr": a1r[t]["corr"], "reannot_boot": a1r[t]["boot"], "reannot_h": a1r[t]["h"],
                "origins_n": a2[(t, "disjoint", "raw")]["n"], "origins_median": a2[(t, "disjoint", "raw")]["median"],
                "origins_wilcox_holm": a2[(t, "disjoint", "raw")]["wilcox_holm"], "origins_npos": a2[(t, "disjoint", "raw")]["n_pos"],
                "origins_radj_holm": a2[(t, "disjoint", "route_adj")]["wilcox_holm"],
                "gene_frac_pos": gn[gmap[t]]["f_obs"], "gene_p_cal": gn[gmap[t]]["p_cal"]})
wtsv("src_composition.tsv", comp, list(comp[0].keys()))
for row in comp:
    put("comp|" + row["trait"], "PGLS %.3f SD (p %.2g, Holm %.2g, p_sim %.3g); deep+Thab %.3f (p %.2g); nochemo %.3f (p %.2g); re-annotation corrected %s; "
        "origins %s/%s positive, Wilcoxon Holm %.3g; gene-level frac %.3f p_cal %.3g" % (
            g(row["all_est"]), g(row["all_p"]), g(row["holm"]), g(row["p_sim"]), g(row["deepT_est"]), g(row["deepT_p"]), g(row["nochemo_est"]), g(row["nochemo_p"]),
            row["reannot_corr"], row["origins_npos"], row["origins_n"], g(row["origins_wilcox_holm"]), g(row["gene_frac_pos"]), g(row["gene_p_cal"])),
        "pgls_controls/controls_summary.tsv; addon A1, A1_reannot, A2")
t1 = open(W("addon", "T1", "report.txt"), encoding="utf-8").read()
m = re.search(r"busco_f_T\n\s+lab\s+deep ([0-9.]+).*?log10 depth \+ Thab ([0-9.]+) \(p ([0-9.]+)\).*?within deep ([0-9.-]+) \(p ([0-9.]+), n (\d+)", t1, re.S)
put("T1_fT", "log10 depth slope %s SD per log10 m (p %s); within deep %s (p %s, n %s)" % m.group(2, 3, 4, 5, 6), "work/addon/T1/report.txt")
a1rep = open(W("addon", "A1", "report.txt"), encoding="utf-8").read()
put("GC3_gene_level", re.search(r"gene-level GC3: (.*)", a1rep).group(1), "work/addon/A1/report.txt")
for top in ("final_v3_ml", "final_v3_unc"):
    r = {x["trait"]: x for x in rtsv(W(top, "pgls_controls", "controls_summary.tsv"))}["busco_f_T"]
    put("fT_" + top, "%.3f (p %.2g)" % (g(r["all_est"]), g(r["all_p"])), "work/%s/pgls_controls/controls_summary.tsv" % top)
for lb in ("pgls_controls_strict", "pgls_controls_1000"):
    r = {x["trait"]: x for x in rtsv(W("final_v3", lb, "controls_summary.tsv"))}["busco_f_T"]
    put("fT_" + lb, "%.3f (p %.2g, n %s)" % (g(r["all_est"]), g(r["all_p"]), r["all_n"]), "work/final_v3/%s/controls_summary.tsv" % lb)
a3 = rtsv(W("addon", "A3", "decision.tsv"))[0]
put("A3_prediction", "M_aa_res Stouffer p %.3f; route baseline >= traits in %s/%s clades; predictable %s" % (
    g(a3["stouffer_p_M_aa_res"]), a3["n_route_ge_M_aa_res"], a3["n_clades"], a3["predictable"]), "work/addon/A3/decision.tsv")

# ---- RER
rer = []
for lb in ("lab", "lab_nochemo", "lab_strict", "lab_1000"):
    for r in rtsv(W("rer_v3", "rer", lb, "genomewide.tsv")):
        rer.append({"labels": lb, "tree": "fixed", **r})
for top in ("rer_v3_unc", "rer_v3_ml"):
    for r in rtsv(W(top, "rer", "lab", "genomewide.tsv")):
        rer.append({"labels": "lab", "tree": top, **r})
wtsv("src_rer.tsv", rer, ["labels", "tree", "test", "statistic", "observed", "null_mean", "null_q025", "null_q975", "p_upper", "n_perm"])
for r in rer:
    if r["test"] == "binary_deep" and r["statistic"] == "n_pos_p05":
        put("rer|%s|%s" % (r["labels"], r["tree"]), "%s accelerated loci vs null %.1f, p %.3g" % (r["observed"], g(r["null_mean"]), g(r["p_upper"])),
            "work/%s/rer/%s/genomewide.tsv" % ("rer_v3" if r["tree"] == "fixed" else r["tree"], r["labels"]))

# ---- CAFE
caf = []
for cl in CLADES:
    f = W("cafe", cl, "cafe_summary.tsv")
    txt = open(f, encoding="utf-8").read()
    rows = rtsv(f)
    dp = [r for r in rows if r["model"] == "deep"][0]
    mm = re.search(r"#matched:.*?= ([0-9.]+)", txt); pm = re.search(r"#permutation:.*?= ([0-9.]+)", txt)
    caf.append({"clade": cl, "deep_gain": dp["gain_vs_base"], "n_sig_families_deep": dp["n_sig_families"],
                "rank_p_matched": mm.group(1) if mm else "NA", "rank_p_perm": pm.group(1) if pm else "NA"})
wtsv("src_cafe.tsv", caf, ["clade", "deep_gain", "n_sig_families_deep", "rank_p_matched", "rank_p_perm"])
put("cafe", "; ".join("%s %s/%s" % (r["clade"], r["rank_p_matched"], r["rank_p_perm"]) for r in caf) + " (rank p matched/permutation)", "work/cafe/*/cafe_summary.tsv")

# ---- function
rec = rtsv(W("func_v2", "combined", "recurrence_summary.tsv"))
wtsv("src_recurrence.tsv", rec, list(rec[0].keys()))
for r in rec:
    put("recurrence|%s|k%s" % (r["test"], r["k"]), "%s HOGs vs null %.1f (97.5%% %s), p %.3g" % (r["obs_OG"], g(r["null_mean"]), r["null_975"], g(r["perm_p"])),
        "work/func_v2/combined/recurrence_summary.tsv")
for t in ("relax", "bph_strict", "hyphy_any", "cs_qual"):
    rr = rtsv(W("addon", "D3_func_competitive", "competitive_%s.tsv" % t))
    qcol = [k for k in rr[0] if k.lower() in ("q", "q_bh", "qval")][0]
    put("D3|" + t, "%d terms tested, %d at q<0.05" % (len(rr), sum(g(r[qcol]) < 0.05 for r in rr)), "work/addon/D3_func_competitive/competitive_%s.tsv" % t)
a5 = rtsv(W("addon", "A5", "combined.tsv"))
calls = Counter(r["call"] for r in a5 if r["call"])
put("A5_calls", "; ".join("%s %d" % kv for kv in calls.most_common()), "work/addon/A5/combined.tsv")
a6 = rtsv(W("addon", "A6", "decision.tsv"))
put("A6_timing", "%d intervals, calls: %s" % (len(a6), "; ".join("%s %d" % kv for kv in Counter(r["call"] for r in a6).items())), "work/addon/A6/decision.tsv")

with open(os.path.join(B, "numbers_final.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["key", "value", "source"]); w.writerows(NUM)
print("%d numbers -> build/numbers_final.tsv" % len(NUM))
