"""Build Supplementary Tables 1-10 (one xlsx, one sheet per sub-table) and Supplementary Data 1.

Usage (from paper/): python scripts/build_si_sty.py submission_ready_20261028
All values are copied from analysis outputs (data_srv/work, data_srv/meta, si_pull); nothing is typed in by hand
except the sheet captions and the preregistration status rows in S7b, which restate Methods.
"""
import csv, gzip, os, re, sys
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

PKG = sys.argv[1] if len(sys.argv) > 1 else "submission_ready_20261029"
D = "data_srv"
PULL = "../si_pull"
OUT = os.path.join(PKG, "supplementary")
os.makedirs(OUT, exist_ok=True)
MAIN8 = ["Actinopteri_m", "Chondrichthyes", "Mollusca", "Annelida", "Arthropoda", "Echinodermata", "Cnidaria",
         "Porifera"]


def rtsv(p, comment="#"):
    with open(p, encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith(comment)]
    r = csv.reader(lines, delimiter="\t")
    rows = [[c.strip('"') for c in x] for x in r]
    return rows[0], rows[1:]


def num(v):
    if not isinstance(v, str) or v in ("", "NA", "nan", "NaN"):
        return v
    try:
        x = float(v)
        return int(x) if re.fullmatch(r"-?\d+", v) else x
    except ValueError:
        return v


wb = Workbook()
toc = wb.active
toc.title = "Contents"
toc.append(["Supplementary Tables for: Recurrent relaxed purifying selection in deep-sea animal lineages"])
toc.append([])
toc.append(["Sheet", "Content", "Source file (analysis repository)"])
for c in toc[3]:
    c.font = Font(bold=True)


def sheet(name, caption, header, rows, source):
    ws = wb.create_sheet(name)
    ws.append([caption])
    ws["A1"].font = Font(bold=True)
    ws.append([f"Source: {source}"])
    ws.append(header)
    for c in ws[3]:
        c.font = Font(bold=True)
    for r in rows:
        ws.append([num(v) for v in r])
    for i, h in enumerate(header, 1):
        ws.column_dimensions[get_column_letter(i)].width = min(40, max(10, len(str(h)) + 2))
    ws.freeze_panes = "A4"
    toc.append([name, caption, source])
    return ws


# ---- S1 species, labels, accessions, QC ----------------------------------
tree = open(f"{D}/work/final_v3/bl.treefile").read()
tips = set(re.findall(r"[(,]([A-Za-z][A-Za-z0-9_]*)", tree))
hl, lab = rtsv(f"{D}/meta/labels_v2.tsv")
hc, cod = rtsv(f"{D}/meta/species_codes.tsv")
hq, qc = rtsv(f"{D}/meta/qc_table.tsv")
acc = {r[0]: r[2] for r in cod}
qcd = {r[0]: r for r in qc}
keep_l = ["code", "species", "phylum", "class", "order", "family", "clade", "role", "source", "lab", "lab_strict",
          "lab_1000", "evidence_level", "depth_median", "frac_gt200", "bg_reason"]
il = [hl.index(k) for k in keep_l]
keep_q = ["n_genes", "genome_C", "prot_C", "prot_S", "prot_D", "prot_F", "prot_M"]
iq = [hq.index(k) for k in keep_q]
rows = []
for r in lab:
    if r[0] in tips:
        q = qcd.get(r[0])
        rows.append([r[i] for i in il] + [acc.get(r[0], "")] + ([q[i] for i in iq] if q else [""] * len(iq)))
rows.sort(key=lambda x: (x[2], x[6], x[1]))
sheet("S1a_species", f"Supplementary Table 1a | The {len(rows)} taxa in the final tree: taxonomy, analysis clade, "
      "habitat labels (lab = main; lab_strict; lab_1000 = deep only if > 1,000 m), depth evidence, protein-set "
      "source, accession and BUSCO (metazoa_odb12) completeness (genome_C, prot_C/S/D/F/M, %).",
      keep_l + ["accession"] + keep_q, rows, "meta/labels_v2.tsv; meta/species_codes.tsv; meta/qc_table.tsv; "
      "work/final_v3/bl.treefile")
sheet("S1b_qc_all", "Supplementary Table 1b | Quality control of all candidate protein sets (status pass/fail, flags).",
      hq, qc, "meta/qc_table.tsv")
he, ex = rtsv(f"{D}/meta/excluded_assemblies.tsv")
sheet("S1c_excluded", "Supplementary Table 1c | Assemblies excluded before quality control, with reason.", he, ex,
      "meta/excluded_assemblies.tsv")
hd, dep = rtsv(f"{D}/meta/depth_basis_v3.tsv")
sheet("S1d_depth_basis", "Supplementary Table 1d | Depth evidence for taxa added in the v3 panel.", hd, dep,
      "meta/depth_basis_v3.tsv")

# ---- S2 loci -------------------------------------------------------------
ht, lt = rtsv(f"{D}/work/tree_v3/loci_table.tsv")
lt = {r[0]: r for r in lt}
rows = []
for l in open(f"{PULL}/work/tree_v3/concat_partition.txt", encoding="utf-8"):
    m = re.match(r"\s*([^,]+),\s*(\S+)\s*=\s*(\d+)-(\d+)", l)
    if m:
        model, loc, a, b = m.groups()
        t = lt.get(loc, [loc, "", "", ""])
        rows.append([loc, model, int(a), int(b), int(b) - int(a) + 1, t[1], t[2], t[3]])
sheet("S2_loci", f"Supplementary Table 2 | The {len(rows)} single-copy BUSCO loci (metazoa_odb12) in the concatenated "
      "alignment: per-locus model, alignment coordinates, trimmed length (sites) and occupancy.",
      ["busco_id", "model", "start", "end", "sites", "n_ingroup", "n_outgroup", "mean_len_aa_untrimmed"], rows,
      "work/tree_v3/concat_partition.txt; work/tree_v3/loci_table.tsv")

# ---- S3 calibrations -----------------------------------------------------
h, r = rtsv(f"{D}/work/final_v3/date/calibrations_used.tsv")
sheet("S3_calibrations", f"Supplementary Table 3 | The {len(r)} node calibrations (minimum and maximum ages, Ma) used "
      "for LSD2 and MCMCTree, with the taxa defining each node and the source node in ref. 3.", h, r,
      "work/final_v3/date/calibrations_used.tsv")

# ---- S4 proteome traits --------------------------------------------------
for k, d, lbl in (("a", "pgls_controls", "main labels"), ("b", "pgls_controls_strict", "strict labels"),
                  ("c", "pgls_controls_1000", "labels > 1,000 m")):
    h, r = rtsv(f"{D}/work/final_v3/{d}/controls_summary.tsv")
    sheet(f"S4{k}_pgls_{d.split('_')[-1]}", f"Supplementary Table 4{k} | Species-level PGLS deep effects for proteome "
          f"traits ({lbl}); each specification gives estimate, s.e., P, BH q and n (all = all taxa; ref = "
          "reference subset; deepT = + habitat temperature; GC = + genome GC; other columns as named).", h, r,
          f"work/final_v3/{d}/controls_summary.tsv")
h, r = rtsv(f"{D}/work/addon/A2/summary.tsv")
sheet("S4d_origin_contrasts", "Supplementary Table 4d | Deep-minus-shallow contrasts at the 65 disjoint deep origins "
      "(A2): median, weighted mean with 95% interval, sign and Wilcoxon tests, Holm-adjusted P.", h, r,
      "work/addon/A2/summary.tsv")
h, r = rtsv(f"{D}/work/addon/A1_reannot/refit.tsv")
sheet("S4e_reannotation", "Supplementary Table 4e | PGLS refit after re-annotating 20 deep and 20 shallow species "
      "through a common route (A1).", h, r, "work/addon/A1_reannot/refit.tsv")
h, r = rtsv(f"{D}/work/final_v3/pgls/pgls_gene_summary.tsv")
sheet("S4f_gene_level", "Supplementary Table 4f | Gene-level PGLS: fraction of loci with a positive deep effect, "
      "sign test and number of loci at q < 0.05.", h, r, "work/final_v3/pgls/pgls_gene_summary.tsv")
h, r = rtsv(f"{D}/work/addon/R1/threonine_gc_R1.tsv")
sheet("S4g_threonine_gc", "Supplementary Table 4g | Threonine control (PREREG_nature R1): origin contrasts in "
      "threonine frequency regressed on GC3 contrasts (gc3) or GC3 plus genome GC (gc3_gc), for all, high-quality "
      "and strict origin sets; intercept and P, GC3 slope, and the adjusted median with Wilcoxon P.", h, r,
      "work/addon/R1/threonine_gc_R1.tsv")

# ---- S5 candidate sets and enrichment ------------------------------------
h, r = rtsv(f"{D}/work/addon/A5/combined.tsv")
sheet("S5a_candidate_sets", "Supplementary Table 5a | Preregistered candidate sets (P1-P9) and the ribosomal negative "
      "control (N1): combined Z across clades, permutation P and Holm-adjusted P per selection statistic.", h, r,
      "work/addon/A5/combined.tsv")
h, r = rtsv(f"{D}/work/addon/A5/perclade.tsv")
sheet("S5b_sets_perclade", "Supplementary Table 5b | Candidate sets per clade (status lt3 = fewer than three HOGs, "
      "not tested).", h, r, "work/addon/A5/perclade.tsv")
rows = []
for t in ("relax", "bph_strict", "hyphy_any", "cs_qual"):
    h, r = rtsv(f"{PULL}/work/addon/D3_func_competitive/competitive_{t}.tsv")
    rows += [[t] + x for x in r]
sheet("S5c_GO_competitive", "Supplementary Table 5c | Post hoc competitive GO enrichment (D3) of per-HOG real-minus-"
      "null selection differences, by statistic; Z, P and BH q across terms, clade z and HOG counts.",
      ["statistic"] + h, rows, "work/addon/D3_func_competitive/competitive_*.tsv")
h, r = rtsv(f"{D}/work/func_v2/combined/recurrence_summary.tsv")
sheet("S5d_recurrence", "Supplementary Table 5d | Recurrence of selected orthogroups across clades against the "
      "label-permutation null.", h, r, "work/func_v2/combined/recurrence_summary.tsv")

# ---- S6 CAFE -------------------------------------------------------------
rows = []
for c in MAIN8:
    for sub in (c, c + "_rescue"):
        p = f"{D}/work/cafe/{sub}/cafe_summary.tsv"
        if os.path.exists(p):
            h, r = rtsv(p)
            rows += [[sub] + x for x in r]
h, r = rtsv(f"{D}/work/addon/R2/cafe_perm/cafe_rank40.tsv")
sheet("S6a_cafe_rank40", "Supplementary Table 6a | CAFE5 deep-rate gain ranked against 40 tip-permuted foregrounds "
      "per clade (PREREG_nature R2): gain = base minus two-rate negative log-likelihood; K = converged permutations; "
      "rank P = (1 + n_ge)/(1 + K); _ok columns exclude permutation fits with negative gain; Actinopteri_m_rescue "
      "deep fit failed (negative gain).", h, r, "work/addon/R2/cafe_perm/cafe_rank40.tsv")
h, r = rtsv(f"{D}/work/addon/R2/cafe_perm/cafe_gains_long.tsv")
sheet("S6b_cafe_gains", "Supplementary Table 6b | Every CAFE5 gain used in Extended Data Fig. 7c (real deep, matched shallow and "
      "permuted foregrounds).", h, r, "work/addon/R2/cafe_perm/cafe_gains_long.tsv")
sheet("S6c_cafe", "Supplementary Table 6c | CAFE5 gene-family fits per clade: base (global lambda + error model), deep "
      "(separate lambda on deep branches) and label-permuted null fits; negative log-likelihood, lambdas, error, gain "
      "over base, LRT P and number of significant families. '_rescue' rows: refits after removing high-range "
      "families (Methods).", ["fit_set"] + h, rows, "work/cafe/<clade>[_rescue]/cafe_summary.tsv")

# ---- S7 topology constraints, preregistration status, deviations ---------
h, r = rtsv(f"{D}/meta/topo_fix_v3.tsv")
sheet("S7a_topology", "Supplementary Table 7a | Literature topological constraints applied to the final tree.", h, r,
      "meta/topo_fix_v3.tsv")
status = [
    ["PREREG_addon A1-A6", "Run as frozen", "work/addon/FROZEN_addon.txt"],
    ["PREREG_addon T2 life-history covariates", "Not run: FishBase/SeaLifeBase not retrievable programmatically on the "
     "analysis server", "Methods"],
    ["combine25c PREREG sensitivity (a) original fish design", "Run (Extended Data Fig. 5)",
     "work/combine25c/si_fishold.*"],
    ["combine25c PREREG sensitivity (b) random deep m-subsets for fish", "Not run", "Methods"],
    ["Actinopteri_m v1 matched design", "Stopped during labelling, before any test; replaced by v2",
     "work/codon/Actinopteri_m/DEVIATIONS.md"],
    ["PREREG_fig4 full-tree MCMCTree", "Stopped (projected 49 days per chain); replaced by preregistered backbone",
     "work/func/DEVIATIONS_fig4.md"],
    ["D3 competitive GO test", "Post hoc, not preregistered; labelled as such", "work/addon/DEVIATIONS_addon.md"],
    ["CAFE5 rescue (Mollusca, Actinopteri_m)", "Post hoc family filtering after base-fit failure",
     "work/addon/DEVIATIONS_addon.md"],
    ["PREREG_nature R1 confounder models M1/M3/M5, PGLS of omega_tip, calibration, dose", "Run as frozen; outcome "
     "scenario B (attenuated) by the prestated rule", "work/addon/R1/SCENARIO.txt"],
    ["PREREG_nature R1 high-quality subset", "Prestated subset empty (0 orthogroups); relaxed definitions reported "
     "as exploratory (Table 8c)", "DEVIATIONS_addon.md D-R1-1"],
    ["PREREG_nature R1 threonine GC3 control", "Run as frozen; claim stands (intercept > 0, P < 0.05 in all sets)",
     "work/addon/R1/threonine_gc_R1.tsv"],
    ["PREREG_nature R1 depth thresholds 500/1,000 m", "Run after other R1 results were seen; descriptive; 1,000 m "
     "uninformative (13 orthogroups)", "DEVIATIONS_addon.md D-R-4"],
    ["PREREG_nature R2 power simulation", "Run; pilot gate passed after per-orthogroup tree-length calibration",
     "DEVIATIONS_addon.md D-R2-1, D-R2-2"],
    ["PREREG_nature R2 GO spike-in", "Run after the power results were seen", "DEVIATIONS_addon.md D-R-4"],
    ["PREREG_nature R2 CAFE 40 permutations", "Run; Arthropoda K = 21 converged", "DEVIATIONS_addon.md D-R2-3"],
    ["PREREG_nature R3 RTF1 selection tests", "Run; Echinodermata null fit failed (NA)", "DEVIATIONS_addon.md D-R3-1"],
    ["PREREG_nature R3 RTF1 fish amino-acid scan", "Run with 2 (not 4) fishes below 3,000 m; null sets of size 2",
     "DEVIATIONS_addon.md D-R-4"],
    ["PREREG_nature R3 genome-scale fish RER and GO rank test", "Run as frozen", "work/addon/R3/fishrer"],
    ["Permulation-calibrated fish GO test", "Post hoc, not preregistered", "DEVIATIONS_addon.md D-R-4"],
    ["PREREG_r5 D1 second null design (null2) and null-versus-null calibration", "Run after all other results were "
     "known; null-design robust (relax, intens, any), calibrated (Tables 8l-n); 123 of 2,348 null2 fits missing "
     "after one retry", "DEVIATIONS_addon.md D-R5-0, D-R5-2"],
    ["PREREG_r5 D2a weak-effect power (K 0.85, 0.90, K75h)", "Run; 18 of 480 replicates missing after one retry "
     "(Tables 10f,g)", "DEVIATIONS_addon.md D-R5-1, D-R5-2"],
    ["PREREG_r5 D2b clade-subset GO spike-in", "Run (Table 10h)", "DEVIATIONS_addon.md D-R5-0"],
    ["PREREG_r5 D3 route-pure subsets", "Run (pruning replaces the empty unified-annotation subset); verdict "
     "single_route_M (Tables 8o-r)", "DEVIATIONS_addon.md D-R5-0, D-R5-2"],
]
sheet("S7b_prereg_status", "Supplementary Table 7b | Status of every preregistered analysis.",
      ["analysis", "status", "record"], status, "PREREG_addon.md; combine25c/PREREG.md; PREREG_fig4.md; "
    "PREREG_nature.md; PREREG_r5.md")
rows = []
for f in ("work_addon_DEVIATIONS_addon.md", "work_codon_Actinopteri_m_DEVIATIONS.md", "work_func_DEVIATIONS_fig4.md"):
    sec = ""
    for l in open(f"{D}/dev/{f}", encoding="utf-8"):
        l = l.rstrip()
        if not l:
            continue
        if l.startswith("#"):
            sec = l.lstrip("# ")
            continue
        rows.append([f.replace("work_", "").replace("_", "/", 1), sec, l.lstrip("- ")])
sheet("S7c_deviations", "Supplementary Table 7c | Verbatim deviation records, with dates, from the analysis repository.",
      ["record", "entry", "text"], rows, "DEVIATIONS files")
wb["S7c_deviations"].column_dimensions["C"].width = 160
for c in wb["S7c_deviations"]["C"]:
    c.alignment = Alignment(wrap_text=True)

# ---- S7d preregistration timeline (2026-10-30 revision) -------------------
# Freeze times: timestamp line of each FROZEN record, or file mtime on the server (probe p661, 2026-10-30).
# Seen / not computed: the timing statements written in each PREREG document.
TL = [
    ["fish matched design (Actinopteri_m)", "work/codon/Actinopteri_m/FROZEN.txt", "2026-10-23 09:56:00 (v1 09:38:33)",
     "not stated in the freeze record", "fish matched-design selection results", "FROZEN record"],
    ["cross-clade combination (25c)", "work/combine25c/PREREG.md, FROZEN.txt", "2026-10-24 16:19:05",
     "25b per-clade results for 6 of 8 clades (rule proposed 2026-10-22 after 4)",
     "Mollusca and fish matched-design results; rescued RELAX merge; any combined statistic", "PREREG timing statement"],
    ["functional, gene-family and dating (Fig. 4, CAFE5, MCMCTree)", "work/func/PREREG_fig4.md, FROZEN_fig4.txt",
     "2026-10-25 02:47:13", "exploratory DEV runs on no-rescue 25b tables",
     "confirmatory run on final 25b tables", "PREREG status statement"],
    ["add-on analyses A1-A6", "work/addon/PREREG_addon.md, FROZEN_addon.txt", "2026-10-26 15:10:28",
     "PGLS, RERconverge, per-clade selection summaries; final 25c test-level tables and last lines of 48.log; "
     "focal traits P6 selected on the same data", "A4 indicators, A5 set statistics, CAFE5 results", "PREREG timing statement"],
    ["priority analyses R1-R3 (PREREG_nature)", "work/addon/PREREG_nature.md, .md5", "2026-10-28 15:49 (file mtime)",
     "all results reported before R1-R3", "R1-R3 outcomes; depth-threshold, GO spike-in and RTF1 scan run after "
     "confounder and power results", "PREREG header; Methods"],
    ["direction-specific RELAX tests and orthologue coverage (post hoc)", "scripts/67_relax_direction_sty.py, "
     "68_og_coverage_sty.py", "2026-10-30 15:21 / 15:27 (file mtime; not preregistered)",
     "all results including the combined RELAX statistic", "-", "deviation record 2026-10-30"],
    ["supplementary computations D1-D3 (PREREG_r5)", "work/addon/PREREG_r5.md (MD5 8209cdef), FROZEN_r5.txt",
     "2026-10-30 16:58", "all results above, including the direction-specific tests",
     "null2 RELAX fits, route-pure fits, weak-effect simulations, clade-subset spike-in", "PREREG_r5 header"],
]
ws7d = sheet("S7d_prereg_timeline", "Supplementary Table 7d | Preregistration timeline. MD5 checksums fix file "
             "content, not freeze time; times are from freeze records or server file modification times and are not "
             "externally certified.", ["stage", "documents", "frozen (CST)", "already seen at freeze",
                                       "not yet computed at freeze", "basis"], TL, "PREREG / FROZEN records")
for col in "BDE":
    ws7d.column_dimensions[col].width = 60

# ---- S8 lineage-wide confounders (PREREG_nature R1) -----------------------
A = f"{D}/work/addon"
for name, cap, p in (
        ("S8a_R1_models", "Supplementary Table 8a | Real-versus-null models of selection calls with lineage-wide "
         "covariates: M1 (baseline), M3 (+ d_omega), M5 (+ d_omega, d_gc, d_buscoC, d_tbl), M5b (GC3 for GC) and "
         "M6 (+ habitat temperature, exploratory); HOG-clustered (CR1) s.e. and P.", "R1/models_R1.tsv"),
        ("S8b_depth_threshold", "Supplementary Table 8b | Depth-threshold sensitivity: M0/M1 refitted with deep "
         "defined at 200, 500 or 1,000 m (descriptive).", "R1/depth_thresh/models_depth_thresh.tsv"),
        ("S8c_hq_explore", "Supplementary Table 8c | Exploratory relaxed high-quality subsets (the prestated subset "
         "was empty): models M0/M1/M5.", "R1/explore_hq/models_hq_explore.tsv"),
        ("S8d_hq_sizes", "Supplementary Table 8d | Sizes of the exploratory high-quality subsets.",
         "R1/explore_hq/hq_sizes.tsv"),
        ("S8e_pgls_meta", "Supplementary Table 8e | Meta-analysed PGLS deep effects on genome-wide terminal dN/dS "
         "(log omega_tip), GC3 and BUSCO completeness.", "R1/pgls_meta.tsv"),
        ("S8f_pgls_clade", "Supplementary Table 8f | Per-clade PGLS deep effects underlying Table 8e.",
         "R1/pgls_tip.tsv"),
        ("S8g_calibration", "Supplementary Table 8g | Calibration: real-foreground calls observed versus expected "
         "from a covariate model fitted on null rows only.", "R1/calibration_R1.tsv"),
        ("S8h_dose", "Supplementary Table 8h | Dose: selection calls on real foregrounds against log mean "
         "foreground depth, with M1 covariates.", "R1/dose_R1.tsv"),
        ("S8i_direction_clade", "Supplementary Table 8i | Direction-specific RELAX tests (post hoc, 2026-10-30; "
         "script 67): per-clade rates of any significant change (any), relaxation (P < 0.05, K < 1) and "
         "intensification (P < 0.05, K > 1) on real and matched null foregrounds, with paired-permutation upper- "
         "and lower-tail P (N = 10,000, seed 20261022). The combined tests are in Table 8j.",
         "R4_direction/perclade.tsv"),
        ("S8j_direction_combined", "Supplementary Table 8j | Weighted Stouffer combination of Table 8i (relaxation: "
         "upper tail; intensification: lower tail; any: upper tail, reproducing the preregistered statistic) with "
         "leave-one-clade-out range; Holm correction over the two post hoc tests is applied in the text.",
         "R4_direction/combined.tsv"),
        ("S8k_og_coverage", "Supplementary Table 8k | Coverage of the recurrence test (post hoc, 2026-10-30; script "
         "68): metazoan eggNOG orthologous groups by the number of clades in which they had a paired real and null "
         "result (n_og) and cumulative count tested in at least that many clades (n_og_ge).",
         "R4_direction/og_coverage_summary.tsv"),
        ("S8l_null2_clade", "Supplementary Table 8l | Alternative null design (PREREG_r5 D1; designed after all "
         "primary results were known): per-clade RELAX rates for any change, relaxation and intensification on "
         "orthogroups with real, null1 and null2 results; C1 real vs null2, C2 real vs null1, C3 null2 vs null1 "
         "(paired sign-flip, N = 10,000, seed 20261022).", "R5/D1/perclade.tsv"),
        ("S8m_null2_combined", "Supplementary Table 8m | Weighted Stouffer combination of Table 8l with "
         "leave-one-clade-out range and pooled rate difference T.", "R5/D1/combined.tsv"),
        ("S8n_null2_verdict", "Supplementary Table 8n | Prestated D1 verdicts (null-design robustness, "
         "null-versus-null calibration).", "R5/D1/verdict.tsv"),
        ("S8o_routepure_models", "Supplementary Table 8o | Route-pure subsets (PREREG_r5 D3): HOG-clustered "
         "binomial models of RELAX calls on trees pruned to one annotation route (A, NCBI/author annotation; M, "
         "miniprot) and on the unpruned trees of the same orthogroups (R0, real + clade; R1, + log tips + gap "
         "difference).", "R5/D3/models.tsv"),
        ("S8p_routepure_perm", "Supplementary Table 8p | Per-clade paired permutation tests on the route-pure "
         "subsets of Table 8o.", "R5/D3/perm_perclade.tsv"),
        ("S8q_routepure_combined", "Supplementary Table 8q | Weighted Stouffer combination of Table 8p.",
         "R5/D3/perm_combined.tsv"),
        ("S8r_routepure_verdict", "Supplementary Table 8r | Prestated D3 verdicts per route and overall.",
         "R5/D3/verdict.tsv"),
        ("S8s_fit_failures", "Supplementary Table 8s | RELAX fits and simulation replicates without a result after "
         "one retry (HyPhy ancestral-reconstruction error; deviation D-R5-2), by analysis, clade, route, frame and "
         "set; missing fits were treated as missing, not imputed.", "R5/fail_counts.tsv")):
    h, r = rtsv(f"{A}/{p}")
    sheet(name, cap, h, r, f"work/addon/{p}")

# ---- S9 reconciliation with published candidates (PREREG_nature R3) ------
for name, cap, p in (
        ("S9a_fishRER_genomewide", "Supplementary Table 9a | Genome-scale RERconverge in ray-finned fishes (75 "
         "species, 46 deep; 10,086 gene trees): accelerated and decelerated loci (binary label) and log-depth test "
         "against 1,000 permulations.", "R3/fishrer/lab/genomewide.tsv"),
        ("S9b_fishRER_GO", "Supplementary Table 9b | Gene Ontology rank test of fish relative rates: AUC, parametric "
         "Mann-Whitney P and BH q, and permulation-calibrated P and BH q (1,000 permulations; post hoc).",
         "R3/fishrer/go/fishrer_go_perm.tsv"),
        ("S9c_RTF1_aascan", "Supplementary Table 9c | RTF1 amino-acid scan in fishes: sites at which the fishes "
         "living below 3,000 m share a residue absent from all shallow fishes, against 1,000 random shallow sets of "
         "the same size.", "R3/rtf1/aascan.tsv"),
        ("S9d_RTF1_sites", "Supplementary Table 9d | The observed RTF1 sites of Table 9c (alignment coordinates).",
         "R3/rtf1/aascan_sites.tsv"),
        ("S9e_RTF1_selection", "Supplementary Table 9e | RELAX, BUSTED-PH and CSUBST on the RTF1 orthogroup per "
         "clade, real and matched-null (null1) labels; K = RELAX intensity parameter.", "R3/rtf1/rtf1_hyphy.tsv")):
    h, r = rtsv(f"{A}/{p}")
    sheet(name, cap, h, r, f"work/addon/{p}")

# ---- S10 power (PREREG_nature R2) ----------------------------------------
for name, cap, p in (
        ("S10a_power_clade", "Supplementary Table 10a | Simulated power and false-positive rate per clade and "
         "scenario (H0r/H0b no change; K75, K50 relaxation; K50h half of foreground lineages; P1, P2 episodic "
         "positive selection), at the nominal P equivalent to q < 0.05 in the real analysis.",
         "R2/main/power_by_clade.tsv"),
        ("S10b_informative", "Supplementary Table 10b | Prestated informativeness (power >= 0.8) per clade and "
         "scenario.", "R2/main/informative.tsv"),
        ("S10c_recurrence_power", "Supplementary Table 10c | Probability that an orthogroup affected in all covered "
         "clades is significant in at least k clades.", "R2/main/recurrence_power.tsv"),
        ("S10d_GO_spikein", "Supplementary Table 10d | Power of the competitive GO test under spike-in, by "
         "scenario, affected fraction and term size.", "R2/go_spikein/spikein_power.tsv"),
        ("S10e_GO_spikein_terms", "Supplementary Table 10e | Per-term spike-in results underlying Table 10d.",
         "R2/go_spikein/spikein_terms.tsv"),
        ("S10f_power_weak", "Supplementary Table 10f | Power for weaker relaxation (PREREG_r5 D2a; K85, K90: K = "
         "0.85, 0.90 on all foreground branches; K75h: K = 0.75 on a random half; rows with source R2 repeat Table "
         "10a for comparison), at the nominal P equivalent to q < 0.05.", "R5/D2/power_by_clade.tsv"),
        ("S10g_recurrence_weak", "Supplementary Table 10g | Probability of significant calls in at least three of k "
         "covered clades under the scenarios of Table 10f (mean, minimum and maximum over clade combinations).",
         "R5/D2/recurrence_power.tsv"),
        ("S10h_GO_spikein_subset", "Supplementary Table 10h | Competitive GO test power (q < 0.05) under spike-in by "
         "scenario, affected fraction of a term's orthogroups and number m of the term's tested clades carrying the "
         "effect (PREREG_r5 D2b, seed 20261030).", "R5/go_spikein/spikein_power.tsv"),
        ("S10i_GO_spikein_subset_terms", "Supplementary Table 10i | Per-term results underlying Table 10h.",
         "R5/go_spikein/spikein_terms.tsv")):
    h, r = rtsv(f"{A}/{p}")
    sheet(name, cap, h, r, f"work/addon/{p}")

toc.column_dimensions["A"].width = 24
toc.column_dimensions["B"].width = 120
toc.column_dimensions["C"].width = 60
wb.save(os.path.join(OUT, "Supplementary_Tables_1-10.xlsx"))

# ---- Supplementary Data 1: per-gene selection results (real rows) --------
n = 0
with gzip.open(os.path.join(OUT, "Supplementary_Data_1_per_gene_selection.tsv.gz"), "wt", encoding="utf-8") as g:
    hdr = None
    for c in MAIN8:
        h, r = rtsv(f"{PULL}/work/codon/{c}/sel_results_v3.tsv")
        if hdr is None:
            hdr = h
            g.write("clade\t" + "\t".join(h) + "\n")
        assert h == hdr, c
        for x in r:
            if x[1] == "real":
                g.write(c + "\t" + "\t".join(x) + "\n")
                n += 1
print("sheets", len(wb.sheetnames) - 1, "per-gene rows", n)
