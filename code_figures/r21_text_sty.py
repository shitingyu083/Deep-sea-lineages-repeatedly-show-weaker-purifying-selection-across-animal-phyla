"""r21 (2026-11-17): author block, end matter, AI section removed, last Results heading, abstract to ~150 words,
stacked caveats condensed in Results/Discussion. Numbers unchanged; every number removed from the main text
remains in the SI or Methods. Run once from paper/: python scripts/r21_text_sty.py submission_ready_20261117_r21"""
import sys, re, hashlib, shutil
P = sys.argv[1]; M = f"{P}/manuscript/main.md"
h = lambda s: hashlib.md5(s.encode()).hexdigest()[:8]
m = open(M, encoding="utf-8").read()
assert h(m) == "c3565319", h(m)
shutil.copy2(M, "backups/r21/main_r20b_c3565319.md")
R = []
# ---- author block --------------------------------------------------------------------------------------------
R.append(("Authors: [AUTHOR TO CONFIRM: names, affiliations, ORCID, corresponding author]\n",
 "Tingyu Shi^1,2^, Meijie Yang^1–3^, Zhen Zhao^1,2^, Zhuoqing Li^1,2^, Pu Shi^1,2^, Yu Wu^1,2^, Xiaoyan Liu^4^, Runqing Wang^1,2^, "
 "Jiameng Zhang^1,2^, Tao Zhang^5,*^ & Hao Song^1–3,*^\n\n"
 "1. CAS Key Laboratory of Marine Ecology and Environmental Sciences, Institute of Oceanology, Chinese Academy of Sciences, Qingdao 266071, China.\n"
 "2. University of Chinese Academy of Sciences, Beijing 100049, China.\n"
 "3. Laboratory for Marine Ecology and Environmental Science, Qingdao Marine Science and Technology Center, Qingdao 266237, China.\n"
 "4. Qingdao Agricultural University, Qingdao 266109, China.\n"
 "5. School of Marine Sciences, Ningbo University, Ningbo 315211, China.\n\n"
 "*Corresponding authors: Tao Zhang (zhangtaocas@nbu.edu.cn) and Hao Song (haosong@qdio.ac.cn).\n"))
# ---- abstract (183 -> ~150 words) ------------------------------------------------------------------------------
a0 = m[m.index("## Abstract\n\n") + 13:]; a0 = a0[:a0.index("\n")]
R.append((a0,
 "Animals have moved from shallow water into the deep sea many times, providing replicated transitions in which to ask "
 "whether evolution repeats through the same genes or through a shared change in the strength of selection. Here we "
 "reconstruct depth habitat for 426 species from 14 phyla and compare deep-sea lineages in eight clades with "
 "shallow-water lineages from the same phylogenies. Significant calls of relaxed selection were more than twice as "
 "frequent on deep-sea branches as on matched shallow branches and arose mainly from weaker purifying selection. The "
 "shift persisted after uniform genome re-annotation, recurred in most independent deep-sea lineages, was reproduced in "
 "a second sample of genes and had a proteome-wide component. By contrast, relaxed genes were shared between clades no "
 "more often than each clade's own rate predicts. Repeated deep-sea transitions were thus associated with weaker "
 "selective constraint without detectable gene reuse, a shift that searches for convergent adaptation must take into account."))
# ---- Results ---------------------------------------------------------------------------------------------------
R.append(("(ratio 0.95, 0.84–1.05; label-by-annotation interaction −0.08, −0.24 to 0.07; Extended Data Fig. 8a)",
          "(label-by-annotation interaction −0.08, −0.24 to 0.07; Extended Data Fig. 8a)"))
R.append(("(1.17–1.53; Extended Data Fig. 8b,c), and intensification was depleted (uniform −0.58, −0.90 to −0.28). Relaxed calls agreed between annotations (κ 0.84 and 0.81 for deep and null branches; Extended Data Fig. 8d). "
          "Fit failures (5.1–7.5% per arm) affected real and null foregrounds alike (62 versus 43 orthogroups failing in only one; exact McNemar P = 0.079; Extended Data Fig. 8e), and across a 6 × 6 grid of outcome assignments for the failed fits, the lowest lower 95% bound of the log odds ratio was +0.72 (Supplementary Note 8). "
          "Refits of 194 orthogroups from new random starting points reproduced the relaxed call in 758 of 776 fits (97.7%), and the log odds ratio changed from 1.135 to 1.215 (Methods).",
          "(1.17–1.53; Extended Data Fig. 8b,c), intensification was depleted (−0.58, −0.90 to −0.28), and relaxed calls agreed between annotations (κ 0.84 and 0.81 for deep and null branches; Extended Data Fig. 8d). "
          "Fit failures (5.1–7.5% per arm) affected real and null foregrounds alike (62 versus 43 orthogroups failing in only one; exact McNemar P = 0.079; Extended Data Fig. 8e; Supplementary Note 8), and refits from new starting points reproduced 97.7% of relaxed calls (Supplementary Methods)."))
R.append(("Substitution saturation did not generate the excess: excluding the most diverged quarter of orthogroups, measured by nucleotide branch length, increased it (log odds ratio 1.30, 1.16–1.43, against 1.10 in all orthogroups), adjusting for foreground and gene-wide divergence gave 1.27 (1.13–1.40), and the excess remained within the most diverged quarter (0.42, 0.19–0.65). Deep foreground branches were also shorter than matched null branches (median paired ratio 0.65; Supplementary Note 8). "
          "Null foregrounds matched to the deep foreground on total branch length gave a similar estimate: in 103 orthogroups the excess was 0.132 (95% CI 0.053–0.213), against 0.111 with the standard null on the same orthogroups. It was significant in cnidarians and annelids (P = 0.004 and 0.015), whereas the combination over all eight clades was not (P = 0.60), because three clades (cartilaginous fishes, echinoderms and sponges) had no discordant pairs (Supplementary Note 8).",
          "Substitution saturation did not generate the excess. Excluding the most diverged quarter of orthogroups increased it (log odds ratio 1.30, 1.16–1.43), the excess remained within that quarter (0.42, 0.19–0.65), adjusting for foreground and gene-wide divergence gave 1.27 (1.13–1.40), and deep foreground branches were shorter than matched null branches (median paired ratio 0.65). "
          "Null foregrounds matched on total branch length gave a similar excess (0.132, 95% CI 0.053–0.213, against 0.111 with the standard null in the same 103 orthogroups), significant in cnidarians and annelids (P = 0.004 and 0.015) but not in the combination over clades (P = 0.60), in which three clades had no discordant pairs (Supplementary Note 8)."))
R.append(("Adding real-minus-null differences in species-level dN/dS, genome GC, genome completeness and tree length attenuated the association by about 42% (odds ratio from 4.27 to 2.31, 95% CI 1.67–3.19, P = 4.0 × 10^−7^; 2,296 orthogroups). Because species-level dN/dS is estimated from the same kind of data as RELAX, this adjustment also absorbs part of the gene-level signal, and 42% is an upper bound on the proteome-wide share (Supplementary Note 4).",
          "Adding real-minus-null differences in species-level dN/dS, genome GC, genome completeness and tree length attenuated the association by about 42% (odds ratio from 4.27 to 2.31, 95% CI 1.67–3.19, P = 4.0 × 10^−7^; 2,296 orthogroups), an upper bound on the proteome-wide share because species-level dN/dS is estimated from the same sequence data (Supplementary Note 4)."))
R.append(("We repeated the analysis on a second sample of orthogroups, drawn per clade with a recorded seed and fitted on the production alignments and trees (Methods). This sample re-draws genes, not colonization events: it covered the same 65 lineages and 16 more, including four sponges.",
          "We repeated the analysis on a second sample of orthogroups (Methods), which re-draws genes rather than colonization events and covered the same 65 lineages and 16 more, including four sponges."))
R.append(("It stayed below zero without the two most sampled clades (−0.045, −0.064 to −0.025) and without each clade's largest lineage (−0.074, −0.095 to −0.050).",
          "It stayed below zero without the two most sampled clades (−0.045, −0.064 to −0.025) or each clade's largest lineage (−0.074, −0.095 to −0.050)."))
R.append(("Arthropods went the other way (clade mean in the second sample +0.038, 0.006 to 0.083). All deep arthropods in our sample, but no shallow ones, have non-chromosome, mostly transcriptome-based data, so this contrast coincides with a difference in data type. Across clades, lineage effects were not associated with assembly or annotation characteristics (within-clade permutation, all P ≥ 0.10), and the pooled effect was stronger among chromosome-level lineages of the second sample (−0.125, −0.198 to −0.073; 33 lineages in seven clades; Supplementary Note 8). Lineage effects became more negative with depth (β = −0.129, −0.191 to −0.008; Holm P = 0.029), but this trend rested on one hadal lineage (without it −0.014, P = 0.40), and neither stem age nor a codon-usage proxy for effective population size predicted the lineage effect.",
          "Arthropods went the other way (clade mean in the second sample +0.038, 0.006 to 0.083), a contrast that coincides with data type: all deep arthropods in our sample, but no shallow ones, have non-chromosome, mostly transcriptome-based data. Across clades, lineage effects were not associated with assembly or annotation characteristics (within-clade permutation, all P ≥ 0.10) and were stronger among chromosome-level lineages (−0.125, −0.198 to −0.073; 33 lineages in seven clades; Supplementary Note 8). They became more negative with depth (β = −0.129, −0.191 to −0.008; Holm P = 0.029), a trend that rested on one hadal lineage (without it −0.014, P = 0.40), whereas stem age and a codon-usage proxy for effective population size did not predict them."))
R.append(("### No robust excess of shared genes between clades\n", "### No excess reuse of relaxed orthologues across clades\n"))
R.append(("(P = 0.42; Fig. 4b). The overlap between clades is what their separate rates of relaxation produce. Only two",
          "(P = 0.42; Fig. 4b). Only two"))
R.append(("more than expected under the first null (29.1, P = 0.003) and under the clade-rate null (29.3, P = 0.021), but not after Holm correction across three secondary tests (P = 0.063), so this is a nominal signal rather than robust evidence of shared targets.",
          "more than expected under the first null (29.1, P = 0.003) and nominally under the clade-rate null (29.3, P = 0.021; Holm P = 0.063 across three secondary tests)."))
R.append((" These RELAX power estimates use a stricter threshold than the calls and are therefore conservative (Methods).", " RELAX power was computed at a stricter threshold than the calls (Methods)."))
# ---- Discussion ------------------------------------------------------------------------------------------------
R.append(("Small, sparse deep-sea populations would reduce the efficacy of selection across the genome^2,37^, consistent with the higher proteome-wide dN/dS of deep species and with the lineage-wide share of the gene-level excess, as in cavefishes, a subterranean catfish and island voles^22–24^. Alternatively, constraint may weaken on functions that matter less in cold, dark, food-poor water^1,2^, in parallel with the reduced morphological constraint of deep-sea anglerfishes^16^. Our data do not yet separate the two: relaxation was not concentrated in visual, circadian or pigmentation genes, the depth trend rested on one hadal lineage and a codon-usage proxy for effective population size did not predict lineage effects. Population-genomic estimates of effective size in paired deep and shallow lineages would provide a direct test.",
          "Small, sparse deep-sea populations would reduce the efficacy of selection across the genome^2,37^, consistent with the higher proteome-wide dN/dS of deep species and with patterns in cavefishes, a subterranean catfish and island voles^22–24^. Alternatively, constraint may weaken on functions that matter less in cold, dark, food-poor water^1,2^, in parallel with the reduced morphological constraint of deep-sea anglerfishes^16^. Our data do not yet separate the two, and population-genomic estimates of effective size in paired deep and shallow lineages would provide a direct test."))
R.append(("Habitat labels rest on occurrence records biased towards shallow sampling, which stricter and 1,000-m label sets were used to bracket; the main results held on trees alternative to the fixed topology; and protein-coding orthologues do not capture regulatory or structural change. Branch-length matching gave an excess similar to the standard null but, in the subset of orthogroups where it was possible, did not by itself establish a significant excess across all clades. Gene-family and colonization-timing analyses had little power (Extended Data Fig. 10). Because gene-level recurrence could be tested in at most four clades, denser sampling within phyla will be needed to resolve weak or partial reuse.",
          "Habitat labels rest on occurrence records biased towards shallow sampling, protein-coding orthologues do not capture regulatory or structural change, and gene-family and colonization-timing analyses had little power (Extended Data Fig. 10). Because gene-level recurrence could be tested in at most four clades, denser sampling within phyla will be needed to resolve weak or partial reuse."))
# ---- end matter ------------------------------------------------------------------------------------------------
R.append(("will be deposited in [AUTHOR TO CONFIRM: Zenodo/Dryad DOI] on acceptance.",
          "will be deposited in Zenodo on acceptance and are available to editors and reviewers on request."))
R.append(("All analysis scripts will be deposited at [AUTHOR TO CONFIRM: GitHub URL and Zenodo DOI].",
          "All analysis scripts and the code used to generate the figures are available at https://github.com/shitingyu083/Deep-sea-lineages-repeatedly-show-weaker-purifying-selection-across-animal-phyla and will be archived in Zenodo on acceptance."))
R.append((" (occurrence records; local copy dated 1 September 2026) [AUTHOR TO CONFIRM: export DOI]", " (occurrence records; accessed 1 September 2026)."))
R.append(("Correspondence and requests for materials should be addressed to [AUTHOR TO CONFIRM].",
          "Correspondence and requests for materials should be addressed to Tao Zhang or Hao Song."))
R.append(("[AUTHOR TO CONFIRM: funding, computing resources, data providers]",
          "This research was supported by the National Natural Science Foundation of China (grant nos. 42522609 and 42376088) and the earmarked fund for the China Agriculture Research System (CARS-49). "
          "H.S. was further supported by the Taishan Scholars Program (grant no. tsqn202312261), the Excellent Young Scientists Fund of Shandong Province (grant no. ZR2024YQ005), the Youth Innovation Promotion Association of the Chinese Academy of Sciences (grant no. 2023217), "
          "the Huaguoshan Outstanding Talents “Innovation and Entrepreneurship Program” of Lianyungang City (grant no. LYGSCRC2025007) and the Ningbo University School Talent Engineering Project (Science): Hard-shell Clam Breeding and Intelligent Seedling Breeding Technology and Equipment (grant no. ZX2026000697). "
          "M.Y. was supported by the Shandong Province Taishan Young Experts Program (grant no. tsqn202408277)."))
R.append(("[AUTHOR TO CONFIRM: CRediT statement]",
          "T.S., T.Z. and H.S. conceived and designed the study. T.S. assembled the data, performed the phylogenetic and selection analyses and wrote the first draft of the manuscript. "
          "M.Y., Z.Z., Z.L., P.S., Y.W., X.L., R.W. and J.Z. contributed to data curation, analysis and interpretation of the results. T.Z. and H.S. supervised the project, and T.Z., H.S. and M.Y. acquired funding. "
          "All authors revised the manuscript and approved the final version."))
R.append(("The authors declare no competing interests. [AUTHOR TO CONFIRM]", "The authors declare no competing interests."))
R.append(("## Use of AI tools\n[AUTHOR TO CONFIRM: describe the use of AI assistance in analysis scripting and drafting, as required by the journal policy]\n\n", ""))
for a, b in R:
    n = m.count(a); assert n == 1, (n, a[:80]); m = m.replace(a, b)
assert "AUTHOR TO CONFIRM" not in m and "VERIFY" not in m, re.findall(r".{40}(?:AUTHOR TO CONFIRM|VERIFY).{20}", m)
open(M, "w", encoding="utf-8", newline="\n").write(m)
ab = m[m.index("## Abstract\n\n") + 13:]; ab = ab[:ab.index("\n")]
print("main", h(m), "abstract words", len(ab.split()))
