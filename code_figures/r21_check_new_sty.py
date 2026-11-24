"""r21: values dropped from the main text in the caveat condensation must still be in the SI; new r21 wording present;
no placeholders, no AI-use section. Retargets r9/r10/r11/r18 expectations whose main-text wording changed in r21."""
import sys
P = sys.argv[1]; M = open(f"{P}/manuscript/main.md", encoding="utf-8").read(); S = open(f"{P}/supplementary/Supplementary_Information.md", encoding="utf-8").read()
C = [("ratio 0.95 in SI", "0.95, 0.84", S), ("grid +0.72 in SI", "+0.72", S), ("refit 758/776 in SI", "758 of 776", S), ("refit logOR in SI", "1.135", S),
     ("refit 97.7 main", "reproduced 97.7% of relaxed calls", M), ("62/43 main", "62 versus 43 orthogroups", M),
     ("sat all 1.10 main", "1.10 (0.98–1.23)", M), ("intens uniform", "(−0.58, −0.90 to −0.28)", M),
     ("bmatch E_B", "(0.132, 95% CI 0.053–0.213, against 0.111", M), ("bmatch combined", "combination over clades (P = 0.60)", M),
     ("bph holm", "Holm P = 0.063 across three secondary tests", M), ("heading", "### No excess reuse of relaxed orthologues across clades", M),
     ("42 upper bound", "an upper bound on the proteome-wide share because species-level dN/dS is estimated from the same sequence data", M),
     ("authors", "Tao Zhang^5,*^ & Hao Song^1–3,*^", M), ("github", "https://github.com/shitingyu083/Deep-sea-lineages-repeatedly-show-weaker-purifying-selection-across-animal-phyla", M),
     ("funding", "42522609 and 42376088", M)]
bad = 0; out = ["| item | text | status |", "|---|---|---|"]
for n, t, w in C:
    ok = t in w; bad += not ok; out.append(f"| {n} | {t[:60]} | {'OK' if ok else 'MISSING'} |")
for n, t in (("no placeholder", "AUTHOR TO CONFIRM"), ("no AI section", "Use of AI"), ("no VERIFY", "VERIFY")):
    ok = t not in M and t not in S; bad += not ok; out.append(f"| {n} | {t} | {'OK' if ok else 'FOUND'} |")
open(f"{P}/QC_number_check_r21_new.md", "w", encoding="utf-8").write("\n".join(out) + "\n"); print(f"{len(C) + 3} checks, {bad} missing"); sys.exit(1 if bad else 0)
