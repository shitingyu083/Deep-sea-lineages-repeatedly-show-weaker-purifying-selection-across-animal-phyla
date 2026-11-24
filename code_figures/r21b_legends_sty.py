"""r21b: Fig. 2-4 legends trimmed below the 250-word Nature limit; Fig. 4 title made descriptive (guide: a brief title
sentence, not the result). No numbers removed. Run once after r21_text_sty.py: python scripts/r21b_legends_sty.py <pkg>"""
import sys, re, hashlib
P = sys.argv[1]; M = f"{P}/manuscript/main.md"
m = open(M, encoding="utf-8").read()
R = [
 ("and of 400 of them on a matched shallow null foreground of the same trees (bottom)", "and of 400 of them on a matched shallow null foreground (bottom)"),
 ("Fits with P < 0.05 are coloured by clade and others are grey;", "Fits with P < 0.05 in clade colour, others grey;"),
 ("(P < 0.05 and K > 1; right); the text gives their mean over clades.", "(P < 0.05 and K > 1; right)."),
 ("densities are drawn between −1.5 and 1.5, bars mark clade medians and the dashed line marks no difference.",
  "densities drawn between −1.5 and 1.5; bars, clade medians; dashed line, no difference."),
 ("resamples of orthogroups within clade, keeping the deep and null fits of an orthogroup together; numbers,",
  "resamples of orthogroups within clade, with paired fits kept together; numbers,"),
 ("Clade colours are the same in Figs 2–4; Annelida and Chondrichthyes are the analysis groups marked † in Fig. 1.",
  "Clade colours are the same in Figs 2–4; † groups as in Fig. 1."),
 ("the 81 of a second sample of orthogroups, which re-draws genes rather than colonization events, and the 83",
  "the 81 of a second sample of orthogroups and the 83"),
 ("**Fig. 4 | No robust excess of orthologues shared between clades after accounting for clade-specific selection rates.**",
  "**Fig. 4 | Reuse of relaxed and positively selected orthologues across clades.**"),
 ("The BUSTED-PH tests are secondary; for two or more clades the panel also gives the Holm-adjusted rate-conditioned P across the three secondary tests (0.063).",
  "For BUSTED-PH (secondary tests), the panel also gives the Holm-adjusted rate-conditioned P for two or more clades (0.063)."),
 ("RELAX power is computed at the nominal P equivalent to q < 0.05, stricter than the P < 0.05 calls, and is therefore conservative (Methods).",
  "RELAX power uses the stricter threshold described in Methods."),
]
for a, b in R:
    n = m.count(a); assert n == 1, (n, a[:70]); m = m.replace(a, b)
open(M, "w", encoding="utf-8", newline="\n").write(m)
leg = m[m.index("## Figure legends"):]
for p in re.findall(r"\*\*(?:Extended Data )?Fig\. \d+ \|.*", leg):
    n = len(p.split()); print(p[:26], n); assert n < 250
print("main", hashlib.md5(m.encode()).hexdigest()[:8])
