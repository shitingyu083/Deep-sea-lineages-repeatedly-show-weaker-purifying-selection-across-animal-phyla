"""Update supplementary/Supplementary_Information.md of the 20261101 package to the R6 figure plan.
1. Remap figure references (old -> new; token pass so no chained replacement):
   Fig. 1a,c -> Fig. 1a, Extended Data Fig. 3c; Fig. 1d -> Fig. 1b; Fig. 2x -> Supplementary Fig. 1x;
   Fig. 3x -> Extended Data Fig. 5x; Fig. 4x -> Supplementary Fig. 2x; Extended Data Fig. 5-8 -> 6-9.
2. Header and contents; route-pure paragraph marked as superseded by the uniform-annotation test.
3. Supplementary Note 7 (R6 results, numbers copied from R6 reports), Supplementary Methods pointer,
   Supplementary Figure legends (old main Fig. 2/4 legends) and Supplementary Data 2 legend.
Usage: python update_si_r6_sty.py <package_dir>"""
import os, re, sys

pkg = sys.argv[1]
p = os.path.join(pkg, "supplementary", "Supplementary_Information.md")
t = open(p, encoding="utf-8").read()
assert "Supplementary Note 7" not in t, "already updated"

# 1. references ---------------------------------------------------------------------------------
ED = {"5": "6", "6": "7", "7": "8", "8": "9"}
t = re.sub(r"Extended Data Fig\. ([5-8])", lambda m: "@@ED" + ED[m.group(1)], t)
assert t.count("Fig. 1a,c") == 1 and t.count("Fig. 1d") == 1
t = t.replace("Fig. 1a,c", "Fig. 1a, @@ED3c").replace("Fig. 1d", "Fig. 1b")
t = re.sub(r"(?<![A-Za-z] )Fig\. 2([a-d])", r"@@SF1\1", t)
t = re.sub(r"(?<![A-Za-z] )Fig\. 3([a-d])", r"@@ED5\1", t)
t = re.sub(r"(?<![A-Za-z] )Fig\. 4([a-d])", r"@@SF2\1", t)
t = t.replace("@@ED", "Extended Data Fig. ").replace("@@SF", "Supplementary Fig. ")
assert not re.search(r"(?<!Data )(?<!Supplementary )Fig\. [2-4]", t), "unmapped main-figure reference"

# 2. header, contents ---------------------------------------------------------------------------
old = "Figure and Extended Data figure numbers refer to the main paper."
assert old in t
t = t.replace(old, "Figure and Extended Data figure numbers refer to the main paper; Supplementary Figures 1 and 2 are "
              "given at the end of this document (files Supplementary_Fig1 and Supplementary_Fig2).")
old = "- Supplementary Notes 1–6\n"
assert old in t
t = t.replace(old, "- Supplementary Notes 1–7\n")
old = "- Legend for Supplementary Data 1 (file Supplementary_Data_1_per_gene_selection.tsv.gz)"
assert old in t
t = t.replace(old, "- Supplementary Figures 1 and 2 with legends\n" + old +
              "\n- Legend for Supplementary Data 2 (file Supplementary_Data_2_R6_results.zip)")

old = "cannot be decided from these data."
assert t.count(old) == 1
t = t.replace(old, old + " This post hoc analysis is superseded by the preregistered uniform-annotation test "
              "(Supplementary Note 7), in which every species was re-annotated by one route instead of being pruned.")

# 3. new sections -------------------------------------------------------------------------------
NOTE7 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "si_r6_text_sty.md"), encoding="utf-8").read()
note7, methods, sfig, sdata = [s.strip() for s in NOTE7.split("<!-- SPLIT -->")]
old = "## Supplementary Methods\n"
assert t.count(old) == 1
t = t.replace(old, note7 + "\n\n" + old)
old = "## Supplementary Table legends"
assert t.count(old) == 1
t = t.replace(old, methods + "\n\n" + old)
t = t.rstrip() + "\n\n" + sdata + "\n\n" + sfig + "\n"
open(p, "w", encoding="utf-8").write(t)
print("updated", p, len(t.split()), "words")
