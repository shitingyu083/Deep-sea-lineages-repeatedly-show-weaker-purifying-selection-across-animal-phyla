"""Splice unchanged Methods, statistics and back matter of the 20261029 main.md into main_r6_20261101.md.
Old ED5-8 -> ED6-9 is not needed here (old Methods cite no figures). ape ref 49 -> 50; DIAMOND = new ref 49
(first cited in the R6 Methods subsection, which precedes the Statistics paragraph)."""
import io, re, sys
OLD = r"./submission_ready_20261029/manuscript/main.md"
NEW = r"./build/main_r6_20261101.md"
old = io.open(OLD, encoding="utf-8").read().split("\n")
new = io.open(NEW, encoding="utf-8").read()


def block(a, b):                                    # 1-based inclusive line range of old main.md
    return "\n".join(old[a - 1:b])


assert old[64] == "## Methods" and old[101].startswith("### Statistics") and old[104] == "## Data availability"
methods = block(67, 101).rstrip()
stats = block(102, 103)
assert "ape 5.8.1^49^" in stats
stats = stats.replace("ape 5.8.1^49^", "ape 5.8.1^50^").replace(
    "TreeShrink 1.4.0,", "TreeShrink 1.4.0, DIAMOND [AUTHOR TO CONFIRM: version],")
back = block(105, 177)
back = back.replace("Source data for Figs. 1–4 and Extended Data Figs. 1–8 are provided with this paper.",
                    "Source data for Figs. 1–4 and Extended Data Figs. 1–9 are provided with this paper. Per-fit results "
                    "and statistics of the uniform-annotation, ω-class, functional-concentration and per-origin "
                    "analyses are in Supplementary Data 2.")
assert "Supplementary Data 2" in back
back = back.replace("[References 34–49 are cited only in Methods.", "[References 34–50 are cited only in Methods.")
ref49 = [l for l in back.split("\n") if l.startswith("49. Paradis")]
assert len(ref49) == 1
back = back.replace(ref49[0], "49. Buchfink, B., Reuter, K. & Drost, H.-G. Sensitive protein alignments at tree-of-life "
                    "scale using DIAMOND. *Nat. Methods* **18**, 366–368 (2021).\n" + "50." + ref49[0][3:])
for tag, txt in (("<!-- OLD-METHODS -->", methods), ("<!-- OLD-STATS -->", stats + "\n\n" + back)):
    assert new.count(tag) == 1
    new = new.replace(tag, txt)
io.open(NEW, "w", encoding="utf-8", newline="\n").write(new)
print("ok", len(new.split()))
