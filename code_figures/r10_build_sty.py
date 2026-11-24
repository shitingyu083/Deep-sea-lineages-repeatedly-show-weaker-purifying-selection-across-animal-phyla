"""r10 (2026-11-07): build submission_ready_20261107/manuscript/main.md from the r9 main.md (md5 b607fe02) and the
author-reviewed draft paper/r10_results_legends_draft_sty.md (md5 asserted below). Text sources: the draft only.
Changes: Abstract, Introduction (Main), Results, Discussion, legends of Figs 2-4 (v4 figures); the old Fig. 3 legend
becomes Extended Data Fig. 8 (old ED Figs 8 and 9 -> 9 and 10); Methods: one sentence on the functional-module test
(Supplementary Note 9) and "Extended Data Figs. 1-10". References are unchanged (asserted: every reference cited, in
first-citation order). Usage: python r10_build_sty.py <r9 main.md> <draft.md> <dst main.md>"""
import re, sys, hashlib

src, drf, dst = sys.argv[1:4]
md5 = lambda s: hashlib.md5(s.encode("utf-8")).hexdigest()
T = open(src, encoding="utf-8").read(); D = open(drf, encoding="utf-8").read()
assert md5(T).startswith("b607fe02"), md5(T)
_dm = hashlib.md5(open(drf, "rb").read()).hexdigest()
assert _dm.startswith(sys.argv[4] if len(sys.argv) > 4 else "bb48dd67"), _dm


def dsec(head, stop):
    i = D.index(head + "\n") + len(head) + 1
    j = D.index(stop, i)
    return D[i:j].strip()


ABS = dsec("## Abstract (r10)", "## Discussion (r10)")
DIS = dsec("## Discussion (r10)", "## Main (Introduction, r10)")
INTRO = dsec("## Main (Introduction, r10)", "## Results")
RES = dsec("## Results", "---\n\n## Figure legends")
LEG = {k: next(b for b in D.split("\n\n") if b.startswith("**Fig. %d |" % k)) for k in (2, 3, 4)}
ED8_TITLE = "**Extended Data Fig. 8 | Uniform-annotation validation and light-related gene sets.**"
for x in (ABS, DIS, INTRO, RES):
    assert "r10" not in x and "DRAFT" not in x


def between(t, a, b):
    i = t.index(a + "\n") + len(a) + 1
    j = t.index(b, i)
    return i, j


def put(a, b, body):
    global T
    assert T.count(a + "\n") == 1, a
    i, j = between(T, a, b)
    T = T[:i] + "\n" + body + "\n\n" + T[j:]


put("## Abstract", "## Main", ABS)
put("## Main", "## Results", INTRO)
put("## Results", "## Discussion", RES)
put("## Discussion", "## Methods", DIS)
T = T.replace("\n\n\n", "\n\n")

B = T.split("\n\n")


def one(pred):
    k = [j for j, b in enumerate(B) if pred(b)]
    assert len(k) == 1, len(k)
    return k[0]


j3 = one(lambda b: b.startswith("**Fig. 3 |"))
old3 = B[j3]
m = re.match(r"\*\*Fig\. 3 \|[^*]*\*\*", old3); assert m
ED8 = ED8_TITLE + old3[m.end():]
for k in (2, 3, 4):
    B[one(lambda b, k=k: b.startswith("**Fig. %d |" % k))] = LEG[k]
j9 = one(lambda b: b.startswith("**Extended Data Fig. 9 |"))
j8 = one(lambda b: b.startswith("**Extended Data Fig. 8 |"))
B[j9] = B[j9].replace("**Extended Data Fig. 9 |", "**Extended Data Fig. 10 |", 1)
B[j8] = B[j8].replace("**Extended Data Fig. 8 |", "**Extended Data Fig. 9 |", 1)
B.insert(j8, ED8)
# Methods
jm = one(lambda b: b.startswith("### Relative rates, functions and gene families\nRelative evolutionary rates"))
F1 = (" In a post hoc test whose modules, statistic and decision rule were frozen by MD5 checksum before it was run, "
      "three cellular modules (membrane and ion homeostasis, proteostasis, and energy metabolism and redox; fixed Gene "
      "Ontology term lists) were tested for relaxation beyond the genome-wide excess: for each orthogroup, the "
      "difference in log10 K between deep and matched null branches was averaged over module members and compared "
      "with the other annotated orthogroups of the clade, clades were pooled with equal weights, significance came "
      "from 10,000 within-clade permutations of module labels, and the test was repeated on the blind orthogroup "
      "sample (Supplementary Note 9).")
assert B[jm].endswith("random sets of shallow fishes.")
B[jm] = B[jm] + F1
jd = one(lambda b: "Source data for Figs. 1–4 and Extended Data Figs. 1–9" in b)
B[jd] = B[jd].replace("Extended Data Figs. 1–9", "Extended Data Figs. 1–10")
T = "\n\n".join(B)

# checks: no stale display references, every reference cited, first-citation order
body = T.split("## References")[0]
assert "Extended Data Fig. 11" not in body
CIT = re.compile(r"\^([\d,–]+)\^")
seq = []
for mm in CIT.finditer(body):
    for part in mm.group(1).split(","):
        if "–" in part:
            x, y = map(int, part.split("–")); seq += list(range(x, y + 1))
        else:
            seq.append(int(part))
seq += [int(x) for x in re.findall(r"ref\. (\d+)", body)]
first = []
for n in seq:
    if n not in first:
        first.append(n)
nref = len(re.findall(r"^\d+\. ", T.split("## References")[1].split("## Acknowledgements")[0], flags=re.M))
assert first == list(range(1, nref + 1)), [(i + 1, n) for i, n in enumerate(first) if n != i + 1][:5]
open(dst, "w", encoding="utf-8", newline="\n").write(T)
print("refs", nref, "md5", md5(T))
