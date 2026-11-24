"""r10b (2026-11-07): trim the main-text Methods of submission_ready_20261107 to <3,000 words.
The full r10 text of five R6 Methods sections (Uniform-annotation validation; omega-class decomposition and functional
concentration; Per-origin RELAX; Lineage covariates; Post hoc robustness audits) moves verbatim into the SI
Supplementary Methods section "Uniform-annotation, omega-class and per-origin analyses (R6)"; main Methods gets the
condensed versions below, and the functional-module sentence moves from "Relative rates..." into "Post hoc analyses".
Usage (from paper/): python scripts/r10b_methods_sty.py submission_ready_20261107
Asserts the r10 inputs (main.md cbb89aa1, SI.md 11d9251b, SI.docx 462f3ee8), so it runs once."""
import re, sys, os, hashlib

PKG = sys.argv[1]
MD, SM, SD = (os.path.join(PKG, *p) for p in (("manuscript", "main.md"), ("supplementary", "Supplementary_Information.md"),
                                               ("supplementary", "Supplementary_Information.docx")))
md5f = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
assert md5f(MD).startswith("cbb89aa1") and md5f(SM).startswith("11d9251b") and md5f(SD).startswith("462f3ee8")
T = open(MD, encoding="utf-8").read()

NEW = {
"### Uniform-annotation validation": """The analyses in this and the next three sections (R6) were specified in a preregistration frozen by MD5 checksum on 31 October 2026, before any R6 fit existed, with an addendum for the per-origin analysis frozen the same day; departures are listed in an append-only record (Supplementary Data 2), and the full protocol is in Supplementary Methods. Two hundred orthogroups per clade were drawn at random from those with production alignments and real and null RELAX fits. Every NCBI- or author-annotated species was re-annotated with the production miniprot route, and original gene models were mapped to miniprot models by DIAMOND^57^ reciprocal best hit (identity ≥ 90%, query coverage ≥ 50%); 106 of 117 species were re-annotated and 78% of their genes mapped. Sequences without a mapped model were dropped from both arms, so original (O) and uniform (U) alignments contain the same species. RELAX was run on four fits per orthogroup (real and null foreground, O and U). The primary estimand was the log odds ratio of a relaxed call for real versus null foreground within U, from logistic regression with clade fixed effects and 95% CIs from 2,000 orthogroup-cluster bootstrap replicates within clade; success required a CI excluding zero and a ratio of the U to the O estimate of at least 0.5. Because RELAX in HyPhy 2.5.96 starts from a single random point, a prespecified random stratified subsample of 194 orthogroups was refitted to assess stability; repeat fits did not replace the primary fits.""",
"### ω-class decomposition and functional concentration": """For the 400 orthogroups per clade with real and null BUSTED-PH fits, the three ω classes of the unconstrained model were summarised for test (T) and background (B) branches as the proportion-weighted mean ω of classes with ω < 1 (w_pur) and the proportion of sites in classes with ω > 1 (p_pos). Δ_pur = log10(w_pur,T + 0.001) − log10(w_pur,B + 0.001) and Δ_pos = p_pos,T − p_pos,B were averaged as real minus null within clades and then across clades with equal weight. Relaxed purifying selection was predicted if the excess Δ_pur had a CI above zero and the same sign in at least six of eight clades. Functional concentration was tested on all production fits with Gene Ontology annotation, using relaxed ~ label × in-set + clade for prespecified light, circadian and pigmentation sets (terms in Supplementary Methods), with orthogroup-cluster bootstrap CIs and Holm correction. Rates in Fig. 2 are means over clades of per-clade fractions in this production sample against its first null foreground; fits without P or K (eight real and two null) were counted as neither relaxed nor intensified.""",
"### Per-origin RELAX": """Each orthogroup of the uniform-annotation set was refitted with RELAX in group mode on its U alignment, with each deep lineage of the real tree, or each shallow pseudo-lineage of the null tree, as its own group. Groups were mapped to canonical lineages, the maximal all-deep (or all-shallow) subtrees of the species tree within the clade, and groups spanning more than one were excluded. For each canonical lineage present in at least five orthogroups, the lineage effect was its median log10 K (K clamped to 0.001–50) minus the median over the null pseudo-lineages of its clade. The pooled effect was the mean of clade means, clades with at least two lineages weighted equally, with CIs from 2,000 bootstrap replicates of orthogroups within clade. Recurrence was claimed only if the pooled effect had a CI below zero, more than half of lineages were below zero and the effect stayed below zero without the two clades with most lineages. Because these fits had been seen before the freeze, the analysis was repeated blind on up to 200 new orthogroups per clade from the production codon set, drawn with a recorded seed and fitted on their production original-route alignments, the only design difference, and the same three conditions were applied. The blind sample re-draws orthogroups rather than origins and contained 81 lineages, 16 more than the first sample.""",
"### Lineage covariates": """Covariates were computed from species data alone: the mean log10 median depth of the lineage's species, log10 stem age, a codon-usage proxy for effective population size (the relative deficit of the effective number of codons against its GC3 expectation) and a chemoautosynthesis indicator. Within each sample, lineage effects were regressed on clade and the z-scored continuous covariates; P values came from 10,000 within-clade permutations of the focal covariate, with Holm correction over the three focal covariates, and intervals from the lineage-effect bootstrap.""",
"### Post hoc robustness audits": """Six analyses were added after the R6 results were known and are post hoc; for four of them the procedure was recorded or frozen by MD5 checksum before any statistic was computed (Supplementary Data 3; full protocols in Supplementary Methods). They assessed (i) fit failures, by reassigning failed fits over a 6 × 6 grid of outcome probabilities; (ii) assembly and annotation quality, by within-clade permutation tests of lineage effects against seven quality measures; (iii) cross-clade recurrence, against a null that keeps each clade's number of relaxed calls; (iv) substitution saturation, using the GTR branch lengths of each production fit as divergence; (v) a cross-clade proxy for shallow-water purifying constraint as a predictor of relaxation; and (vi) relaxation in three cellular modules (membrane and ion homeostasis, proteostasis, and energy metabolism and redox), measured as the difference in log10 K between deep and null branches of module members against the other annotated orthogroups, pooled over clades and repeated on the blind sample (Supplementary Note 9).""",
}
B = T.split("\n\n")
moved = []
for h, body in NEW.items():
    k = [j for j, b in enumerate(B) if b.startswith(h + "\n")]
    assert len(k) == 1, h
    moved.append(B[k[0]])
    B[k[0]] = (h if h != "### Post hoc robustness audits" else "### Post hoc analyses") + "\n" + body
jm = [j for j, b in enumerate(B) if b.startswith("### Relative rates, functions and gene families\n")]; assert len(jm) == 1
i = B[jm[0]].index(" In a post hoc test whose modules"); B[jm[0]] = B[jm[0]][:i]
assert B[jm[0]].endswith("random sets of shallow fishes.")
T = "\n\n".join(B)
# references: still all cited, in first-citation order
body = T.split("## References")[0]
seq = []
for mm in re.finditer(r"\^([\d,–]+)\^", body):
    for part in mm.group(1).split(","):
        if "–" in part:
            x, y = map(int, part.split("–")); seq += list(range(x, y + 1))
        else:
            seq.append(int(part))
seq += [int(x) for x in re.findall(r"ref\. (\d+)", body)]
first = []
[first.append(n) for n in seq if n not in first]
assert first == list(range(1, 59)), first
wc = lambda x: len(re.findall(r"\S+", re.sub(r"\^[^^]*\^", "", x)))
nm = wc(T.split("## Methods")[1].split("## Data availability")[0])
assert nm < 3000, nm
open(MD, "w", encoding="utf-8", newline="\n").write(T)

# SI: the full sections go under the R6 Supplementary Methods heading
S = open(SM, encoding="utf-8").read()
OLD1 = "These analyses are described in full in Methods (Uniform-annotation validation; ω-class decomposition and functional concentration; Per-origin RELAX). "
NEW1 = "The full protocols of the R6 analyses and of the post hoc analyses (condensed in Methods) follow. "
assert S.count(OLD1) == 1
S = S.replace(OLD1, NEW1)
head = "### Uniform-annotation, ω-class and per-origin analyses (R6)"
i = S.index(head); j = S.index("\n## Supplementary Table legends", i)
sec = S[i:j].rstrip("\n")
add = []
for blk in moved:
    h, txt = blk.split("\n", 1)
    add.append("*%s.* %s" % (h[4:], txt))
S = S[:i] + sec + "\n\n" + "\n\n".join(add) + "\n" + S[j:]
open(SM, "w", encoding="utf-8", newline="\n").write(S)

from docx import Document
TOK = re.compile(r"(\*\*.+?\*\*|\*[^*\s][^*]*?\*|\^[^^\s][^^]*?\^|\b[A-Za-z]_(?:paired|real|null|deep|shallow|pur|pos|[A-Z]{1,2})\b)")


def add_runs(p, s):
    for part in TOK.split(s):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            r = p.add_run(part[2:-2]); r.bold = True
        elif part.startswith("^") and part.endswith("^") and len(part) > 2:
            r = p.add_run(part[1:-1]); r.font.superscript = True
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            r = p.add_run(part[1:-1]); r.italic = True
        elif re.fullmatch(r"[A-Za-z]_\w+", part):
            x, y = part.split("_", 1)
            r = p.add_run(x); r.italic = True
            r = p.add_run(y); r.font.subscript = True
        else:
            p.add_run(part)


doc = Document(SD)
hit = [p for p in doc.paragraphs if p.text.startswith("These analyses are described in full in Methods")]
assert len(hit) == 1
p0 = hit[0]
for r in p0.runs:
    if OLD1 in r.text:
        r.text = r.text.replace(OLD1, NEW1); break
else:
    raise SystemExit("sentence split over runs")
anchor = p0._p
for a in add:
    el = doc.add_paragraph(); add_runs(el, a); anchor.addnext(el._p); anchor = el._p
doc.save(SD)
print("Methods words", nm, "| main", md5f(MD)[:8], "SI.md", md5f(SM)[:8], "SI.docx", md5f(SD)[:8])
