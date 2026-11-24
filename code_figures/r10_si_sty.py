"""r10 (2026-11-07): Supplementary Information for submission_ready_20261107.
md: python r10_si_sty.py md  <r9 SI.md (md5 asserted)> <dst SI.md>
docx: python r10_si_sty.py docx <r10 SI.md> <SI.docx, patched in place; r9 docx md5 asserted>
Changes: Extended Data Fig. 9 -> 10 and 8 -> 9 (renumbering), main Fig. 3e -> Extended Data Fig. 8e, Contents
"Notes 1-9", Supplementary Note 9 (functional modules; F1, numbers from work/addon/F1_modules/out/F1_report.txt
md5 22232a0c), Supplementary Fig. 2 legend points to Fig. 4b,c, Supplementary Data 3 legend adds F1 files."""
import re, sys, hashlib

mode, a, b = sys.argv[1:4]
md5f = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()

NOTE9 = """### Supplementary Note 9 | Functional modules

This test asked whether different orthogroups in different clades point to the same cellular functions. It was added after the main results were known and is post hoc; its modules, statistic and decision rule were frozen by MD5 checksum before any module statistic was computed, and the one change, recorded before computation, kept an obsolete Gene Ontology identifier (GO:0051082) and added the terms that the ontology lists for it (Supplementary Data 3). The modules were membrane and ion homeostasis (M1: ion transmembrane transport, phospholipid and fatty-acid metabolism, lipid transport, membrane organization, Na^+^/K^+^-ATPase activity and response to osmotic stress), proteostasis (M2: protein folding, chaperones, the endoplasmic-reticulum unfolded-protein response, ER-associated degradation, ubiquitin-dependent catabolism, the proteasome and autophagy) and energy metabolism and redox (M3: oxidative phosphorylation, electron transport, ATP synthesis, the tricarboxylic acid cycle, fatty-acid β-oxidation and redox homeostasis); light, circadian, pigmentation and ribosome sets were controls. For each orthogroup, d was log10 K on deep minus log10 K on matched null branches in the first sample (production fits, 400 orthogroups per clade with both fits) and the difference between deep-lineage and null pseudo-lineage medians in the blind sample. In each clade, θ was the mean d of module members minus that of the other annotated orthogroups, so θ < 0 means relaxation beyond the genome-wide excess. Clades with at least ten members were testable; clades were pooled with equal weights and tested by 10,000 within-clade permutations of module labels (one-sided). A module counted as supported only if the first-sample pooled θ was below zero at P < 0.0167, θ was below zero in at least five of at least six testable clades, every leave-one-clade-out estimate was below zero at P < 0.05, and the blind-sample θ was below zero at P < 0.05.

No module met the rule. For M1, θ was 0.006 in the first sample (P = 0.57; four of eight clades below zero) and −0.046 in the blind sample (P = 0.10). For M2, θ was 0.044 in the first sample (P = 0.86; two of eight clades below zero) and −0.105 in the blind sample (P = 0.012; five of six testable clades below zero), so the two samples disagree in sign. M3 was testable in two clades in the first sample and in none in the blind sample. Among the controls, light-related orthogroups gave θ = 0.065 (P = 0.92), and the circadian and pigmentation sets were too small to test. Permuting within deciles of background purifying ω left the first-sample results unchanged. Relaxed module orthogroups were mostly clade-specific: of 111 metazoan orthologous groups with a relaxed M1 member, 107 were relaxed in one clade only, and of 71 for M2, 66."""

DATA3_ADD = (" The archive also holds the frozen plan, deviation record, script and outputs of the functional-module test "
             "(Supplementary Note 9; folder F1_modules).")
SF2_ADD = " Panels a and d are redrawn in more detail in Fig. 4b,c of the main paper."


def renum(t):
    n9 = t.count("Extended Data Fig. 9"); n8 = t.count("Extended Data Fig. 8")
    t = t.replace("Extended Data Fig. 9", "Extended Data Fig. 10").replace("Extended Data Fig. 8", "Extended Data Fig. 9")
    t2 = re.sub(r"(?<!Data )(?<!Supplementary )Fig\. 3e", "Extended Data Fig. 8e", t)
    return t2, n9, n8, t2 != t


if mode == "md":
    assert md5f(a).startswith("e625f373"), md5f(a)
    T = open(a, encoding="utf-8").read()
    T, n9, n8, f3 = renum(T)
    assert (n9, n8, f3) == (1, 5, True), (n9, n8, f3)
    assert T.count("- Supplementary Notes 1–8") == 1
    T = T.replace("- Supplementary Notes 1–8", "- Supplementary Notes 1–9")
    assert T.count("\n## Supplementary Methods") == 1
    T = T.replace("\n## Supplementary Methods", "\n" + NOTE9 + "\n\n## Supplementary Methods")
    m = re.search(r"\*\*Supplementary Data 3 \|\*\*[^\n]*", T); assert m
    T = T[:m.end()] + DATA3_ADD + T[m.end():]
    m = re.search(r"\*\*Supplementary Fig\. 2 \|[^\n]*", T); assert m
    T = T[:m.end()] + SF2_ADD + T[m.end():]
    open(b, "w", encoding="utf-8", newline="\n").write(T)
    print("md", md5f(b))
else:
    from docx import Document
    sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
    assert md5f(b).startswith("fe2b00b7"), md5f(b)
    blocks = [x.strip() for x in re.split(r"\n\s*\n", open(a, encoding="utf-8").read()) if x.strip()]
    doc = Document(b)
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

    nrep = 0
    for p in doc.paragraphs:          # renumbering inside runs (each target string sits in one run in the r9 docx)
        for r in p.runs:
            t, n9, n8, f3 = renum(r.text)
            if t != r.text:
                r.text = t; nrep += n9 + n8 + int(f3)
    full = "\n".join(p.text for p in doc.paragraphs)
    assert "Extended Data Fig. 8" not in full.replace("Extended Data Fig. 8e", "") , "renumber missed"
    assert nrep >= 6, nrep

    def one(pred):
        h = [p for p in doc.paragraphs if pred(p.text)]
        assert len(h) == 1, len(h)
        return h[0]

    c = one(lambda t: t.startswith("- Supplementary Notes 1–8")); c.runs[0].text = c.runs[0].text.replace("1–8", "1–9")
    assert c.text.startswith("- Supplementary Notes 1–9")
    one(lambda t: t.startswith("Supplementary Data 3 |")).add_run(DATA3_ADD)
    one(lambda t: t.startswith("Supplementary Fig. 2 |")).add_run(SF2_ADD)
    anchor = one(lambda t: t == "Supplementary Methods")._p
    for x in NOTE9.split("\n\n"):
        if x.startswith("### "):
            el = doc.add_heading(level=3); add_runs(el, x[4:].strip())
        else:
            el = doc.add_paragraph(); add_runs(el, x)
        anchor.addprevious(el._p)
    doc.save(b)
    print("docx", md5f(b), "renumbered", nrep)
