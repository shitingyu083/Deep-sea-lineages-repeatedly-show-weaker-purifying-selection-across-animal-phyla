"""r10: as build_docx_sty.py, plus significance stars before P kept literal (ED5 / Fig. 2 legends);
r10c: Δ_pur, Δ_pos, w_pur and p_pos set with subscripts.
Build manuscript/main.docx from manuscript/main.md (python-docx).
Markdown subset: #/##/### headings, paragraphs, numbered reference lines, **bold**, *italic*, ^superscript^,
x_T style subscripts for trait symbols (f_T, f_F, n_paired ...), and HTML comments (dropped)."""
import re, sys, os
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_LINE_SPACING

root = sys.argv[1] if len(sys.argv) > 1 else "."
src = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "manuscript", "main.md")
dst = sys.argv[3] if len(sys.argv) > 3 else os.path.join(root, "manuscript", "main.docx")
text = re.sub(r"<!--.*?-->", "", open(src, encoding="utf-8").read(), flags=re.S)

doc = Document()
st = doc.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(11)
st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
st.paragraph_format.space_after = Pt(0)
for s in doc.sections:
    s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Cm(2.5)

TOK = re.compile(r"(\*\*.+?\*\*|\*[^*\s][^*]*?\*|\^[^^\s][^^]*?\^|\b[A-Za-zΔ]_(?:paired|real|null|deep|shallow|pur|pos|[A-Z]{1,2})\b)")


STAR = {"": "***", "": "**", "": "*"}   # r10: significance stars before P (e.g. ***P < 0.01)


def _protect(s):
    return re.sub(r"(?<![*\w])(\*{1,3})(?=P <)", lambda m: ""[:0] + chr(0xE000 + len(m.group(1))), s)


def _unprotect(x):
    for k, v in STAR.items():
        x = x.replace(k, v)
    return x


def add_runs(p, s):
    _add = p.add_run
    p.add_run = lambda t=None, *a, **k: _add(_unprotect(t) if t else t, *a, **k)
    for part in TOK.split(_protect(s)):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            r = p.add_run(part[2:-2]); r.bold = True
        elif part.startswith("^") and part.endswith("^") and len(part) > 2:
            r = p.add_run(part[1:-1]); r.font.superscript = True
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            r = p.add_run(part[1:-1]); r.italic = True
        elif re.fullmatch(r"[A-Za-zΔ]_\w+", part):
            a, b = part.split("_", 1)
            r = p.add_run(a); r.italic = a != "Δ"       # r10c: Δ_pur / Δ_pos upright, w_pur / p_pos italic
            r = p.add_run(b); r.font.subscript = True
        else:
            p.add_run(part)


# line numbers are requested by Nature at submission
sect = doc.sections[0]._sectPr
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); sect.append(ln)

for block in re.split(r"\n\s*\n", text.strip()):
    lines = [l for l in block.split("\n") if l.strip()]
    if not lines:
        continue
    if lines[0].startswith("#"):
        lvl = len(lines[0]) - len(lines[0].lstrip("#"))
        h = doc.add_heading(level=min(lvl, 3) if lvl > 1 else 0)
        add_runs(h, lines[0].lstrip("#").strip())
        lines = lines[1:]
    # reference lists: one paragraph per numbered line; otherwise join soft-wrapped lines
    while lines and not re.match(r"^\d+\. ", lines[0]) and len(lines) > 1 and all(re.match(r"^\d+\. ", l) for l in lines[1:]):
        add_runs(doc.add_paragraph(), lines[0].strip()); lines = lines[1:]
    if lines and all(re.match(r"^\d+\. ", l) for l in lines):
        for l in lines:
            p = doc.add_paragraph(); p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            p.paragraph_format.space_after = Pt(3); add_runs(p, l)
    elif lines:
        add_runs(doc.add_paragraph(), " ".join(l.strip() for l in lines))

for i, s in enumerate(doc.sections):
    f = s.footer.paragraphs[0]
    fld = OxmlElement("w:fldSimple"); fld.set(qn("w:instr"), "PAGE"); f._p.append(fld)
doc.save(dst)
print("wrote", dst)
