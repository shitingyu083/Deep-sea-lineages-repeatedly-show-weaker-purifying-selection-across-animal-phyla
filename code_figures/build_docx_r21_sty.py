"""r21 Word builder for Nature Ecology & Evolution submission.
main:  python build_docx_r21_sty.py main  <pkg>   manuscript/main.md -> manuscript/main.docx
cover: python build_docx_r21_sty.py cover <pkg>   cover_letter.md    -> cover_letter.docx
si:    python build_docx_r21_sty.py si    <pkg>   restyles supplementary/Supplementary_Information.docx in place
Format: A4, 2.5 cm margins, Times New Roman 12 pt, black headings, double-spaced text (cover letter single),
continuous line numbers and page numbers (manuscript), each main figure embedded on its own page with its legend.
Markdown subset as build_docx_r10_sty.py: #/##/### headings, **bold**, *italic*, ^superscript^, x_T subscripts."""
import re, sys, os, io
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_LINE_SPACING, WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from PIL import Image

FONT = "Times New Roman"
TOK = re.compile(r"(\*\*.+?\*\*|\*[^*\s][^*]*?\*|\^[^^\s][^^]*?\^|\b[A-Za-zΔ]_(?:paired|real|null|deep|shallow|pur|pos|[A-Z]{1,2})\b)")
STAR = {chr(0xE000 + n): "*" * n for n in (3, 2, 1)}


def _protect(s):   # significance stars before P (***P < 0.001) and the corresponding-author asterisk stay literal
    s = re.sub(r"(?<![*\w])(\*{1,3})(?=P <)", lambda m: chr(0xE000 + len(m.group(1))), s)
    return re.sub(r"(?<=,)\*(?=\^)|^\*(?=Corresponding)", "", s)


def _unprotect(x):
    for k, v in STAR.items():
        x = x.replace(k, v)
    return x


def add_runs(p, s, size=None):
    def run(t):
        r = p.add_run(_unprotect(t))
        if size: r.font.size = Pt(size)
        return r
    for part in TOK.split(_protect(s)):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            run(part[2:-2]).bold = True
        elif part.startswith("^") and part.endswith("^") and len(part) > 2:
            run(part[1:-1]).font.superscript = True
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            run(part[1:-1]).italic = True
        elif re.fullmatch(r"[A-Za-zΔ]_\w+", part):
            a, b = part.split("_", 1)
            run(a).italic = a != "Δ"
            run(b).font.subscript = True
        else:
            run(part)


def set_font(style, size, bold=None, italic=None):
    f = style.font; f.name = FONT; f.size = Pt(size); f.color.rgb = RGBColor(0, 0, 0)
    if bold is not None: f.bold = bold
    if italic is not None: f.italic = italic
    rpr = style.element.get_or_add_rPr(); rf = rpr.find(qn("w:rFonts"))
    if rf is None: rf = OxmlElement("w:rFonts"); rpr.insert(0, rf)
    for k in list(rf.attrib): del rf.attrib[k]           # drop theme fonts (Calibri/Cambria) from built-in styles
    for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"): rf.set(qn(a), FONT)


def setup(doc, spacing):
    st = doc.styles["Normal"]; set_font(st, 12)
    pf = st.paragraph_format; pf.line_spacing_rule = spacing; pf.space_after = Pt(0); pf.space_before = Pt(0)
    for name, size in (("Title", 14), ("Heading 1", 12), ("Heading 2", 12), ("Heading 3", 12)):
        s = doc.styles[name]; set_font(s, size, bold=True, italic=False)
        s.paragraph_format.line_spacing_rule = spacing; s.paragraph_format.space_before = Pt(12 if name != "Title" else 0)
        s.paragraph_format.space_after = Pt(0); s.paragraph_format.keep_with_next = True
        ppr = s.element.get_or_add_pPr()
        for b in ppr.findall(qn("w:pBdr")): ppr.remove(b)   # Title style carries a bottom border
    for s in doc.sections:
        s.page_width, s.page_height = Cm(21.0), Cm(29.7)
        s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Cm(2.5)


def page_numbers(doc):
    for s in doc.sections:
        f = s.footer.paragraphs[0]; f.alignment = WD_ALIGN_PARAGRAPH.CENTER
        fld = OxmlElement("w:fldSimple"); fld.set(qn("w:instr"), "PAGE"); f._p.append(fld)


def line_numbers(doc):
    for s in doc.sections:
        ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous")
        ln.set(qn("w:distance"), "283"); s._sectPr.append(ln)


SVG_EXT = "{96DAC541-7B7A-43D3-8B79-37D633B846F1}"
SVG_NS = "http://schemas.microsoft.com/office/drawing/2016/SVG/main"


def add_svg(run_inline, doc, pdf_path):
    """r22: attach a vector SVG (from the figure PDF) to an inline picture as Word's svgBlip, keeping the PNG as the
    fallback for readers without SVG support. Word 2016+/365 displays and prints the vector version."""
    import pymupdf
    from docx.opc.part import Part
    from docx.opc.packuri import PackURI
    from docx.opc.constants import RELATIONSHIP_TYPE as RT
    svg = pymupdf.open(pdf_path)[0].get_svg_image(text_as_path=True).encode("utf-8")
    pkg = doc.part.package
    n = 1 + sum(1 for p in pkg.iter_parts() if str(p.partname).endswith(".svg"))
    part = Part(PackURI(f"/word/media/vector{n}.svg"), "image/svg+xml", svg, pkg)
    rid = doc.part.relate_to(part, RT.IMAGE)
    blip = run_inline._inline.xpath(".//a:blip")[0]
    ext_lst = OxmlElement("a:extLst"); ext = OxmlElement("a:ext"); ext.set("uri", SVG_EXT)
    from lxml import etree
    sb = etree.SubElement(ext, "{%s}svgBlip" % SVG_NS, nsmap={"asvg": SVG_NS})
    sb.set(qn("r:embed"), rid)
    ext_lst.append(ext); blip.append(ext_lst)


def picture(doc, path, width_cm, pdf_path=None):
    """r22: full-resolution PNG (no downsampling) plus the vector SVG from pdf_path when given."""
    p = doc.add_paragraph(); p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    shp = p.add_run().add_picture(path, width=Cm(width_cm))
    if pdf_path:
        add_svg(shp, doc, pdf_path)
    return shp


def blocks(text):
    for block in re.split(r"\n\s*\n", text.strip()):
        lines = [l for l in block.split("\n") if l.strip()]
        if lines: yield lines


def build_main(P):
    text = re.sub(r"<!--.*?-->", "", open(f"{P}/manuscript/main.md", encoding="utf-8").read(), flags=re.S)
    doc = Document(); setup(doc, WD_LINE_SPACING.DOUBLE); line_numbers(doc)
    section = None
    for lines in blocks(text):
        if lines[0].startswith("#"):
            lvl = len(lines[0]) - len(lines[0].lstrip("#")); title = lines[0].lstrip("#").strip()
            if lvl == 2 and title in ("Methods", "References", "Figure legends", "Extended Data figure legends"):
                doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
            section = title if lvl <= 2 else section
            add_runs(doc.add_heading(level=0 if lvl == 1 else min(lvl - 1, 3)), title)
            lines = lines[1:]
            if not lines: continue
        if section == "Deep-sea lineages repeatedly show weaker purifying selection across animal phyla" and all(re.match(r"^\d+\. ", l) for l in lines):
            for l in lines:      # affiliations: superscript number + address
                n, rest = l.split(". ", 1); p = doc.add_paragraph(); r = p.add_run(n); r.font.superscript = True; add_runs(p, rest)
            continue
        while lines and not re.match(r"^\d+\. ", lines[0]) and len(lines) > 1 and all(re.match(r"^\d+\. ", l) for l in lines[1:]):
            add_runs(doc.add_paragraph(), lines[0].strip()); lines = lines[1:]
        if lines and all(re.match(r"^\d+\. ", l) for l in lines):
            for l in lines: add_runs(doc.add_paragraph(), l)
            continue
        para = " ".join(l.strip() for l in lines)
        mfig = re.match(r"^\*\*Fig\. (\d) \|", para)
        if section == "Figure legends" and mfig:
            k = int(mfig.group(1))
            if k > 1: doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
            picture(doc, f"{P}/figures/Fig{k}.png", 15.0, pdf_path=f"{P}/figures/Fig{k}.pdf")
            p = doc.add_paragraph(); p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            p.paragraph_format.space_before = Pt(6); add_runs(p, para, size=11)
            continue
        add_runs(doc.add_paragraph(), para)
    page_numbers(doc)
    cp = doc.core_properties; cp.title = text.split("\n", 1)[0].lstrip("# ").strip(); cp.author = "Tingyu Shi"
    cp.last_modified_by = "Tingyu Shi"; cp.comments = ""; cp.keywords = ""; cp.subject = ""
    doc.save(f"{P}/manuscript/main.docx"); print("wrote main.docx")


def build_cover(P):
    text = open(f"{P}/cover_letter.md", encoding="utf-8").read()
    doc = Document(); setup(doc, WD_LINE_SPACING.SINGLE)
    for lines in blocks(text):
        p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(10); p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
        add_runs(p, " ".join(l.strip() for l in lines))
    cp = doc.core_properties; cp.title = "Cover letter"; cp.author = "Hao Song"; cp.last_modified_by = "Hao Song"; cp.comments = ""
    doc.save(f"{P}/cover_letter.docx"); print("wrote cover_letter.docx")


def restyle_si(P):
    f = f"{P}/supplementary/Supplementary_Information.docx"; doc = Document(f)
    setup(doc, WD_LINE_SPACING.ONE_POINT_FIVE)
    for name in ("Caption", "Table Grid"):
        if name in [s.name for s in doc.styles] and doc.styles[name].type == 1: set_font(doc.styles[name], 10)
    for t in doc.tables:      # table text 9 pt, single spaced, so the inserted tables keep their layout
        for row in t.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
                    for r in p.runs: r.font.size = Pt(9); r.font.name = FONT
    for p in doc.paragraphs:  # pictures single spaced so they are not clipped by an exact/1.5 line height
        if p._p.xpath(".//w:drawing"): p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    cp = doc.core_properties; cp.author = "Tingyu Shi"; cp.last_modified_by = "Tingyu Shi"; cp.comments = ""
    doc.save(f); print("restyled SI.docx")


if __name__ == "__main__":
    mode, P = sys.argv[1], sys.argv[2]
    {"main": build_main, "cover": build_cover, "si": restyle_si}[mode](P)


def si_contents(P):
    """r21: the SI 'Contents' list was one run-on paragraph of ' - ' items; split into one line per item (docx and md)."""
    import copy
    f = f"{P}/supplementary/Supplementary_Information.docx"; doc = Document(f)
    for p in doc.paragraphs:
        if p.text.startswith("- Supplementary Notes 1"):
            items = [x.strip() for x in re.split(r"(?:^|\s)-\s(?=[A-Z])", p.text) if x.strip()]
            for it in reversed(items):
                q = copy.deepcopy(p._p); p._p.addnext(q)
            nxt = p._p.getnext()
            for it in items:
                from docx.text.paragraph import Paragraph
                para = Paragraph(nxt, p._parent)
                for r in para.runs[1:]: r._r.getparent().remove(r._r)
                para.runs[0].text = it; nxt = nxt.getnext()
            p._p.getparent().remove(p._p); print("contents items", len(items)); break
    doc.save(f)
