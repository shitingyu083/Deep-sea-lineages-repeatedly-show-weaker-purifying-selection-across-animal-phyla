"""Shared helpers for the v2 figures: every panel is redrawn in the form of one published benchmark panel.

Benchmarks (local copies in stage_sty/_ref_sty; panel-by-panel map in paper/figure_benchmark_map.md):
  terr - Nature 2025, doi:10.1038/s41586-025-09722-4: Fig. 1 (time tree with collapsed clades, phylum bands,
         silhouettes and node boxes), Extended Data Fig. 8 (null histogram with observed line), Extended Data Fig. 9
         (annotated heatmap)
  cell - Cell 2025 deep-sea fishes (PMID 40054449): Fig. 4A (framed scatter, two series with OLS fits and bands),
         Fig. 4G/H (adjacent bars, condition box, square-patch legend)
  sack - Sackton et al. 2019, Science 364:74: Fig. 3A (null histogram), Fig. 3C (Manhattan), Fig. 4A (dots)
Style, palette and IO come from figcommon_sty (paper/scripts, imported read-only, no bytecode written); outputs go
to viz_v2_sty/out, never to the package.
"""
import os, re, sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "paper", "scripts"))
import numpy as np
from figcommon_sty import *  # noqa: F401,F403  palette, rcParams, rtsv, fnum, W, DATA, BUILD, MM
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.colors import to_rgb
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(V2, exist_ok=True)


def save_v2(fig, name):
    fig.savefig(os.path.join(V2, name + ".pdf"))
    fig.savefig(os.path.join(V2, name + ".png"), dpi=300)
    fig.savefig(os.path.join(V2, name + ".tif"), dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def write_facts(name, rows):
    with open(os.path.join(V2, name), "w", encoding="utf-8") as f:
        f.write("".join("\t".join(str(c) for c in r) + "\n" for r in rows))


# ---- colours ------------------------------------------------------------------------------------------------
# terr Fig. 1, sampled from the published PNG (scratch x06/x08/x09). The benchmark time bar starts at 635 Ma, so
# Tonian and Cryogenian use the ICS 2023 chart colours.
ERA = [("Neoproterozoic", "Neop.", 1000.0, 538.8, "#f4a858"), ("Paleozoic", "Paleozoic", 538.8, 251.9, "#80ac94"),
       ("Mesozoic", "Mesozoic", 251.9, 66.0, "#08b8e0"), ("Cenozoic", "Cen.", 66.0, 0.0, "#e0e440")]
PERIOD = [("Tonian", "Ton.", "", 1000.0, 720.0, "#FEBF4E"), ("Cryogenian", "Cry.", "", 720.0, 635.0, "#FECC5C"),
          ("Ediacaran", "Edi.", "", 635.0, 538.8, "#f4c878"), ("Cambrian", "Camb.", "Є", 538.8, 485.4, "#709c6c"),
          ("Ordovician", "Ord.", "O", 485.4, 443.8, "#00947c"), ("Silurian", "Sil.", "S", 443.8, 419.2, "#94c4b0"),
          ("Devonian", "Dev.", "D", 419.2, 358.9, "#c08c50"), ("Carboniferous", "Carb.", "C", 358.9, 298.9, "#48989c"),
          ("Permian", "Perm.", "P", 298.9, 251.9, "#e05c44"), ("Triassic", "Tri.", "T", 251.9, 201.4, "#744888"),
          ("Jurassic", "Jur.", "J", 201.4, 145.0, "#00a8e0"), ("Cretaceous", "Cret.", "K", 145.0, 66.0, "#68b46c"),
          ("Paleogene", "Palg.", "Pg", 66.0, 23.03, "#ec9464"), ("Neogene", "Ng.", "N", 23.03, 2.58, "#f8d82c"),
          ("Quaternary", "Q.", "Q", 2.58, 0.0, "#fcf88c")]
# phylum bands: white on the left to colour on the right (terr 1280-px frame x positions 150..930)
_BX = np.array([150, 200, 260, 330, 400, 480, 560, 640, 700, 760, 820, 900, 930], float)
BAND_POS = (_BX - _BX[0]) / (_BX[-1] - _BX[0])
BAND_BLUE = ["#fcfcfc", "#f4f8fc", "#f0f4f8", "#e8f0f4", "#e4ecf4", "#dce4f0", "#d4e0ec", "#ccdcec", "#c8d8e8",
             "#c4d8e8", "#bcd4e4", "#b8d0e4", "#b4d0e0"]
BAND_GREY = ["#fcfcfc", "#f4f4f4", "#f0f0f0", "#e8e8e8", "#e4e4e4", "#dcdcdc", "#d8d8d8", "#d0d0d0", "#cccccc",
             "#c8c8c8", "#c4c4c4", "#bcc0c0", "#b8bcbc"]
# collapsed-clade gradient: pale apex -> state colour at the base; terr's mid stop lies 40% of the way (apex
# #d0dcdc, mid #90b8a4, base #2c8450); the apex here is the same white blend of our deep colour.
TRI_APEX, TRI_BASE, TRI_MIDF = "#d0d9e2", DEEP, 0.40
HPD_FILL, HPD_EDGE, AGE_TXT = "#8cb0d0", "#222428", "#280000"
BOX_EDGE, ARROW, LEADER, OUTGRP = "#606868", "#606060", "#a0a0a0", "#909090"
SHAL_TXT, DEEP2, DEEP3 = "#a07020", "#5b84ad", "#8fb3d4"
# terr Extended Data Fig. 8: bars (benchmark #5b93c9 -> our secondary blue), observed line and text, grid
ED8_BAR, ED8_RED, ED8_GRID = DEEP2, "#e8262a", "#ebebeb"


def ramp(stops, pos, n=256):
    """n x 3 RGB ramp through colour stops at positions pos (0..1)."""
    rgb = np.array([to_rgb(c) for c in stops]); t = np.linspace(0, 1, n)
    return np.stack([np.interp(t, pos, rgb[:, k]) for k in range(3)], axis=1)


TRI_RAMP = ramp([TRI_APEX, tuple(np.array(to_rgb(TRI_APEX)) + TRI_MIDF * (np.array(to_rgb(TRI_BASE)) -
                 np.array(to_rgb(TRI_APEX)))), TRI_BASE], [0, 0.5, 1])


# ---- dated tree (LSD2 nexus: ages in Ma, CI_height = 95% interval) --------------------------------------------
class N2:
    __slots__ = ("name", "length", "children", "parent", "age", "ci", "tips", "x", "y", "id")

    def __init__(self, parent=None):
        self.name, self.length, self.children, self.parent = "", 0.0, [], parent
        self.age, self.ci, self.tips, self.x, self.y, self.id = None, None, None, 0.0, 0.0, None


_TOK = re.compile(r"\(|\)|,|;|\[&[^\]]*\]|:[-+0-9.eE]+|[^(),:;\[\]\s]+")


def parse_nex_tree(path):
    s = open(path, encoding="utf-8").read()
    t = re.search(r"tree\s+[^=]+=\s*(.*?;)", s, re.S | re.I).group(1)
    root = N2(); cur = root
    for tok in _TOK.findall(t):
        if tok == "(":
            c = N2(cur); cur.children.append(c); cur = c
        elif tok == ",":
            c = N2(cur.parent); cur.parent.children.append(c); cur = c
        elif tok == ")":
            cur = cur.parent
        elif tok == ";":
            break
        elif tok.startswith("[&"):
            d = re.search(r'(?<![A-Za-z_])date="?(-?[0-9.eE+]+)"?', tok)
            if d:
                cur.age = -float(d.group(1)) + 0.0
            ci = re.search(r"CI_height=\{([-0-9.eE+]+),([-0-9.eE+]+)\}", tok)
            if ci:
                cur.ci = (float(ci.group(1)), float(ci.group(2)))
        elif tok.startswith(":"):
            cur.length = float(tok[1:])
        else:
            cur.name = tok
    for n in postorder(root):
        n.tips = frozenset([n.name]) if not n.children else frozenset().union(*[c.tips for c in n.children])
    return root


def mrca(root, names):
    names = frozenset(names); best = root
    for n in postorder(root):
        if names <= n.tips and len(n.tips) < len(best.tips):
            best = n
    return best

# ---- terr Fig. 1 drawing helpers: one axes in mm, x to the right, y downwards ------------------------------
def mm_axes(fig, w, h):
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, w); ax.set_ylim(h, 0); ax.axis("off")
    return ax

def text_mm(ax, t):
    """Rendered width and height of a text artist in mm."""
    bb = t.get_window_extent(ax.figure.canvas.get_renderer())
    return bb.width / ax.figure.dpi * 25.4, bb.height / ax.figure.dpi * 25.4

def hgrad(ax, x0, x1, y0, y1, rgb, clip=None, z=1, n=160):
    """Horizontal ramp (n x 3 RGB, left to right) over x0..x1 and y0 (top)..y1 (bottom), optionally clipped.
    r22: drawn as n abutting vector strips (no embedded image), colours linearly interpolated from the ramp."""
    from matplotlib.collections import PatchCollection
    rgb = np.asarray(rgb, float); t = np.linspace(0, 1, len(rgb)); tc = (np.arange(n) + 0.5) / n
    cols = np.stack([np.interp(tc, t, rgb[:, k]) for k in range(3)], 1)
    xs = x0 + (x1 - x0) * np.arange(n + 1) / n; ov = (x1 - x0) / n * 0.6
    pc = PatchCollection([Rectangle((xs[k], min(y0, y1)), xs[k + 1] - xs[k] + (ov if k < n - 1 else 0), abs(y1 - y0))
                          for k in range(n)], facecolors=cols, edgecolors="none", lw=0, zorder=z)
    ax.add_collection(pc)
    if clip is not None:
        pc.set_clip_path(clip)
    return pc


def vector_mask(ax, path, x, y, color, w=None, h=None, ha="left", z=6, alpha=1.0):
    """r22: PhyloPic silhouette as a filled vector outline (contour of the PNG alpha at 0.5, holes kept),
    at its true aspect ratio and centred on y. Returns (x0, x1, h)."""
    from PIL import Image
    import contourpy
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    a = np.asarray(Image.open(path).convert("LA"))[:, :, 1].astype(float) / 255.0
    ys, xs = np.nonzero(a > 0.02); a = a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    a = np.pad(a, 1)
    H, W = a.shape; ar = (W - 2) / (H - 2)
    w = h * ar if w is None else w; h = w / ar
    x0 = x - w if ha == "right" else x - w / 2 if ha == "center" else x
    gen = contourpy.contour_generator(z=a, fill_type=contourpy.FillType.OuterCode)
    polys, codes = gen.filled(0.5, 2.0)
    sx, sy = w / (W - 2), h / (H - 2)
    verts = np.concatenate([np.c_[x0 + (p[:, 0] - 0.5) * sx, (y - h / 2) + (p[:, 1] - 0.5) * sy] for p in polys])
    pp = PathPatch(Path(verts, np.concatenate(codes)), fc=color, ec="none", lw=0, alpha=alpha, zorder=z)
    ax.add_patch(pp)
    return x0, x0 + w, h

def tri(ax, xa, ya, xb, y0, y1, fill, lw=0.5, z=3):
    """Collapsed clade: apex at the crown (xa, ya), base on the tip line x = xb from y0 to y1. Focal-state
    clades carry the pale-to-colour ramp, the others are outlines only (terr Fig. 1)."""
    p = Polygon([(xa, ya), (xb, y0), (xb, y1)], closed=True, fc="none", ec="black", lw=lw, joinstyle="miter",
                zorder=z + 0.2)
    ax.add_patch(p)
    if fill:
        hgrad(ax, xa, xb, y0, y1, TRI_RAMP, clip=p, z=z)
    return p

def tipsq(ax, xb, y, s, lw=0.4, z=3.5):
    """Single focal-state species: square on the tip end with the same ramp (terr Fig. 1, Onychophora)."""
    r = Rectangle((xb - s, y - s / 2), s, s, fc="none", ec="black", lw=lw, zorder=z + 0.2)
    ax.add_patch(r)
    hgrad(ax, xb - s, xb, y - s / 2, y + s / 2, TRI_RAMP, clip=r, z=z)
    return r

def band(ax, x0, x1, y0, y1, grey, z=0):
    """Phylum band, white on the left to blue or grey on the right (terr Fig. 1)."""
    hgrad(ax, x0, x1, y0, y1, ramp(BAND_GREY if grey else BAND_BLUE, BAND_POS), z=z)

def silhouette(ax, path, x, y, w=None, h=None, ha="left", z=6):
    """PhyloPic silhouette in black at its true aspect ratio, centred on y (r22: vector outline, see vector_mask).
    Give the width w or the height h in mm; returns (x0, x1, h)."""
    return vector_mask(ax, path, x, y, "black", w=w, h=h, ha=ha, z=z)

def arrow(ax, x0, x1, y, z=5):
    """Grey connector with a filled head from a focal lineage to its box (terr Fig. 1)."""
    ax.annotate("", xy=(x1, y), xytext=(x0, y), zorder=z, annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.22", color=ARROW, lw=0.5,
                                shrinkA=0, shrinkB=0))
def hpd_bar(ax, x_old, x_young, y, h=0.7, z=7):
    """Age-interval bar at a crown (terr Fig. 1); the age label is placed by the caller, clear of the branches.
    Returns the bar box (x0, y0, x1, y1) in mm."""
    ax.add_patch(Rectangle((x_old, y - h / 2), x_young - x_old, h, fc=HPD_FILL, ec=HPD_EDGE, lw=0.3, zorder=z))
    return x_old, y - h / 2, x_young, y + h / 2


def node_box(ax, x, ytop, w, title, row1, row2, hh=2.6, rh=2.6, fs=5.5, w1=None, z=6):
    """terr node box: white with BOX_EDGE frame and separators, header, three cells, two cells. Rows are lists of
    (text, colour); w1 optionally gives the row-1 cell widths (the key box is sized to its text). Returns height."""
    Hb = hh + 2 * rh
    ax.add_patch(Rectangle((x, ytop), w, Hb, fc="white", ec=BOX_EDGE, lw=0.5, zorder=z))
    for yy in (ytop + hh, ytop + hh + rh):
        ax.plot([x, x + w], [yy, yy], color=BOX_EDGE, lw=0.4, zorder=z + 0.1, solid_capstyle="butt")
    ax.text(x + w / 2, ytop + hh / 2, title, ha="center", va="center_baseline", fontsize=fs, zorder=z + 0.2)
    for row, y0, ws in ((row1, ytop + hh, w1), (row2, ytop + hh + rh, None)):
        ws = ws or [w / len(row)] * len(row)
        xs = x + np.concatenate([[0.0], np.cumsum(ws)])
        for k, (s, c) in enumerate(row):
            if k:
                ax.plot([xs[k]] * 2, [y0, y0 + rh], color=BOX_EDGE, lw=0.4, zorder=z + 0.1, solid_capstyle="butt")
            ax.text((xs[k] + xs[k + 1]) / 2, y0 + rh / 2, s, ha="center", va="center_baseline", fontsize=fs,
                    color=c, zorder=z + 0.2)
    return Hb


def key_icons(ax, x, y, tri_w=6.0, tri_h=3.0, sq=1.2, gap=1.0):
    """Key glyphs of terr Fig. 1: a collapsed-clade triangle and, to its right, a single-species square."""
    tri(ax, x, y, x + tri_w, y - tri_h / 2, y + tri_h / 2, fill=True, lw=0.5)
    tipsq(ax, x + tri_w + gap + sq, y, sq)


def axes_mm(fig, Wf, Hf, x0, y0, w, h):
    """Axes from a box in figure mm (x0, y0 = top-left corner, y downwards)."""
    return fig.add_axes([x0 / Wf, 1 - (y0 + h) / Hf, w / Wf, h / Hf])
# ---- terr Fig. 1 geological time bar -------------------------------------------------------------------------
ICS_AGE = [(720.0, "720"), (635.0, "635"), (538.8, "538.8"), (485.4, "485.4"), (443.8, "443.8"), (419.2, "419.2"),
           (358.9, "358.9"), (298.9, "298.9"), (251.9, "251.9"), (201.4, "201.4"), (145.0, "145.0"), (66.0, "66.0"),
           (23.03, "23.03"), (2.58, "2.58")]
ICS_FIRST = (538.8, 251.9, 66.0, 720.0, 2.58)  # era boundaries and both ends are placed before the others


_INK = {}


def ink_mm(s, fs):
    """Outline ink extents (mm) of a string about its baseline origin, y up: (x0, y0, x1, y1)."""
    k = (s, fs)
    if k not in _INK:
        p = TextPath((0, 0), s, size=fs, prop=FontProperties(family=plt.rcParams["font.family"])).get_extents()
        kk = 25.4 / 72.0
        _INK[k] = (p.x0 * kk, p.y0 * kk, p.x1 * kk, p.y1 * kk)
    return _INK[k]


def ink_text(ax, x, yb, s, fs, align="center", **kw):
    """Text on baseline yb whose ink is centred on x (align='center') or starts at x (align='left')."""
    i0, _, i1, _ = ink_mm(s, fs)
    x = x - 0.5 * (i0 + i1) if align == "center" else x - i0
    return ax.text(x, yb, s, ha="left", va="baseline", fontsize=fs, **kw)


def _fit_label(ax, xc, yb, width, cands, fs):
    """Longest candidate label that fits its cell (terr writes names in full and abbreviates narrow cells)."""
    for s in cands:
        if not s:
            continue
        t = ink_text(ax, xc, yb, s, fs, zorder=3)
        if text_mm(ax, t)[0] <= width - 0.2:
            return s
        t.remove()
    return ""


def ics_bar(ax, X, t_old, y0, rh=2.43, fs=5, scale_x=None, axis_from=None, gap=0.3):
    """terr Fig. 1 time bar under the tree: abutting era and period rows, a tick line with a tick at every boundary
    and 'Present', the boundary ages (Ma), a numeric axis (majors 50, minors 25 Myr) and a '50 Myr' scale bar.
    X maps Ma to mm; the numeric axis starts at axis_from (default t_old), as terr's starts at a round age; boundary
    ages keep >= 2*gap mm of ink between them. Vertical offsets are terr's, as fractions of its row height rh.
    Returns the (shown, full name) pairs of shortened labels, for the legend, and the levels used (mm)."""
    axis_from = t_old if axis_from is None else axis_from
    cap = ink_mm("H", fs)[3]
    used = []
    for rows, yt in ((ERA, y0), (PERIOD, y0 + rh)):
        yb = yt + 0.455 * rh + 0.5 * cap             # cap centre 0.045 rh above the row middle (terr)
        for r in rows:
            nm, a0, a1, col = r[0], r[-3], r[-2], r[-1]
            if a1 >= t_old:
                continue
            a0 = min(a0, t_old)
            xa, xb = X(a0), X(a1)
            ax.add_patch(Rectangle((xa, yt), xb - xa, rh, fc=col, ec="none", zorder=2))
            if nm == "Quaternary":                   # terr: 'Q.' starts at the cell's left edge and runs past it
                ink_text(ax, xa, yb, "Q.", fs, align="left", zorder=3)
                used.append(("Q.", nm))
                continue
            s = _fit_label(ax, (xa + xb) / 2, yb, xb - xa, r[:-3], fs)
            if s and s != nm:
                used.append((s, nm))
    yr = y0 + 2 * rh                                 # bottom of the period row: tick line top, boundary ticks
    ax.plot([X(t_old), X(0)], [yr + 0.25 * 25.4 / 72] * 2, color="black", lw=0.5, zorder=3, solid_capstyle="butt")
    for a in [a for a, _ in ICS_AGE if a < t_old] + [0.0]:
        ax.plot([X(a)] * 2, [yr, yr + 0.35 * rh], color="black", lw=0.4, zorder=3, solid_capstyle="butt")
    ink_text(ax, X(0) + 0.30 * rh, yr + 0.13 * rh + 0.5 * cap, "Present", fs, align="left")
    yb = yr + 0.67 * rh + ink_mm("0", fs)[3]         # boundary ages: digit tops 0.67 rh under the tick line top
    ages = [(a, s) for a, s in ICS_AGE if a < t_old]
    order = [p for p in ages if p[0] in ICS_FIRST] + [p for p in ages if p[0] not in ICS_FIRST]
    spans, last = [], None
    for a, s in order:                               # ink centred on the tick, nudged by <= 0.6 mm or left out
        i0, _, i1, _ = ink_mm(s, fs)
        for dx in (0.0, -0.15, 0.15, -0.3, 0.3, -0.45, 0.45, -0.6, 0.6):
            lo, hi = X(a) + dx - 0.5 * (i1 - i0) - gap, X(a) + dx + 0.5 * (i1 - i0) + gap
            if all(hi <= l or lo >= h for l, h in spans):
                ink_text(ax, X(a) + dx, yb, s, fs); spans.append((lo, hi))
                break
        else:
            continue
        if a == 2.58:
            last = hi - gap
    if last is not None:                             # '(Ma)' ink starts 0.35 rh after the last age (terr)
        ink_text(ax, last + 0.35 * rh, yb, "(Ma)", fs, align="left")
    ya = yr + 1.5 * rh                               # numeric axis: majors every 50, minors every 25 Myr
    ax.plot([X(axis_from), X(0)], [ya, ya], color="black", lw=0.5, zorder=3, solid_capstyle="butt")
    yn = ya + 0.81 * rh + ink_mm("0", fs)[3]         # number tops 0.81 rh under the axis (terr)
    for a in np.arange(0.0, axis_from + 1e-9, 25.0):
        major = a % 50 == 0
        ax.plot([X(a)] * 2, [ya, ya + (0.30 if major else 0.15) * rh], color="black", lw=0.4, zorder=3,
                solid_capstyle="butt")
        if major:
            ink_text(ax, X(a), yn, "%d" % a, fs)
    if scale_x is not None:                          # bar on the axis line, '50 Myr' ink centred 0.50 rh under it
        x1 = scale_x + X(0) - X(50)
        ax.plot([scale_x, x1], [ya, ya], color="black", lw=0.5, zorder=3, solid_capstyle="butt")
        ink_text(ax, 0.5 * (scale_x + x1), ya + 0.50 * rh + ink_mm("50 Myr", fs)[3], "50 Myr", fs)
    return used, {"rh": rh, "yr": yr, "yb": yb, "ya": ya}


# ---- terr Extended Data Fig. 8 null histogram -----------------------------------------------------------------
def ed8_hist(ax, v, bins, at, lines, xlabel=None, ylabel="Count", log=False, head=1.45, fs=5, lab_fs=5.5,
             tick_fs=5, color=None):
    """terr ED Fig. 8a: blue bins with thin white edges, a full-height red line at the observed value with red text
    right-aligned just left of it at the top, horizontal grid only, black left and bottom axes, outward ticks.
    color overrides the bin colour (Fig. 4c v2: grey bins for the null panel). Returns the red text so the caller
    can check it stays inside the axes."""
    ax.hist(v, bins=bins, color=ED8_BAR if color is None else color, ec="white", lw=0.25, zorder=2)
    ax.axvline(at, color=ED8_RED, lw=0.8, zorder=3)
    t = ax.annotate("\n".join(lines), xy=(at, 0.97), xycoords=("data", "axes fraction"), xytext=(-2.5, 0),
                    textcoords="offset points", ha="right", va="top", fontsize=fs, color=ED8_RED, linespacing=1.15)
    if log:
        ax.set_xscale("log")
    ax.set_ylim(0, ax.get_ylim()[1] * head)
    ax.yaxis.grid(True, color=ED8_GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(direction="out", labelsize=tick_fs, length=2, width=0.5, pad=1.5)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=lab_fs, labelpad=1.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=lab_fs, labelpad=1.5)
    return t


# ---- Figs 2-4 v2 --------------------------------------------------------------------------------------------------
# Forms: grouped_bars + cond_box + sq_legend = Cell Fig. 4G/H; dot_whisker = Sackton et al. 2019 Fig. 4A; ref_line =
# Sackton Fig. 3C threshold line; null_hist = Sackton Fig. 3A; scatter_fit = Cell Fig. 4A; ed9_heat = terr Extended
# Data Fig. 9. Statistics: multinomial orthogroup bootstrap, Haldane 2x2 log OR, clade fixed-effects logistic. qa_v2.
import functools
import matplotlib.text as mtext
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.colors import LinearSegmentedColormap

PHYLO = ["Porifera", "Cnidaria", "Annelida", "Mollusca", "Arthropoda", "Echinodermata", "Chondrichthyes",
         "Actinopteri_m"]  # = Fig. 1 v2 FOCAL order (data name of Actinopteri)
ALL_LAB = "All clades"


def clab(c):
    return CLADE_LABEL.get(c, c)


def mt_arial():  # $P$, $K$, $n$ in Arial italic, embedded as TrueType like the rest of the text
    plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                         "mathtext.bf": "Arial:bold"})


def letter(fig, Wf, Hf, s, x_mm, y_mm):  # Fig. 1 v2 panel letter: 8 pt bold, baseline y_mm below the top edge
    return fig.text(x_mm / Wf, 1 - y_mm / Hf, s, fontsize=8, fontweight="bold", ha="left", va="baseline")


def expit(x):
    return 1 / (1 + np.exp(-x))


def logor2(e1, n1, e0, n0):  # 2x2 log OR deep vs null; +0.5 to every cell if any cell is 0 (= r6h lor, one clade)
    a, b, c, d = (np.asarray(v, float) for v in (e1, np.subtract(n1, e1), e0, np.subtract(n0, e0)))
    z = 0.5 * (np.minimum(np.minimum(a, b), np.minimum(c, d)) <= 0)
    return np.log((a + z) * (d + z) / ((b + z) * (c + z)))


def boot_w(rng, n, B=2000):  # resampling counts of n units drawn with replacement (= tabulated sample(n, n, TRUE))
    return rng.multinomial(n, np.full(n, 1.0 / n), size=B)


def pci(b, q=(.025, .975)):  # percentile CI along axis 0 (numpy 'linear' = R quantile type 7)
    return np.nanquantile(b, q[0], axis=0), np.nanquantile(b, q[1], axis=0)


def bp(b):  # two-sided bootstrap P of r6g/r6h: min(1, 2 min(mean(b <= 0), mean(b >= 0))), NA removed
    b = np.asarray(b, float); ok = ~np.isnan(b); k = ok.sum(0)
    return np.minimum(1, 2 * np.minimum((ok & (b <= 0)).sum(0) / k, (ok & (b >= 0)).sum(0) / k))


def pfmt(p, B=2000):  # exact P to two significant digits; a bootstrap P of 0 is shown as P < 0.001
    return "$P$ < 0.001" if p < 1.0 / B else ("$P$ = %.2g" % p if p < 0.1 else "$P$ = %.2f" % p)


def fe_logit(e1, n1, e0, n0, it=100, tol=1e-10):
    """Common log OR beta of glm(y ~ real + factor(clade), binomial) from grouped counts: logit p1 = a_c + beta,
    logit p0 = a_c. Arrays (..., C); returns beta (...). Newton steps with the clade intercepts profiled out (Schur
    complement); clades with no events or only events carry no information on beta and are dropped, as in glm."""
    e1, n1, e0, n0 = (np.asarray(v, float) for v in (e1, n1, e0, n0))
    ev, nn = e1 + e0, n1 + n0
    inf = (ev > 0) & (ev < nn)
    a = np.log((ev + .5) / (nn - ev + .5))
    b = np.zeros(e1.shape[:-1])
    for _ in range(it):
        p1, p0 = expit(a + b[..., None]), expit(a)
        h1, h0 = n1 * p1 * (1 - p1), n0 * p0 * (1 - p0)
        g1, g0 = e1 - n1 * p1, e0 - n0 * p0
        hs = np.where(inf, h1 + h0, 1.0)
        S = np.where(inf, h1 * h0 / hs, 0).sum(-1)
        gb = np.where(inf, g1 - h1 * (g1 + g0) / hs, 0).sum(-1)
        with np.errstate(divide="ignore", invalid="ignore"):
            db = np.clip(gb / S, -5, 5)
        da = np.where(inf, np.clip((g1 + g0 - h1 * db[..., None]) / hs, -5, 5), 0)
        a, b = a + da, b + db
        if np.nanmax(np.abs(db)) < tol:
            break
    return b


# ---- text QA on drawn artists (display px) -------------------------------------------------------------------------
def _quad(bb):
    return np.array([[bb.x0, bb.y0], [bb.x1, bb.y0], [bb.x1, bb.y1], [bb.x0, bb.y1]], float)


def _tbox(t, r):  # true (rotated) box of a drawn Text; a bbox patch (condition box) is the box
    if t.get_bbox_patch() is not None:
        return _quad(t.get_bbox_patch().get_window_extent(r))
    bb, rot = t.get_window_extent(r), t.get_rotation()
    if min(rot % 180, 180 - rot % 180) < .01:
        return _quad(bb)
    t.set_rotation(0); u = t.get_window_extent(r); t.set_rotation(rot)
    th = np.deg2rad(rot); R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    c = np.array([(bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2])
    return c + (np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * [u.width / 2, u.height / 2]) @ R.T


def _sat(P, Q, pad=0.0):  # convex quads overlap deeper than pad px (separating axis theorem)
    if (P[:, 0].max() <= Q[:, 0].min() + pad or Q[:, 0].max() <= P[:, 0].min() + pad or
            P[:, 1].max() <= Q[:, 1].min() + pad or Q[:, 1].max() <= P[:, 1].min() + pad):
        return False
    for M in (P, Q):
        for i in range(4):
            e = M[(i + 1) % 4] - M[i]; L = np.hypot(*e)
            if L > 0:
                nrm = np.array([-e[1], e[0]]) / L; p, q = P @ nrm, Q @ nrm
                if p.max() <= q.min() + pad or q.max() <= p.min() + pad:
                    return False
    return True


def _seg(u, v, rad):  # discs of radius rad every <= 1 px along segment u-v
    n = int(np.hypot(*(v - u))) + 2
    return np.column_stack([np.linspace(u[0], v[0], n), np.linspace(u[1], v[1], n), np.full(n, rad)])


def _hits(Q, ink, pad):  # text quad Q (shrunk by pad px) on bar rectangles or on point discs (x, y, radius)
    rects, pts = ink
    if any(_sat(Q, R, pad) for R in rects):
        return True
    e1, e2 = Q[1] - Q[0], Q[3] - Q[0]; L1, L2 = np.hypot(*e1), np.hypot(*e2)
    if not len(pts) or L1 == 0 or L2 == 0:
        return False
    d = pts[:, :2] - Q.mean(0)
    dx = np.maximum(np.abs(d @ (e1 / L1)) - (L1 / 2 - pad), 0)
    dy = np.maximum(np.abs(d @ (e2 / L2)) - (L2 / 2 - pad), 0)
    return bool(np.any(np.hypot(dx, dy) < pts[:, 2]))


def _ink(ax, r):  # data ink of an axes: bar rectangles; lines, caps, markers, error bars and scatter points as discs
    pxpt = ax.figure.dpi / 72.0
    rects = [_quad(p.get_window_extent(r)) for p in ax.patches if isinstance(p, Rectangle) and p.get_visible()]
    pts = []
    for ln in ax.lines:
        if not ln.get_visible():
            continue
        xy = ln.get_transform().transform(np.asarray(ln.get_xydata(), float))
        xy = xy[np.isfinite(xy).all(1)]
        if not len(xy):
            continue
        if ln.get_linestyle() not in ("None", "none", "", " ") and len(xy) > 1:
            pts += [_seg(xy[i], xy[i + 1], ln.get_linewidth() * pxpt / 2) for i in range(len(xy) - 1)]
        mk = ln.get_marker()
        if mk in ("_", "|"):  # error-bar caps: short segments, not discs
            h = ln.get_markersize() * pxpt / 2; w = max(ln.get_markeredgewidth(), .5) * pxpt / 2
            dv = np.array([h, 0.0]) if mk == "_" else np.array([0.0, h])
            pts += [_seg(p - dv, p + dv, w) for p in xy]
        elif mk not in (None, "None", "none", "", " "):
            pts.append(np.column_stack([xy, np.full(len(xy), ln.get_markersize() * pxpt / 2)]))
    for col in ax.collections:
        if not col.get_visible():
            continue
        if isinstance(col, LineCollection):
            lw = np.atleast_1d(col.get_linewidths())
            for k, s in enumerate(col.get_segments()):
                s = col.get_transform().transform(np.asarray(s, float))
                pts += [_seg(s[i], s[i + 1], lw[k % len(lw)] * pxpt / 2) for i in range(len(s) - 1)]
        elif isinstance(col, PathCollection):
            off = col.get_offset_transform().transform(np.asarray(col.get_offsets(), float))
            rad = np.sqrt(np.resize(np.atleast_1d(col.get_sizes()), len(off))) * pxpt / 2
            pts.append(np.column_stack([off, rad]))
    return rects, (np.vstack(pts) if pts else np.zeros((0, 3)))


def qa_v2(fig, min_fs=5.0, pad_pt=0.3, edge_mm=0.3):
    """Text QA on what is drawn: font size >= min_fs; every text inside the figure (and >= edge_mm from its edge);
    no two true (rotated) text boxes overlap (a legend, or a text with a bbox, is one box); no text or legend on the
    data ink of any axes. Texts flagged t._qa_skip (numbers inside heatmap cells) skip only the ink test."""
    seen, orig = [], mtext.Text.draw

    @functools.wraps(orig)
    def rec(self, renderer, *a, **k):
        if self.get_visible() and self.get_text().strip():
            seen.append(self)
        return orig(self, renderer, *a, **k)
    mtext.Text.draw = rec
    try:
        fig.canvas.draw()
    finally:
        mtext.Text.draw = orig
    r, pad = fig.canvas.get_renderer(), pad_pt * fig.dpi / 72.0
    T = list({id(t): t for t in seen}.values())
    legs = [l for l in [a.get_legend() for a in fig.axes] + list(fig.legends) if l is not None]
    inleg = {id(t) for l in legs for t in l.get_texts()}
    it = [(t.get_text(), _tbox(t, r), bool(getattr(t, "_qa_skip", False))) for t in T if id(t) not in inleg]
    it += [("<legend %s>" % "|".join(x.get_text() for x in l.get_texts()), _quad(l.get_window_extent(r)), False)
           for l in legs]
    fb, em = fig.bbox, edge_mm / 25.4 * fig.dpi
    res = {"texts": len(T), "legends": len(legs), "small": [], "outside": [], "near_edge": [], "overlap": [],
           "on_ink": []}
    res["small"] = ["%s (%.2f pt)" % (t.get_text(), t.get_fontsize()) for t in T if t.get_fontsize() < min_fs - 1e-6]
    for lab, Q, _ in it:
        lo, hi = Q.min(0), Q.max(0)
        if lo[0] < fb.x0 - .5 or lo[1] < fb.y0 - .5 or hi[0] > fb.x1 + .5 or hi[1] > fb.y1 + .5:
            res["outside"].append(lab)
        elif lo[0] < fb.x0 + em or lo[1] < fb.y0 + em or hi[0] > fb.x1 - em or hi[1] > fb.y1 - em:
            res["near_edge"].append(lab)
    for i in range(len(it)):
        for j in range(i + 1, len(it)):
            if _sat(it[i][1], it[j][1], pad):
                res["overlap"].append((it[i][0], it[j][0]))
    inks = [(k, _ink(a, r)) for k, a in enumerate(fig.axes)]
    for lab, Q, sk in it:
        hit = [k for k, ink in inks if not sk and _hits(Q, ink, pad)]
        if hit:
            res["on_ink"].append((lab, hit))
    return res


def qa_report(res, name):  # ASCII counts on the console, full lists (utf-8) in out/<name>_qa_v2.txt; True if clean
    print(name, "QA", " ".join("%s=%s" % (k, len(v) if isinstance(v, list) else v) for k, v in res.items()))
    with open(os.path.join(V2, name + "_qa_v2.txt"), "w", encoding="utf-8") as f:
        f.write("".join("%s\t%s\n" % (k, x) for k, v in res.items() if isinstance(v, list) for x in v))
    return not any(isinstance(v, list) and v for v in res.values())


# ---- panel forms ---------------------------------------------------------------------------------------------------
def grouped_bars(ax, labels, vals, lo, hi, colors, gap=0.6, extra=(), extra_gap=0.8, cap_frac=0.36, lab_fs=5.5):
    """Cell Fig. 4G/H: S adjacent edge-less bars per group, groups gap bar widths apart, capped error bars (0.5 pt),
    centred horizontal group labels without x ticks, full-width bottom spine. vals/lo/hi (G, S). Group centres out."""
    vals, lo, hi = (np.asarray(v, float) for v in (vals, lo, hi))
    G, S = vals.shape; xs, x = [], 0.0
    for g in range(G):
        x += extra_gap if g in extra else 0.0
        xs.append(x + np.arange(S) + .5); x += S + gap
    xs = np.array(xs); ax.set_xlim(-.6, x - gap + .6)
    ppd = ax.get_position().width * ax.figure.get_figwidth() * 72 / (x - gap + 1.2)  # pt per bar width
    for s in range(S):
        ax.bar(xs[:, s], vals[:, s], width=1.0, color=colors[s], linewidth=0, zorder=2)
        ax.errorbar(xs[:, s], vals[:, s], yerr=[vals[:, s] - lo[:, s], hi[:, s] - vals[:, s]], fmt="none",
                    ecolor=INK, elinewidth=.5, capsize=cap_frac / 2 * ppd, capthick=.5, zorder=3)
    ax.set_xticks(xs.mean(1)); ax.set_xticklabels(labels, fontsize=lab_fs)
    ax.tick_params(axis="x", length=0, pad=2.0)
    ax.tick_params(axis="y", direction="out", labelsize=5, length=2, width=.5, pad=1.5)
    return xs.mean(1)


def cond_box(ax, s, x=1.0, y=1.0, ha="right", va="top", fs=5.5, pad=.35):  # Cell Fig. 4G: dark box, white text
    return ax.annotate(s, xy=(x, y), xycoords="axes fraction", fontsize=fs, color="white", ha=ha, va=va,
                       xytext=((-1 if ha == "right" else 1) * pad * fs, (-1 if va == "top" else 1) * pad * fs),
                       textcoords="offset points", bbox=dict(boxstyle="square,pad=%g" % pad, fc=INK, ec="none"))


def sq_legend(fig, x_frac, y_frac, labels, colors, fs=5.5):  # Cell Fig. 4H: square patches, no frame, outside right
    from matplotlib.patches import Patch
    return fig.legend(handles=[Patch(color=c, lw=0, label=l) for l, c in zip(labels, colors)], loc="center left",
                      bbox_to_anchor=(x_frac, y_frac), bbox_transform=fig.transFigure, fontsize=fs, frameon=False,
                      handlelength=1.25, handleheight=.9, handletextpad=.6, labelspacing=.45, borderpad=0,
                      borderaxespad=0)


def dot_whisker(ax, xs, est, lo, hi, labels, color=INK, ecolor=ARROW, ms=2.8, cap=2.1, rot=45, lab_fs=5.5):
    """Sackton et al. 2019 Fig. 4A: filled dots, grey capped whiskers (0.5 pt; cap ~1.5 dot widths), labels at 45
    degrees anchored at their right end, no x ticks."""
    est, lo, hi = (np.asarray(v, float) for v in (est, lo, hi))
    dots = ax.errorbar(xs, est, yerr=[est - lo, hi - est], fmt="o", ms=ms, mfc=color, mec=color, ecolor=ecolor,
                       elinewidth=.5, capsize=cap, capthick=.5, zorder=3)[0]
    dots.set_markeredgewidth(0)  # mew= in errorbar() would also zero the caps (it overrides capthick)
    ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=rot, ha="right", rotation_mode="anchor", fontsize=lab_fs)
    ax.tick_params(axis="x", length=0, pad=2.0)
    ax.tick_params(axis="y", direction="out", labelsize=5, length=2, width=.5, pad=1.5)
    ax.set_xlim(xs[0] - .7, xs[-1] + .7)


def ref_line(ax, y, label=None, color=NULL, x=0.015, fs=5):  # Sackton Fig. 3C: dashed grey line, label above at left
    ax.axhline(y, color=color, lw=.5, ls=(0, (3, 2)), zorder=1)
    if label:
        ax.annotate(label, xy=(x, y), xycoords=("axes fraction", "data"), xytext=(0, .8), textcoords="offset points",
                    ha="left", va="bottom", fontsize=fs, color=color)


# ---- terr Extended Data Fig. 9 annotated heatmap ------------------------------------------------------------------
# colour bar of terr ED Fig. 9 sampled from the published JPEG (scratch z_ed9_cbar), from the bottom (low, red) to
# the top (high, blue); white is the middle of the scale
ED9_STOPS = ["#bc2a36", "#ce5449", "#dc7c66", "#e9a486", "#eecab2", "#f3e8dd", "#f5f5f3", "#e1f0f5", "#c7e3ee",
             "#a3cfe2", "#79b5ce", "#5899c3", "#487db2"]
ED9_CMAP = LinearSegmentedColormap.from_list("ed9", ED9_STOPS)
ED9_GRID = "#ebebeb"


def ed9_axes(fig, Wf, Hf, x0, y0, nr, nc, tw, th, stub=0.1):
    """Axes for ed9_heat: nr x nc abutting tiles of tw x th mm whose block has its top-left corner at (x0, y0) mm.
    The axes reach stub tile widths past the block on every side for terr's grid-line stubs. Returns the axes and
    the stub (x, y) in tile units."""
    s = stub * tw
    ax = axes_mm(fig, Wf, Hf, x0 - s, y0 - s, nc * tw + 2 * s, nr * th + 2 * s)
    return ax, (stub, s / th)


def ed9_heat(ax, fill, text, rows, cols, vmax, stub, fs=5.5, lab_fs=5.5, rot=35, cmap=ED9_CMAP):
    """terr ED Fig. 9: abutting tiles (no gaps) coloured by fill on [-vmax, vmax] (white = 0), black numbers centred
    in the tiles, light-grey grid lines through the row and column centres that show only as short stubs past the
    block on all four sides, row labels right-aligned at the left, column labels at rot degrees ending under the tile
    centres, no spines or ticks. The numbers are flagged _qa_skip (they sit on the tile colour by design)."""
    fill = np.asarray(fill, float); nr, nc = fill.shape; sx, sy = stub
    ax.pcolormesh(np.arange(nc + 1) - .5, np.arange(nr + 1) - .5, fill, cmap=cmap, vmin=-vmax, vmax=vmax,
                  edgecolors="face", linewidth=.2, zorder=2)
    for i in range(nr):
        ax.plot([-.5 - sx, nc - .5 + sx], [i, i], color=ED9_GRID, lw=.7, zorder=1, solid_capstyle="butt")
    for j in range(nc):
        ax.plot([j, j], [-.5 - sy, nr - .5 + sy], color=ED9_GRID, lw=.7, zorder=1, solid_capstyle="butt")
    for i in range(nr):
        for j in range(nc):
            t = ax.text(j, i, text[i][j], ha="center", va="center_baseline", fontsize=fs, color="black", zorder=3)
            t._qa_skip = True
    ax.set_xlim(-.5 - sx, nc - .5 + sx); ax.set_ylim(nr - .5 + sy, -.5 - sy)
    ax.set_yticks(range(nr)); ax.set_yticklabels(rows, fontsize=lab_fs)
    ax.set_xticks(range(nc))
    ax.set_xticklabels(cols, rotation=rot, ha="right", rotation_mode="anchor", fontsize=lab_fs)
    ax.tick_params(length=0, pad=1.5)
    for s in ax.spines.values():
        s.set_visible(False)


def ed9_cbar(cax, vmax, ticks, labels, title, fs=5, title_fs=5.5, cmap=ED9_CMAP):
    """terr ED Fig. 9 colour bar: vertical, flat ends, no outline, short white tick dashes inside both long edges,
    tick labels to the right, title left-aligned above the bar."""
    v = np.linspace(-vmax, vmax, 257)
    cax.pcolormesh([0, 1], v, ((v[:-1] + v[1:]) / 2)[:, None], cmap=cmap, vmin=-vmax, vmax=vmax, edgecolors="face",
                   linewidth=.2, zorder=1)
    cax.set_xlim(0, 1); cax.set_ylim(-vmax, vmax)
    for y in ticks:
        for xa, xb in ((0, .22), (.78, 1)):
            cax.plot([xa, xb], [y, y], color="white", lw=.5, zorder=2, solid_capstyle="butt")
    cax.set_xticks([]); cax.yaxis.tick_right(); cax.set_yticks(ticks); cax.set_yticklabels(labels, fontsize=fs)
    cax.tick_params(axis="y", length=0, pad=1.2)
    for s in cax.spines.values():
        s.set_visible(False)
    return cax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 2.5), textcoords="offset points",
                        ha="left", va="bottom", fontsize=title_fs)


# ---- Cell Fig. 4A framed scatter -----------------------------------------------------------------------------------
def cell_frame(ax, lw=.8):  # heavy frame on all four sides, outward ticks on the left and bottom only
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(lw); s.set_color(INK)
    ax.tick_params(direction="out", top=False, right=False, length=2, width=.6, labelsize=5, pad=1.5)


def scatter_fit(ax, x, y, color, line_color, fit, s=1.7, alpha=.6, band_alpha=.25, z=2):
    """Cell Fig. 4A, one series: edge-less dots; an OLS line in a darker colour over the fitted points' x range and a
    translucent 95% confidence band of the mean, both drawn over all dots. fit = mask of the points used for the line.
    Returns n, Pearson r, slope and intercept of the fitted points."""
    from scipy import stats
    x, y, fit = np.asarray(x, float), np.asarray(y, float), np.asarray(fit, bool)
    ax.scatter(x, y, s=s, c=color, alpha=alpha, linewidths=0, zorder=z)
    xf, yf = x[fit], y[fit]; n = len(xf)
    b1, b0 = np.polyfit(xf, yf, 1); r = np.corrcoef(xf, yf)[0, 1]
    g = np.linspace(xf.min(), xf.max(), 200); e = yf - (b0 + b1 * xf)
    se = np.sqrt(e @ e / (n - 2) * (1 / n + (g - xf.mean()) ** 2 / ((xf - xf.mean()) ** 2).sum()))
    tq = stats.t.ppf(.975, n - 2)
    ax.fill_between(g, b0 + b1 * g - tq * se, b0 + b1 * g + tq * se, color=line_color, alpha=band_alpha, lw=0,
                    zorder=z + 10)
    ax.plot(g, b0 + b1 * g, color=line_color, lw=.8, zorder=z + 20, solid_capstyle="butt")
    return n, r, b1, b0


def cell_legend(ax, labels, colors, fs=5, loc="upper left"):  # Cell Fig. 4A: circle markers, no frame, top left
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls="none", marker="o", ms=3, mfc=c, mec="none", label=l) for l, c in zip(labels, colors)]
    return ax.legend(handles=h, loc=loc, fontsize=fs, frameon=False, handlelength=.8, handletextpad=.3,
                     labelspacing=.3, borderpad=.3, borderaxespad=.2)
#@@APPEND@@
