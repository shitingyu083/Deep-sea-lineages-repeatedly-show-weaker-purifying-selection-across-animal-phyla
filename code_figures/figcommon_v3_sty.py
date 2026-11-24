"""Shared helpers for the v3 figures (Figs 2-4): every panel copies the form of one viewed published panel.

Benchmarks (crops in stage_sty/_scratch_sty/ref_crops, originals in _ref_sty; panel map in
paper/figure_benchmark_map.md, section "Fig. 2-4 v3"):
  berg 1c - Bergeron et al. 2023 Nature, doi:10.1038/s41586-023-05752-y, Fig. 1c: time-scaled tree, upright names,
            dotted leaders, silhouettes, black bars from the reference value with grey capped 95% intervals,
            half-step grey stripes, thin reference line, two-line header with one double-headed arrow
  berg 1a - same paper, Fig. 1a: time tree with age axis, one dot per estimate, black-edged mean square with capped
            interval, alternating stripes, count boxes
  berg 2  - same paper, Fig. 2a,b: scatter with large edgeless dots, dashed regression line, centred statistics text
  cell 2B - Xu et al. 2025 Cell, doi:10.1016/j.cell.2025.01.002, Fig. 2B: tree, two framed stacked-bar panels with
            stacked box legends
  terr    - Wei et al. 2026 Nature, doi:10.1038/s41586-025-09722-4, Extended Data histogram (terr_ED8.jpg) and
            annotated heatmap (terr_ED9.jpg), drawn with the v2 helpers ed8_hist / ed9_heat
Everything of v2 (palette, statistics, QA) is imported unchanged; save_v2 / write_facts / qa_report write to
viz_v3_sty/out (their module global V2 is redirected here), never to the package.
"""
import os, sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "viz_v2_sty"))
import figcommon_v2_sty as _v2
from figcommon_v2_sty import *  # noqa: F401,F403
from PIL import Image

V3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(V3, exist_ok=True)
_v2.V2 = V3

STRIPE = "#efefef"   # berg 1c background stripes
CI = ARROW           # berg 1c grey interval lines

# per-clade categorical colour (berg 1a/1c: every clade owns one hue, carried through the silhouette, the dotted
# leader, the bar or estimate mark and the count box). Eight hues walked round the circle so that neighbouring rows
# never share a hue; all are mid-lightness and print- and colour-vision-safe at 5-7 pt. All clades stays near-black.
CLADE_COL = {"Porifera": "#8C6239", "Cnidaria": "#D6603F", "Annelida": "#A277B8", "Mollusca": "#E3A93B",
             "Arthropoda": "#5C9B4A", "Echinodermata": "#D0578F", "Chondrichthyes": "#4C7EA8",
             "Actinopteri": "#2E9E96", "All clades": "#2B2B2B"}
assert set(CLADE_COL) == {clab(c) for c in PHYLO} | {ALL_LAB}, sorted(CLADE_COL)
CLASS_FOCAL = ("Chondrichthyes", "Actinopteri")
PP = os.path.join(ROOT, "figures_src_sty", "phylopic")
PNG = {f.split("__")[0]: os.path.join(PP, f) for f in os.listdir(PP) if f.endswith(".png")}
assert set(PNG) == {clab(c) for c in PHYLO}, sorted(PNG)


# ---- time tree pruned to the focal crowns (berg 1a/1c, cell 2B: ultrametric, tips aligned at the present) ------
class PN:
    def __init__(self, age, kids=(), clade=None):
        self.age, self.kids, self.clade, self.y = age, list(kids), clade, None


def leaves(p):
    return [p.clade] if not p.kids else [c for k in p.kids for c in leaves(k)]


def focal_tree(clades):
    """LSD time tree (final_v3/date/lsd.timetree.nex) pruned to the crown nodes of the given data clades: a crown is
    a leaf drawn to the present, empty branches are dropped and single-child nodes pass through. Children are
    rotated (topology unchanged) so that the leaves come in the order of `clades`."""
    tree = parse_nex_tree(W("final_v3", "date", "lsd.timetree.nex"))
    lab = {r["code"]: r for r in rtsv(os.path.join(DATA, "meta", "labels_v2.tsv"))}
    for n in postorder(tree):
        if not n.children:
            n.age = 0.0 if n.age is None else n.age
        assert n.age is not None, "undated internal node"
    cr = {}
    for c in clades:
        nm = clab(c)
        ts = frozenset(t for t in tree.tips if lab[t]["class" if nm in CLASS_FOCAL else "phylum"] == nm)
        n = mrca(tree, ts)
        assert n.tips == ts, nm   # the focal taxon is monophyletic in the tree
        cr[id(n)] = c

    def prune(n):
        if id(n) in cr:
            return PN(n.age, (), cr[id(n)])
        ks = [k for k in (prune(c) for c in n.children) if k is not None]
        return None if not ks else ks[0] if len(ks) == 1 else PN(n.age, ks)

    def rot(p):
        for k in p.kids:
            rot(k)
        p.kids.sort(key=lambda k: min(clades.index(c) for c in leaves(k)))

    root = prune(tree)
    rot(root)
    assert leaves(root) == list(clades), leaves(root)
    return root


def draw_tree(ax, root, x0, x1, ymap, age0=None, lw=0.5, color="k"):
    """Rectangular time tree, root (or age0, Ma) at x0, present at x1, no root stem; leaf y from ymap.
    Returns the age -> x function."""
    A = root.age if age0 is None else age0
    X = lambda a: x0 + (A - a) / A * (x1 - x0)

    def sety(p):
        if not p.kids:
            p.y = ymap[p.clade]
        else:
            for k in p.kids:
                sety(k)
            p.y = 0.5 * (p.kids[0].y + p.kids[-1].y)

    def dr(p):
        xp = X(p.age)
        ax.plot([xp, xp], [p.kids[0].y, p.kids[-1].y], color=color, lw=lw, solid_capstyle="projecting", zorder=4)
        for k in p.kids:
            ax.plot([xp, x1 if not k.kids else X(k.age)], [k.y, k.y], color=color, lw=lw, solid_capstyle="butt",
                    zorder=4)
            if k.kids:
                dr(k)

    sety(root)
    dr(root)
    return X


def nodes(p):
    return [p] + [n for k in p.kids for n in nodes(k)]


# ---- berg 1c data column --------------------------------------------------------------------------------------
MINUS = "−"   # true minus sign in hand-drawn tick labels


def berg_col(ax, x0, lo, hi, unit, rows, yaxis, yh1, yh2, ref=0.0, head1=None, head2=None, arrow_to=0.8,
             xlabel=None, ticks=None, tfs=5.0, fs=5.5, lab_fs=6.0, spine_top=None):
    """berg 1c data column. x0, spine at value lo; stripes from rows[0] (top) to rows[1] (bottom) on half a tick
    step from each tick; thin black reference line at `ref` under header line 1; two-line header with one
    double-headed arrow from the spine through the reference line to `arrow_to`. Returns the value -> x function."""
    u = lambda v: x0 + (v - lo) / (hi - lo) * (hi - lo) * unit
    ticks = list(range(int(np.ceil(lo)), int(np.floor(hi)) + 1)) if ticks is None else list(ticks)
    step = ticks[1] - ticks[0]
    ytop, ybot = rows
    for t in ticks:
        ax.fill_betweenx([ybot, ytop], u(t), u(min(t + step / 2, hi)), color=STRIPE, lw=0, zorder=0.3)
    if head1 is not None:
        ax.plot([u(ref), u(ref)], [yaxis, yh1 + 0.6], color="k", lw=0.5, solid_capstyle="butt", zorder=1.5)
        ax.text(u(ref), yh1, head1, ha="center", va="baseline", fontsize=fs, color=INK, zorder=10)
    ax.plot([u(lo), u(lo)], [yaxis, yh2 - 0.6 if spine_top is None else spine_top], color="k", lw=0.7,
            solid_capstyle="butt", zorder=1.6)   # r21: spine_top lets a caller stop the spine below its own header
    if head2 is not None:
        ax.annotate("", xy=(u(arrow_to), yh2), xytext=(u(lo), yh2), annotation_clip=False,
                    arrowprops=dict(arrowstyle="<|-|>,head_length=0.4,head_width=0.2", lw=0.5, color="k",
                                    shrinkA=0, shrinkB=0, mutation_scale=4), zorder=10)
        ax.text(u(lo) - 0.8, yh2, head2[0], ha="right", va="center", fontsize=fs, color=INK, zorder=10)
        ax.text(u(arrow_to) + 0.8, yh2, head2[1], ha="left", va="center", fontsize=fs, color=INK, zorder=10)
    ax.plot([u(lo), u(hi)], [yaxis, yaxis], color="k", lw=0.6, solid_capstyle="butt", zorder=1.6)
    for t in ticks:
        ax.plot([u(t), u(t)], [yaxis, yaxis + 1.0], color="k", lw=0.5, solid_capstyle="butt", zorder=1.6)
        ax.text(u(t), yaxis + 1.4, (MINUS + "%g" % -t) if t < 0 else "%g" % t, ha="center", va="top", fontsize=tfs,
                color=INK, zorder=10)
    if xlabel is not None:
        ax.text(0.5 * (u(lo) + u(hi)), yaxis + 1.4 + tfs * 0.46, xlabel, ha="center", va="top", fontsize=lab_fs,
                color=INK, zorder=10)
    return u


def berg_bars(ax, u, ys, est, lo95, hi95, base=0.0, pitch=6.0, bar=0.45, cap=0.6, color="k", ci_color=None):
    """berg 1c bars: solid bars from `base` (origin when base == lo of the axis), grey capped intervals on top."""
    h = bar * pitch
    cc = CI if ci_color is None else ci_color
    for y, e, l, hh in zip(ys, est, lo95, hi95):
        xb, xe = u(base), u(e)
        ax.add_patch(Rectangle((min(xb, xe), y - h / 2), abs(xe - xb), h, fc=color, ec="none", lw=0, zorder=2))
        if l is not None and hh is not None:
            ax.plot([u(l), u(hh)], [y, y], color=cc, lw=0.6, solid_capstyle="butt", zorder=3)
            for xx in (u(l), u(hh)):
                ax.plot([xx, xx], [y - cap * h / 2, y + cap * h / 2], color=cc, lw=0.6, solid_capstyle="butt",
                        zorder=3)


def name_row(ax, x, y, s, xspine, png=None, h=None, fs=6.0, color=INK):
    """berg 1c row label: upright name, dotted leader, silhouette right-aligned before the spine. A row without a
    silhouette (All clades) gets the leader alone, ending 1 mm before the spine."""
    t = ax.text(x, y, s, ha="left", va="center", fontsize=fs, color=color, zorder=10)
    x1 = x + text_mm(ax, t)[0]
    if png is None:
        ax.plot([x1 + 0.8, xspine - 1.0], [y, y], color=LEADER, lw=0.6, ls=(0, (1, 1.2)), solid_capstyle="butt",
                zorder=5)
        return x1
    sw = h * aspect(png)
    x0 = xspine - 1.0 - sw
    ax.plot([x1 + 0.8, x0 - 0.8], [y, y], color=LEADER, lw=0.6, ls=(0, (1, 1.2)), solid_capstyle="butt", zorder=5)
    silhouette(ax, png, x0, y, h=h, z=6)
    return x1


def aspect(path):
    a = np.asarray(Image.open(path).convert("LA"))[:, :, 1]
    ys, xs = np.nonzero(a > 5)
    return (xs.max() - xs.min() + 1) / (ys.max() - ys.min() + 1)


# ---- per-clade colour layer (berg 1a/1c + cell 2B) -------------------------------------------------------------
def ccol(c):
    """Categorical hue of a clade: the data name, its display label or ALL_LAB."""
    return CLADE_COL[clab(c)]


def silhouette_c(ax, path, x, y, color, w=None, h=None, ha="left", z=6, alpha=1.0):
    """v2 silhouette tinted to `color` (berg 1c colours its silhouettes by group). r22: vector outline via
    figcommon_v2_sty.vector_mask when available (no embedded image)."""
    if hasattr(_v2, "vector_mask"):
        return _v2.vector_mask(ax, path, x, y, color, w=w, h=h, ha=ha, z=z, alpha=alpha)
    a = np.asarray(Image.open(path).convert("LA"))[:, :, 1].astype(float) / 255.0
    ys, xs = np.nonzero(a > 0.02)
    a = a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    ar = a.shape[1] / a.shape[0]
    w = h * ar if w is None else w; h = w / ar
    rgba = np.zeros(a.shape + (4,)); rgba[:, :, 3] = a * alpha
    rgba[:, :, :3] = np.asarray(to_rgb(color))[None, None, :]
    x0 = x - w if ha == "right" else x - w / 2 if ha == "center" else x
    ax.imshow(rgba, extent=(x0, x0 + w, y + h / 2, y - h / 2), origin="upper", aspect="auto",
              interpolation="antialiased", zorder=z)
    return x0, x0 + w, h


def name_row_c(ax, x, y, s, xspine, clade, png=None, h=None, fs=6.0, color=INK, lead=True):
    """berg 1c row label in the clade hue: name in `color` (ink by default, as in the benchmark), the dotted leader
    and the silhouette both in the clade colour. A row without a silhouette gets the leader alone."""
    cc = ccol(clade)
    t = ax.text(x, y, s, ha="left", va="center", fontsize=fs, color=color, zorder=10)
    x1 = x + text_mm(ax, t)[0]
    if png is None:
        if lead:
            ax.plot([x1 + 0.8, xspine - 1.0], [y, y], color=cc, lw=0.6, ls=(0, (1, 1.2)),
                    solid_capstyle="butt", zorder=5)
        return x1
    sw = h * aspect(png)
    x0 = xspine - 1.0 - sw
    if lead:
        ax.plot([x1 + 0.8, x0 - 0.8], [y, y], color=cc, lw=0.6, ls=(0, (1, 1.2)), solid_capstyle="butt", zorder=5)
    silhouette_c(ax, png, x0, y, cc, h=h, z=6)
    return x1


def berg_bars_c(ax, u, ys, est, lo95, hi95, colors, base=0.0, pitch=6.0, bar=0.45, cap=0.6, ci_color=None):
    """berg_bars with one fill colour per row (a scalar fills every row, as in v2); intervals stay grey."""
    cs = colors if not isinstance(colors, str) else [colors] * len(ys)
    h = bar * pitch
    cc = CI if ci_color is None else ci_color
    for y, e, l, hh, fc in zip(ys, est, lo95, hi95, cs):
        xb, xe = u(base), u(e)
        ax.add_patch(Rectangle((min(xb, xe), y - h / 2), abs(xe - xb), h, fc=fc, ec="none", lw=0, zorder=2))
        if l is not None and hh is not None:
            ax.plot([u(l), u(hh)], [y, y], color=cc, lw=0.6, solid_capstyle="butt", zorder=3)
            for xx in (u(l), u(hh)):
                ax.plot([xx, xx], [y - cap * h / 2, y + cap * h / 2], color=cc, lw=0.6, solid_capstyle="butt",
                        zorder=3)


def chip_row(ax, x, y, s, clade, w=1.9, h=1.9, gap=1.0, fs=6.0, color=INK):
    """cell 2B row label: a filled colour chip immediately left of the name, one hue per row."""
    ax.add_patch(Rectangle((x, y - h / 2), w, h, fc=ccol(clade), ec="none", lw=0, zorder=6))
    return ax.text(x + w + gap, y, s, ha="left", va="center", fontsize=fs, color=color, zorder=10)


def end_val(ax, u, y, txt, est, lo=None, hi=None, base=0.0, fs=5.0, color=INK, gap=0.7):
    """`txt` just outside the ink of one berg row (bar end or outer CI cap), on the side away from `base`."""
    right = est >= base
    tip = max([est] + ([hi] if (hi is not None and right) else [])) if right else         min([est] + ([lo] if (lo is not None and not right) else []))
    return ax.text(u(tip) + (gap if right else -gap), y, txt, ha="left" if right else "right", va="center",
                   fontsize=fs, color=color, zorder=10)


# ---- cell 2B framed stacked-bar panel -------------------------------------------------------------------------
def lum(hexcol):
    """Rec.709 relative luminance of a hex colour, used to pick white or ink text on a filled box."""
    r, g, b = (int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def val_box(ax, xr, y, txt, color, h=3.1, pad=0.7, fs=5.0, fill=True):
    """berg 1a count box, right-aligned to `xr`: a box in the row colour holding one value, its text white on a dark
    fill and ink on a light one. `fill=False` draws an open box (white fill, coloured edge and ink text)."""
    t = ax.text(0, y, txt, ha="right", va="center", fontsize=fs, color=INK, zorder=12)
    w = text_mm(ax, t)[0] + 2 * pad
    t.remove()
    x0 = xr - w
    if fill:
        ax.add_patch(Rectangle((x0, y - h / 2), w, h, fc=color, ec="none", lw=0, zorder=10))
        fc = "w" if lum(color) < 0.45 else INK
    else:
        ax.add_patch(Rectangle((x0, y - h / 2), w, h, fc="white", ec=color, lw=0.5, zorder=10))
        fc = INK
    return ax.text(x0 + w / 2, y, txt, ha="center", va="center", fontsize=fs, color=fc, zorder=11)


def cell2b_frame(fig, Wf, Hf, x0, y0, w, h, xlim, xticks, title, tfs=5.0):
    """Cell 2B frame: thin grey full rectangle, outward ticks and labels under the bottom spine, no y ticks,
    title centred above the frame. Returns the axes; (x0, y0) is the top-left corner in mm."""
    ax = axes_mm(fig, Wf, Hf, x0, y0, w, h)
    ax.set_xlim(*xlim)
    ax.set_ylim(h, 0)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color(NULL)
        s.set_linewidth(0.5)
    ax.set_yticks([])
    ax.set_xticks(list(xticks))
    ax.tick_params(axis="x", direction="out", length=2, width=0.5, pad=1.2, labelsize=tfs, colors=INK)
    ax.set_title(title, fontsize=7, color=INK, pad=2.0)
    ax.set_facecolor("none")
    return ax
