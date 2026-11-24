"""Fig. 1 v2 in the form of terr Fig. 1 (Wei et al., Nature 649, 638-646, 2026; doi:10.1038/s41586-025-09722-4).
a, time tree with collapsed single-state clades, phylum bands, silhouettes with arrows to node boxes, age bars at the
focal crowns, higher-taxon names with grey leaders and the ICS time bar (terr Fig. 1). b, stem ages of the disjoint
deep origins in the form of terr Extended Data Fig. 8a (histogram with the median as a red line and red text). Same
data and numbers as fig1_sty.py, asserted against build/fig1_facts.tsv. Output: viz_v2_sty/out/Fig1.{pdf,png,tif}
and fig1_facts_v2.tsv (never the package)."""
from collections import Counter, defaultdict
from figcommon_v2_sty import *  # noqa: F401,F403
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

FW, FH = 180.0, 162.0                                # figure (mm)
Y0, Y1 = 6.0, 147.0                                  # tree rows (mm, y downwards)
XTIP, XB1, XLAB, XBOX, WBOX = 90.0, 115.0, 114.5, 118.0, 17.0
LW, MINB, GAP = 0.4, 1.9, 0.3                        # branch width (pt), minimum band (mm), gap between triangles (rows)
FOCAL = ["Porifera", "Cnidaria", "Annelida", "Mollusca", "Arthropoda", "Echinodermata", "Chondrichthyes",
         "Actinopteri"]
CLASS_FOCAL = ("Chondrichthyes", "Actinopteri")
NA = {"", "NA", "na", "None", "nan", "-"}


def XT(a):                                           # Ma -> mm: 850 Ma at 15.5 mm, present at 90 mm (= XTIP)
    return 15.5 + (850.0 - a) * 74.5 / 850.0


tree = parse_nex_tree(W("final_v3", "date", "lsd.timetree.nex"))
lab = {r["code"]: r for r in rtsv(os.path.join(DATA, "meta", "labels_v2.tsv"))}
NODES = postorder(tree)
for n in NODES:
    if not n.children:
        n.age = 0.0 if n.age is None else n.age
    assert n.age is not None, "undated internal node"
    n.children = n.children[::-1]                    # Porifera on top, Chordata at the bottom (terr order)


def phy(t):
    return lab[t]["phylum"]


def st(n):
    s = {lab[t]["lab"] for t in n.tips}
    return s.pop() if len(s) == 1 else "M"


def has(n, p):
    return any(phy(t) == p for t in n.tips)


prot = mrca(tree, [t for t in tree.tips if phy(t) in ("Mollusca", "Arthropoda")])
prot.children.sort(key=lambda c: 0 if has(c, "Mollusca") else 1)
assert len(tree.tips) == 426 and {lab[t]["lab"] for t in tree.tips} <= {"D", "S", "U"}
rows = []


def walk(n):                                         # a row = tip, or maximal single-state single-phylum clade
    if not n.children or (st(n) in ("D", "S") and len({phy(t) for t in n.tips}) == 1):
        rows.append(n)
        return
    for c in n.children:
        walk(c)


walk(tree)
assert sum(len(r.tips) for r in rows) == 426
dsets = [frozenset(r["tips"].split(",")) for r in rtsv(W("final_v3", "asr", "deep_origins.tsv"))]
assert len(dsets) == 97 and len(set(dsets)) == 97
assert set().union(*[r.tips for r in rows if st(r) == "D"]) == set().union(*dsets)
blocks = []                                          # contiguous phylum runs: [phylum, first row, last row]
for i, r in enumerate(rows):
    p = phy(next(iter(r.tips)))
    if blocks and blocks[-1][0] == p:
        blocks[-1][2] = i
    else:
        blocks.append([p, i, i])
assert len({b[0] for b in blocks}) == len(blocks), [b[0] for b in blocks]
units = [max(1.0, 0.6 * len(r.tips)) for r in rows]
bu = [sum(units[b[1]:b[2] + 1]) for b in blocks]
ann = mrca(tree, [t for t in tree.tips if phy(t) == "Annelida"])
first = [min(i for i, r in enumerate(rows) if r.tips <= c.tips) for c in ann.children]
assert len(first) == 2 and first[0] < first[1], first
FORK = {first[1]: 1.2}                               # space at the Annelida crown fork for its age label (mm)
gb = [sum(FORK.get(i, 0.0) for i in range(b[1], b[2] + 1)) for b in blocks]
lo, hi = 0.01, 5.0
for _ in range(60):                                  # unit height so that the phylum bands fill Y0..Y1
    s = 0.5 * (lo + hi)
    lo, hi = (s, hi) if sum(max(MINB, s * u) + g for u, g in zip(bu, gb)) < Y1 - Y0 else (lo, s)
ext, yb, y = [0.0] * len(rows), [], Y0
for b, u, g in zip(blocks, bu, gb):
    h = max(MINB, lo * u)
    yb.append((y, y + h + g))
    for i in range(b[1], b[2] + 1):
        y += FORK.get(i, 0.0)
        ext[i] = units[i] * h / u
        rows[i].y = y + 0.5 * ext[i]
        y += ext[i]
rowid = {id(r) for r in rows}
inrow = {id(m) for r in rows for m in postorder(r) if m is not r}
for n in NODES:
    if id(n) in inrow:
        continue
    n.x = XT(n.age)
    if n.children and id(n) not in rowid:
        n.y = 0.5 * (n.children[0].y + n.children[-1].y)
print("rows %d blocks %d unit %.3f mm" % (len(rows), len(blocks), lo))
print(" ".join("%s:%d" % (b[0][:7], b[2] - b[1] + 1) for b in blocks))
FOC = {cl: frozenset(t for t in tree.tips if lab[t]["clade"] == cl) for cl in FOCAL}
CNT = {cl: Counter(lab[t]["lab"] for t in FOC[cl]) for cl in FOCAL}
NCR = {cl: sum(1 for d in dsets if d <= FOC[cl]) for cl in FOCAL}
ORG = [r for r in rtsv(W("addon", "A2", "origins.tsv")) if r["set"] == "disjoint"]
NOR = Counter(r["clade"] for r in ORG)
assert sum(NCR.values()) == 97 and len(ORG) == 65 and sum(NOR[cl] for cl in FOCAL) == 65


def rank_tips(col, name):
    return frozenset(t for t in tree.tips if lab[t][col] == name)


CROWN = {}                                           # class crowns for the two fish sets, phylum crowns otherwise
for cl in FOCAL:
    ts = rank_tips("class", cl) if cl in CLASS_FOCAL else rank_tips("phylum", cl)
    CROWN[cl] = mrca(tree, ts)
    assert len(CROWN[cl].tips) == len(ts), (cl, len(ts), len(CROWN[cl].tips))
    print("crown %-14s focal %3d taxon %3d D/S/U %s crowns %2d origins %2d age %.3f ci %s" % (
        cl, len(FOC[cl]), len(ts), "/".join(str(CNT[cl][k]) for k in "DSU"), NCR[cl], NOR[cl], CROWN[cl].age,
        CROWN[cl].ci))


def ptips(*ps):
    return [t for t in tree.tips if phy(t) in ps]


BIL = mrca(tree, ptips("Mollusca", "Chordata"))
LOPH = mrca(tree, ptips("Mollusca", "Annelida", "Nemertea", "Brachiopoda", "Bryozoa"))
assert any(c is prot for c in BIL.children)
ECD = [c for c in prot.children if not has(c, "Mollusca")][0]
DEU = [c for c in BIL.children if not has(c, "Mollusca")][0]
HIGHER = [("Metazoa", tree), ("Bilateria", BIL), ("Spiralia" if has(LOPH, "Rotifera") else "Lophotrochozoa", LOPH),
          ("Ecdysozoa", ECD), ("Deuterostomia", DEU)]
print("phyla", " ".join("%s:%d" % kv for kv in sorted(Counter(phy(t) for t in tree.tips).items())))
for nm, n in HIGHER:
    print("higher %-15s tips %3d age %.2f y %.1f phyla %s" % (nm, len(n.tips), n.age, n.y,
                                                               ",".join(sorted({phy(t) for t in n.tips}))))
fig = plt.figure(figsize=(FW * MM, FH * MM), dpi=300)
ax = mm_axes(fig, FW, FH)
R = fig.canvas.get_renderer()
FS, CELL = 5, 4.0                                    # panel text (pt; terr Fig. 1 uses one size) and grid cell (mm)
ITEMS, GRID, TXT = [], defaultdict(list), []         # drawn ink for label placement: (kind, geometry, box)


def reg(kind, g):                                    # kind "rect": (x0, y0, x1, y1); "seg" and "tri": points
    b = tuple(g) if kind == "rect" else (min(p[0] for p in g), min(p[1] for p in g), max(p[0] for p in g),
                                          max(p[1] for p in g))
    ITEMS.append((kind, g, b))
    for c in _cells(b):
        GRID[c].append(len(ITEMS) - 1)
    return len(ITEMS) - 1


def _cells(b):
    for i in range(int(b[0] // CELL), int(b[2] // CELL) + 1):
        for j in range(int(b[1] // CELL), int(b[3] // CELL) + 1):
            yield i, j


def near(b):
    return {k for c in _cells(b) for k in GRID.get(c, ())}


def seg_box(p, q, b):                                # Liang-Barsky: does segment pq meet box b?
    dx, dy, t0, t1 = q[0] - p[0], q[1] - p[1], 0.0, 1.0
    for pp, qq in ((-dx, p[0] - b[0]), (dx, b[2] - p[0]), (-dy, p[1] - b[1]), (dy, b[3] - p[1])):
        if pp == 0:
            if qq < 0:
                return False
            continue
        t = qq / pp
        if pp < 0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return False
    return True


def _side(a, b, p):
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])


def in_tri(p, T):
    s = [_side(T[k], T[(k + 1) % 3], p) for k in range(3)]
    return min(s) >= 0 or max(s) <= 0


def box_tri(b, T):                                   # does box b meet triangle T?
    if any(in_tri(p, T) for p in ((b[0], b[1]), (b[2], b[1]), (b[0], b[3]), (b[2], b[3]))):
        return True
    if any(b[0] <= p[0] <= b[2] and b[1] <= p[1] <= b[3] for p in T):
        return True
    return any(seg_box(T[k], T[(k + 1) % 3], b) for k in range(3))


def hit_items(b, pad=0.3, ign=()):                   # the drawn items that the padded box b meets
    bp = (b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)
    out = []
    for k in sorted(near(bp)):
        if k in ign:
            continue
        kind, g, e = ITEMS[k]
        if e[2] < bp[0] or e[0] > bp[2] or e[3] < bp[1] or e[1] > bp[3]:
            continue
        if kind == "rect" or (kind == "seg" and seg_box(g[0], g[1], bp)) or (kind == "tri" and box_tri(bp, g)):
            out.append(k)
    return out


def hits(b, pad=0.3, ign=()):                        # number of drawn items that the padded box b meets
    return len(hit_items(b, pad, ign))


def seg_seg(p, q, a, b):                             # proper crossing of segments pq and ab
    return _side(a, b, p) * _side(a, b, q) < 0 and _side(p, q, a) * _side(p, q, b) < 0


def crossings(p, q, ign=()):                         # (branches crossed, vertical 1 / horizontal 2; boxes and triangles met)
    e = (min(p[0], q[0]), min(p[1], q[1]), max(p[0], q[0]), max(p[1], q[1]))
    nb, ng = 0.0, 0
    for k in near(e):
        if k in ign:
            continue
        kind, g, _ = ITEMS[k]
        if kind == "seg" and seg_seg(p, q, g[0], g[1]):
            nb += 1.0 if abs(g[0][0] - g[1][0]) < 1e-9 else 2.0
        elif kind == "rect" and seg_box(p, q, g):
            ng += 1
        elif kind == "tri" and (in_tri(p, g) or in_tri(q, g) or any(seg_seg(p, q, g[j], g[(j + 1) % 3])
                                                                    for j in range(3))):
            ng += 1
    return nb, ng


_TS = {}


def tsize(s, ha, va, fs=FS, weight="normal", style="normal"):    # ink box offsets (mm) from the text anchor: x from
    k = (s, ha, va, fs, weight, style)                           # the Agg layout box (it covers the hinted ink), y
    if k not in _TS:                                             # from the glyph outlines about the baseline
        inv = ax.transData.inverted()
        t = ax.text(50.0, 50.0, s, ha=ha, va=va, fontsize=fs, fontweight=weight, style=style)
        tb = ax.text(50.0, 50.0, s, ha=ha, va="baseline", fontsize=fs, fontweight=weight, style=style)
        e, eb = (u.get_window_extent(R).transformed(inv) for u in (t, tb))
        t.remove()
        tb.remove()
        base = 50.0 + max(e.y0, e.y1) - max(eb.y0, eb.y1)        # baseline of t (mm, y downwards)
        p = TextPath((0, 0), s, size=fs, prop=FontProperties(family=plt.rcParams["font.family"], weight=weight,
                                                             style=style)).get_extents()
        assert np.isfinite([p.y0, p.y1]).all() and p.y1 > p.y0, ascii(s)
        kk = 25.4 / 72.0                                         # pt -> mm
        _TS[k] = (min(e.x0, e.x1) - 50.0, base - p.y1 * kk - 50.0, max(e.x0, e.x1) - 50.0, base - p.y0 * kk - 50.0)
    return _TS[k]


def tbox(x, y, o):                                   # ink box (mm) of a text anchored at (x, y)
    return x + o[0], y + o[1], x + o[2], y + o[3]


def put(s, x, y, ha, va, fs=FS, weight="normal", style="normal", color="black", z=8, **kw):
    ax.text(x, y, s, ha=ha, va=va, fontsize=fs, fontweight=weight, style=style, color=color, zorder=z, **kw)
    b = tbox(x, y, tsize(s, ha, va, fs, weight, style))
    TXT.append((s, reg("rect", b)))                  # placed texts, for the overlap check before saving
    return b


def seg(x0, y0, x1, y1, **kw):
    ax.plot([x0, x1], [y0, y1], **kw)
    return reg("seg", [(x0, y0), (x1, y1)])


def place(cands, pad=0.3, ign=()):                   # cands (box, penalty, payload): lowest 100 x hits + penalty
    best = None
    for b, pen, pl in sorted(cands, key=lambda c: c[1]):
        sc = 100 * hits(b, pad, ign) + pen
        if best is None or sc < best[0]:
            best = (sc, b, pl)
        if sc == pen:                                # first free candidate in penalty order (penalties < 100)
            break
    return best


YL = 3.1                                             # baseline of the panel letters and of the b title (terr ED8a)
put("a", 1.0, YL, "left", "baseline", fs=8, weight="bold")
for k, (b, (ya, yz)) in enumerate(zip(blocks, yb)):  # phylum bands, alternately blue and grey (terr)
    band(ax, XT(850.0), XB1, ya, yz, grey=k % 2 == 1)
seg(tree.x - 1.0, tree.y, tree.x, tree.y, color="black", lw=LW, solid_capstyle="butt")
for n in NODES:                                      # branches; unassigned single tips dashed grey (terr outgroups)
    if id(n) in inrow or id(n) in rowid or not n.children:
        continue
    ys_ = [c.y for c in n.children]
    seg(n.x, min(ys_), n.x, max(ys_), color="black", lw=LW, solid_capstyle="projecting")
    for c in n.children:
        kw = dict(color="black", lw=LW, solid_capstyle="butt")
        if not c.children and lab[next(iter(c.tips))]["lab"] == "U":
            kw.update(color=OUTGRP, dashes=(5.5, 3.5))
        seg(n.x, c.y, c.x, c.y, **kw)
for i, r in enumerate(rows):                         # deep clades: ramp triangles; shallow: outlines; deep singles
    mpu = ext[i] / units[i]
    if len(r.tips) > 1:
        h = ext[i] - GAP * mpu
        tri(ax, r.x, r.y, XTIP, r.y - h / 2, r.y + h / 2, fill=st(r) == "D", lw=LW)
        reg("tri", [(r.x, r.y), (XTIP, r.y - h / 2), (XTIP, r.y + h / 2)])
    elif st(r) == "D":
        s_ = min(1.2, max(0.3, ext[i] - GAP * mpu))
        tipsq(ax, XTIP, r.y, s_, lw=0.3)
        reg("rect", (XTIP - s_, r.y - s_ / 2, XTIP, r.y + s_ / 2))
AGE = {}
for cl in FOCAL:                                     # 95% interval bars at the focal crowns (LSD2)
    c = CROWN[cl]
    a, b = sorted(abs(v) for v in c.ci)
    assert a - 1e-6 <= c.age <= b + 1e-6
    AGE[cl] = (c.age, a, b, reg("rect", hpd_bar(ax, XT(b), XT(a), c.y)))
def sight(t, bx, k):                                 # does label box t see its bar bx: at least 0.5 mm of the two
    x0, x1 = max(t[0], bx[0]) + 0.05, min(t[2], bx[2]) - 0.05    # overlap in x and no other ink lies in between?
    if x1 - x0 < 0.4:
        return False
    y0, y1 = (t[3], bx[1]) if t[3] <= bx[1] else (bx[3], t[1])
    return not hits((x0, y0 + 0.05, x1, y1 - 0.05), 0.0, {k})


AGEBAR = {}                                          # age-label text item -> its bar item (contact allowed in the QA)
for cl in FOCAL:                                     # "%.2f Ma" just above its bar (terr), clear of the branches
    age, a, b, k = AGE[cl]
    s, bx = "%.2f Ma" % age, ITEMS[k][2]
    o, xc = tsize(s, "center", "baseline"), 0.5 * (bx[0] + bx[2])
    cands = []
    for gp in np.arange(0.15, 1.401, 0.05):          # gap between the ink and the bar (mm; terr about 0.35)
        for dx in np.arange(-12.0, 12.001, 0.25):
            pen = abs(dx) / 2 + 3 * abs(gp - 0.35)
            yu, yd = bx[1] - gp - o[3], bx[3] + gp - o[1]    # baselines above and below the bar
            for y_, side, p_ in ((yu, "above", pen), (yd, "below", 20 + pen)):   # below only if no place above
                t = tbox(xc + dx, y_, o)
                cands.append((t, p_ + (0 if sight(t, bx, k) else 40), (xc + dx, y_, side, gp, dx)))
    for pd, lim in ((0.3, 40), (0.22, 40), (0.3, 100), (0.22, 100), (0.15, 100)):   # sight first, then 0.3 mm
        sc, bb, (x_, y_, side, gp, dx) = place(cands, pad=pd, ign={k})              # to the line centres (terr
        if sc < lim:                                                                 # spacing), less if crowded
            break
    if sc >= 100:
        print("WARN age label overlaps:", cl, " ".join("%s%s" % (ITEMS[j][0], tuple(round(v, 1) for v in ITEMS[j][2]))
                                                      for j in hit_items(bb, pd, {k})[:4]))
    print("age %-15s %s gap %.2f dx %+.2f pad %.2f score %.1f%s" % (cl, side, gp, dx, pd, sc,
                                                                     "" if sc < 40 else " (no sight)"))
    put(s, x_, y_, "center", "baseline", color=AGE_TXT)
    AGEBAR[TXT[-1][1]] = k
put("Metazoa", tree.x - 1.4, tree.y, "right", "center_baseline")     # left of the stem, no leader (terr)
def lead_txt(p, q, pad=0.3):                         # does leader pq pass within pad mm of a placed text?
    e_ = (min(p[0], q[0]) - pad, min(p[1], q[1]) - pad, max(p[0], q[0]) + pad, max(p[1], q[1]) + pad)
    for j in near(e_):
        if j in TXK:
            t = ITEMS[j][2]
            if seg_box(p, q, (t[0] - pad, t[1] - pad, t[2] + pad, t[3] + pad)):
                return True
    return False


for nm, n in sorted(HIGHER[1:], key=lambda h: h[1].y):  # name beside its node with a short grey leader (terr)
    best, TXK = None, {k for _, k in TXT}
    for side, vd in ((-1, 1), (-1, -1), (1, 1), (1, -1)):  # left or right; label below (1) or above (-1)
        ha, va = ("right" if side < 0 else "left"), ("top" if vd > 0 else "bottom")
        o = tsize(nm, ha, va)
        for dx in np.arange(1.5, 20.01, 0.5):        # far enough to clear the deeper verticals (terr Lophotrochozoa)
            for dy in np.arange(0.8, 4.01, 0.4):
                x_, y_ = n.x + side * dx, n.y + vd * dy
                b = tbox(x_, y_, o)
                if b[0] < 0.5 or b[2] > XTIP:
                    continue
                p = (b[2] + 0.5 if side < 0 else b[0] - 0.5, 0.5 * (b[1] + b[3]))   # leader 0.5 mm off the text
                q = (n.x + 0.15 * side, n.y + 0.15 * vd)
                nb, ng = crossings(p, q)
                ang = np.degrees(np.arctan2(abs(q[1] - p[1]), abs(q[0] - p[0]) + 1e-9))
                ln = np.hypot(q[0] - p[0], q[1] - p[1])
                sc = (100 * hits(b, 0.2) + 3 * nb + 10 * ng + abs(ang - 30) / 10 + abs(ln - 3) +
                      (6 if side > 0 else 0) + (5 if ang < 12 or ang > 70 else 0))
                if best is not None and sc >= best[0]:
                    continue
                sc += (50 * (hits((p[0], p[1], p[0], p[1]), 0.4) > 0) +   # leader starts clear of other ink
                       20 * any(j in TXK for j in hit_items(b, 1.0)) +   # 1 mm from the other texts
                       20 * lead_txt(p, q))                              # leader clear of the other texts
                if best is None or sc < best[0]:
                    best = (sc, x_, y_, ha, va, p, q, b)
    sc, x_, y_, ha, va, p, q, b = best
    if hits(b, 0.2):
        print("WARN higher-taxon label overlaps:", nm)
    print("higher %-16s score %.1f" % (nm, sc))
    put(nm, x_, y_, ha, va)
    seg(p[0], p[1], q[0], q[1], color=LEADER, lw=0.4, solid_capstyle="butt")
SPAN = {}
for cl in FOCAL:                                     # the rows of each focal crown bound its arrow
    ys_ = [(rows[i].y - 0.5 * ext[i], rows[i].y + 0.5 * ext[i]) for i in range(len(rows))
           if rows[i].tips <= CROWN[cl].tips]
    SPAN[cl] = (min(a for a, _ in ys_), max(z for _, z in ys_))
HB, HGAP = 7.8, 0.8                                  # node box height (2.6 + 2 x 2.6) and the gap between boxes
order = sorted(FOCAL, key=lambda c: CROWN[c].y)
lo_ = {c: max(SPAN[c][0] + 0.9, 1.0 + 2.6 + 1.3) for c in order}   # every box top stays inside the figure
hi_ = {c: SPAN[c][1] - 0.9 for c in order}
ya_ = {c: min(max(CROWN[c].y, lo_[c]), hi_[c]) for c in order}
for _ in range(400):                                 # push the arrows apart within their crown rows
    for p, q in zip(order, order[1:]):
        ya_[q] = min(max(ya_[q], ya_[p] + HB + HGAP), hi_[q])
    for p, q in zip(order[::-1][1:], order[::-1]):
        ya_[p] = max(min(ya_[p], ya_[q] - HB - HGAP), lo_[p])
assert all(ya_[q] - ya_[p] >= HB + HGAP - 1e-6 for p, q in zip(order, order[1:])), ya_
from PIL import Image
PP = os.path.join(ROOT, "figures_src_sty", "phylopic")
PNG = {f.split("__")[0]: os.path.join(PP, f) for f in os.listdir(PP) if f.endswith(".png")}
assert set(PNG) == set(FOCAL), sorted(PNG)


def aspect(path):                                    # width / height of the silhouette ink
    a = np.asarray(Image.open(path).convert("LA"))[:, :, 1]
    ys, xs = np.nonzero(a > 5)
    return (xs.max() - xs.min() + 1) / (ys.max() - ys.min() + 1)


BOXES = {}
for cl in order:                                     # silhouette, "Clade (n)", grey arrow, node box (terr Rotifera)
    ya = ya_[cl]
    b = put("%s (%d)" % (cl, len(FOC[cl])), XLAB, ya, "right", "center_baseline")
    x0s, x1s, hs = silhouette(ax, PNG[cl], b[0] - 0.6, ya, h=min(3.4, 6.0 / aspect(PNG[cl])), ha="right")
    reg("rect", (x0s, ya - 0.5 * hs, x1s, ya + 0.5 * hs))
    arrow(ax, b[2] + 0.4, XBOX, ya)
    reg("seg", [(b[2] + 0.4, ya), (XBOX, ya)])
    top = ya - 2.6 - 1.3                             # the arrow meets the middle of the first row (terr)
    node_box(ax, XBOX, top, WBOX, cl, [("%d" % CNT[cl]["D"], DEEP), ("%d" % CNT[cl]["S"], SHAL_TXT),
                                       ("%d" % CNT[cl]["U"], NULL)], [("%d" % NCR[cl], DEEP2), ("%d" % NOR[cl], INK)],
             fs=FS)
    BOXES[cl] = (top, top + HB)
for cl in FOCAL:                                     # analysis-clade tips outside the drawn crown (for the legend)
    xs_ = FOC[cl] - CROWN[cl].tips
    if xs_:
        print("extra %-14s %d: %s" % (cl, len(xs_), " ".join(sorted(Counter(
            "%s/%s" % (phy(t), lab[t]["class"]) for t in xs_).elements()))))
FOCB = {b[0] for b in blocks if b[0] in FOCAL}       # focal phyla: their "Clade (n)" label names the band
for b, (ya, yz) in zip(blocks, yb):                  # "Phylum (n)" right-aligned at the band end (terr)
    if b[0] in FOCB:
        continue
    n = sum(len(rows[i].tips) for i in range(b[1], b[2] + 1))
    s, yc = "%s (%d)" % (b[0], n), 0.5 * (ya + yz)
    o = tsize(s, "right", "center_baseline")
    lo2, hi2 = (ya + 1.0, yz - 1.0) if yz - ya > 3.0 else (ya - 0.5, yz + 0.5)
    yy = sorted(set([yc] + list(np.arange(lo2, hi2 + 1e-9, 0.25))), key=lambda v: abs(v - yc))
    sc, bb, y_ = place([(tbox(XLAB, v, o), abs(v - yc), v) for v in yy], pad=0.2)
    if sc >= 100:
        print("WARN band label overlaps:", b[0], "%.2f-%.2f" % (bb[1], bb[3]), " ".join(
            "%s%s" % (ITEMS[j][0], tuple(round(v, 2) for v in ITEMS[j][2])) for j in hit_items(bb, 0.2)[:4]))
    put(s, XLAB, y_, "right", "center_baseline")


def gen(t):
    s_ = lab[t]["species"].replace("_", " ").split()
    return s_[0] if s_ else ""


RK = ("class", "order", "family", "genus")
SETS = {k: defaultdict(set) for k in RK}
for t in tree.tips:
    for k in RK:
        SETS[k][gen(t) if k == "genus" else lab[t][k]].add(t)


def tipname(r):                                      # highest rank whose sampled species are exactly the row,
    shared = (None, False)                           # else the lowest rank the row shares; genus in italics
    for k in RK:
        vs = {gen(t) if k == "genus" else lab[t][k] for t in r.tips}
        v = vs.pop() if len(vs) == 1 else None
        if v is None or v in NA:
            continue
        if SETS[k][v] == set(r.tips):
            return v, k == "genus"
        shared = (v, k == "genus")
    return shared


BAND = {i: k for k, b in enumerate(blocks) for i in range(b[1], b[2] + 1)}
SKIP = set(FOCAL) | {b[0] for b in blocks}
NAMED = []
for i in sorted(range(len(rows)), key=lambda i: ({"D": 0, "S": 1}.get(st(rows[i]), 9), -len(rows[i].tips), i)):
    r = rows[i]
    if st(r) not in ("D", "S"):
        continue                                     # unassigned tips stay unnamed (dashed grey, terr outgroups)
    nm, it = tipname(r)
    k = BAND[i]
    if nm is None or nm in SKIP or blocks[k][1] == blocks[k][2]:
        continue                                     # a one-row band is named by its band label; names unique (terr)
    col, cnt = (DEEP if st(r) == "D" else "black"), len(r.tips) > 1 or st(r) == "D"
    sty = "italic" if it else "normal"
    o1, o2 = tsize(nm, "left", "center_baseline", style=sty), tsize("(%d)" % len(r.tips), "left", "center_baseline")
    x1, x2 = XTIP + 0.6, XTIP + 0.6 + o1[2] + 0.45   # name just right of the tips, count one space after it
    yb0, yb1 = (min(o1[1], o2[1]), max(o1[3], o2[3])) if cnt else (o1[1], o1[3])
    b = (x1 + o1[0], r.y + yb0, (x2 + o2[2]) if cnt else (x1 + o1[2]), r.y + yb1)
    TXK = {kk for _, kk in TXT}                      # terr stacks its tip names with nearly touching ink,
    if b[2] > XB1 - 0.3 or any(j in TXK for j in hit_items(b, 0.15)) or any(j not in TXK for j in hit_items(b, 0.4)):
        continue                                     # and keeps them clear of the icons and branches
    put(nm, x1, r.y, "left", "center_baseline", style=sty, color=col)
    if cnt:
        put("(%d)" % len(r.tips), x2, r.y, "left", "center_baseline", color=col)
    NAMED.append((st(r), nm, len(r.tips)))
    SKIP.add(nm)                                     # the largest row of a taxon carries its name
for top, bot in BOXES.values():                      # the node boxes, for the key and the checks below
    reg("rect", (XBOX, top, XBOX + WBOX, bot))
YT = Y1 + 0.5                                        # top of the time bar
nt0, nl0 = len(ax.texts), len(ax.lines)
i0, _, i1, _ = ink_mm("50 Myr", FS)                  # terr: '50 Myr' bar at the far left, its text under the
SX = 1.0 + 0.5 * (i1 - i0) - 0.5 * (XT(0.0) - XT(50.0))   # bar, ink starting under the panel letter
used, LV = ics_bar(ax, XT, 850.0, YT, rh=2.43, fs=FS, scale_x=SX)
BAR = reg("rect", (XT(850.0), YT, XT(0.0), LV["yr"]))    # era and period rows
ROWLAB = set()
for t in ax.texts[nt0:]:                             # the time-bar texts and lines join the checks below
    x_, y_ = t.get_position()
    k_ = reg("rect", tbox(x_, y_, tsize(t.get_text(), t.get_ha(), t.get_va(), t.get_fontsize())))
    TXT.append((t.get_text(), k_))
    if ITEMS[k_][2][3] < LV["yr"]:                   # row labels sit on the rows by design
        ROWLAB.add(k_)
for l_ in ax.lines[nl0:]:
    xs, ys = l_.get_xdata(), l_.get_ydata()
    reg("seg", [(xs[0], ys[0]), (xs[-1], ys[-1])])


def tw(s):                                           # ink width (mm) of a panel-a text
    o = tsize(s, "left", "baseline")
    return o[2] - o[0]


XR = 178.5                                           # right edge of the key
r1 = [("Deep", DEEP), ("Shallow", SHAL_TXT), ("Unassigned", NULL)]
r2 = [("Deep crowns", DEEP2), ("Disjoint origins", INK)]
w1 = [tw(s) + 1.6 for s, _ in r1]
wk = max(sum(w1), 2 * max(tw(s) + 1.6 for s, _ in r2))   # row 2 splits the box into two equal cells
w1 = [w * wk / sum(w1) for w in w1]
ktop = YT - 7.3 - HB                                 # key bottom 7.3 mm above the time bar (terr)
kb = (XR - wk, ktop, XR, ktop + HB)
if hits(kb, 0.3):
    print("WARN key box overlaps")
node_box(ax, XR - wk, ktop, wk, "Focal clade", r1, r2, w1=w1, fs=FS)
reg("rect", kb)
xi, yi = XR - wk + 5.4, YT + LV["rh"]                # key icons centred on the era/period boundary (terr)
CAP = "Deep-sea lineages (collapsed clades and single species) indicated in blue"
cb = put(CAP, xi + 7.35 + 2.06, LV["yb"], "right", "baseline")   # boundary-age baseline, ending at the square
if hits(cb, 0.2, ign={len(ITEMS) - 1}):
    print("WARN key caption overlaps")
key_icons(ax, xi, yi, tri_w=7.35, tri_h=5.06, sq=1.76, gap=2.06)
reg("tri", [(xi, yi), (xi + 7.35, yi - 2.53), (xi + 7.35, yi + 2.53)])
reg("rect", (xi + 9.41, yi - 0.88, xi + 11.17, yi + 0.88))
LX, AX0, AY0, AW, AH = 138.0, 145.5, 5.2, 30.7, 26.0  # b: letter, axes (mm; the '1,000' label ends at the key edge)
put("b", LX, YL, "left", "baseline", fs=8, weight="bold")
TB = put("Distribution of deep-origin ages", AX0, YL, "left", "baseline", fs=6.5)   # title at the y axis (ED8a)
assert TB[2] <= FW - 0.5, TB
ages = np.array([fnum(r["stem_age"]) for r in ORG])
assert len(ages) == 65 and np.isfinite(ages).all() and 2.0 <= ages.min() and ages.max() <= 800.0, ages
rng = "%.1f–%.0f" % (ages.min(), ages.max())
assert rng == "2.3–679", ascii(rng)
axb = axes_mm(fig, FW, FH, AX0, AY0, AW, AH)
red = ed8_hist(axb, ages, np.geomspace(2.0, 800.0, 21), np.median(ages),
               ["Median = %.1f Ma" % np.median(ages), "Range = %s Ma" % rng],
               xlabel="Stem age (Ma)", log=True, fs=5.5, lab_fs=5.5, tick_fs=5)
axb.set_xlim(1.0, 1000.0)
axb.set_xticks([1, 10, 100, 1000])
axb.set_xticklabels(["1", "10", "100", "1,000"])
axb.minorticks_off()
axb.yaxis.set_major_locator(plt.MaxNLocator(4, integer=True))
for sp in axb.spines.values():                       # black axes and labels (ED8a)
    sp.set_color("black")
axb.tick_params(colors="black", labelcolor="black")
axb.xaxis.label.set_color("black")
axb.yaxis.label.set_color("black")
fig.canvas.draw()
e1, e2 = red.get_window_extent(), axb.get_window_extent()
assert e1.x0 >= e2.x0 and e1.x1 <= e2.x1 and e1.y0 >= e2.y0 and e1.y1 <= e2.y1, (e1, e2)
tb = axb.get_tightbbox(fig.canvas.get_renderer())
BX = (tb.x0 / fig.dpi / MM, FH - tb.y1 / fig.dpi / MM, tb.x1 / fig.dpi / MM, FH - tb.y0 / fig.dpi / MM)
assert BX[0] >= XBOX + WBOX + 1.0 and BX[2] <= XR + 0.1 and BX[3] <= FH, BX   # b ends with the key (XR)
print("panel b ink (mm): %.1f %.1f %.1f %.1f; title ends %.1f" % (BX + (TB[2],)))
reg("rect", BX)
RT = {r.tips for r in rows}
NEW = {"tips": "%d" % len(tree.tips), "crowns_drawn": "%d" % sum(1 for d in dsets if d in RT),
       "origins": "%d" % len(ORG), "origin_age_median": "%.1f" % np.median(ages),
       "origin_age_range": "%.2f-%.1f" % (ages.min(), ages.max())}
OLD = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(BUILD, "fig1_facts.tsv"), encoding="utf-8")
           if "\t" in l)
bad = [k for k in NEW if OLD.get(k) != NEW[k]]
assert not bad, {k: (OLD.get(k), NEW[k]) for k in bad}  # same numbers as fig1_sty.py
out = [[k, v] for k, v in NEW.items()]
out += [["rows", "%d" % len(rows)], ["blocks", "%d" % len(blocks)], ["tip_labels", "%d" % len(NAMED)],
        ["tip_labels_deep", "%d" % sum(1 for s_, _, _ in NAMED if s_ == "D")]]
for cl in FOCAL:
    out.append([cl, "n=%d D=%d S=%d U=%d crowns=%d origins=%d crown_tips=%d crown_age=%.2f ci=%.2f-%.2f" % (
        len(FOC[cl]), CNT[cl]["D"], CNT[cl]["S"], CNT[cl]["U"], NCR[cl], NOR[cl], len(CROWN[cl].tips),
        *AGE[cl][:3])])
out.append(["higher", " ".join("%s=%.2f" % (nm, n.age) for nm, n in HIGHER)])
out.append(["ics_short", " ".join("%s=%s" % p for p in used)])
out.append(["named", " ".join("%s|%s|%d" % x for x in NAMED)])
TS_ = {k: s for s, k in TXT}


def meets(b, j):                                     # does item j meet box b (the test of hits, for one item)?
    kind, g, e = ITEMS[j]
    if e[2] < b[0] or e[0] > b[2] or e[3] < b[1] or e[1] > b[3]:
        return False
    return kind == "rect" or (kind == "seg" and seg_box(g[0], g[1], b)) or (kind == "tri" and box_tri(b, g))


nt = no = 0                                          # placed texts that meet another text or drawn ink
for s, k in TXT:
    b = ITEMS[k][2]
    bs = (b[0] + 0.05, b[1] + 0.05, b[2] - 0.05, b[3] - 0.05)   # ink boxes shrunk by 0.05 mm: real contact only
    m = sorted(j for j in set(near(bs)) if j != k and meets(bs, j) and not (k in ROWLAB and j == BAR))
    if m:
        nt += any(j in TS_ for j in m)
        no += 1
        print("OVERLAP %s -> %s" % (ascii(s), ", ".join("%s%s" % (ascii(TS_.get(j, ITEMS[j][0])), tuple(
            round(v, 1) for v in ITEMS[j][2])) for j in m[:3])))
edge = [ITEMS[j][2] for j in range(len(ITEMS)) if min(ITEMS[j][2][0], ITEMS[j][2][1], FW - ITEMS[j][2][2],
                                                       FH - ITEMS[j][2][3]) < 0.3]
print("texts %d, overlapping %d (with texts %d); items near the edge %d %s" % (len(TXT), no, nt, len(edge),
                                                                            edge[:3]))
out.append(["qa", "texts=%d overlapping=%d near_edge=%d" % (len(TXT), no, len(edge))])


def gap_pt_seg(p, a, b):                             # distance (mm) from point p to segment ab
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = 0.0 if dx == dy == 0 else min(1.0, max(0.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)))
    return float(np.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy))


def gap_box(b, e):                                   # distance (mm) between boxes b and e
    return float(np.hypot(max(e[0] - b[2], 0.0, b[0] - e[2]), max(e[1] - b[3], 0.0, b[1] - e[3])))


def gap_seg(b, a, c):                                # distance (mm) between box b and segment ac that do not meet
    cor = ((b[0], b[1]), (b[2], b[1]), (b[0], b[3]), (b[2], b[3]))
    return min([gap_box(b, (p[0], p[1], p[0], p[1])) for p in (a, c)] + [gap_pt_seg(p, a, c) for p in cor])


def gap(b, j):                                       # geometric gap between text box b and item j (0 if they meet)
    kind, g, e = ITEMS[j]
    if meets(b, j):
        return 0.0
    if kind == "rect":
        return gap_box(b, e)
    if kind == "seg":
        return gap_seg(b, g[0], g[1])
    return min(gap_seg(b, g[i], g[(i + 1) % 3]) for i in range(3))


CLR = []                                             # ink gaps below 1 mm of every placed text (strokes: 0.07 mm)
for s, k in TXT:
    b = ITEMS[k][2]
    for j in set(near((b[0] - 1.0, b[1] - 1.0, b[2] + 1.0, b[3] + 1.0))):
        if j == k or (j in TS_ and j < k) or AGEBAR.get(k) == j or (k in ROWLAB and (j == BAR or j in ROWLAB)):
            continue
        d = gap(b, j) - (0.0 if j in TS_ else 0.07)
        if d < 1.0:
            CLR.append((max(d, 0.0), ascii(s), ascii(TS_.get(j, ITEMS[j][0])), j in TS_))
CLR.sort()
for d, s, t, tt in CLR[:14]:
    print("clear %.2f mm  %s -> %s%s" % (d, s, t, " (text)" if tt else ""))
ci = [d for d, _, _, tt in CLR if not tt] or [1.0]
ct = [d for d, _, _, tt in CLR if tt] or [1.0]
out.append(["clearance", "text_ink_min=%.2f text_text_min=%.2f text_ink_below_0.2=%d text_text_below_0.2=%d" % (
    min(ci), min(ct), sum(d < 0.2 for d in ci), sum(d < 0.2 for d in ct))])
print("clearance:", out[-1][1])
write_facts("fig1_facts_v2.tsv", out)
save_v2(fig, "Fig1")
print("saved Fig1; named %d (deep %d); box tops %s; key %.1f-%.1f x %.1f-%.1f" % (
    len(NAMED), sum(1 for s_, _, _ in NAMED if s_ == "D"), " ".join("%.1f" % BOXES[c][0] for c in order),
    kb[0], kb[2], kb[1], kb[3]))
