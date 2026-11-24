"""Shared style, IO and a minimal Newick parser for the manuscript figures.

Nature print specification: 180 mm (double column) or 89 mm (single), Arial, 5-7 pt text, 0.5 pt lines.
Fixed palette (colour-blind safe, prints to distinct greys): deep #1f4e79, shallow #d4a03a, null #888888,
real/ink #333333, light #cccccc.
"""
import csv, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data_srv")
BUILD = os.path.join(ROOT, "build")
OUT = os.path.join(ROOT, "submission_ready_20261029")
W = lambda *p: os.path.join(DATA, "work", *p)
MM = 1 / 25.4
DEEP, SHAL, NULL, INK, LIGHT = "#1f4e79", "#d4a03a", "#888888", "#333333", "#cccccc"
CLADES = ["Mollusca", "Cnidaria", "Echinodermata", "Annelida", "Actinopteri_m", "Arthropoda", "Chondrichthyes", "Porifera"]
CLADE_LABEL = {"Actinopteri_m": "Actinopteri", "Actinopteri": "Actinopteri"}
TRAIT_LABEL = {"busco_f_T": "Thr frequency ($f_T$)", "whole_pI": "Proteome pI", "whole_R_RK": "Arg/(Arg+Lys)",
               "whole_acid_base": "Acidic/basic ratio", "whole_net_charge7_per_res": "Net charge (pH 7)",
               "busco_f_F": "Phe frequency ($f_F$)"}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5, "xtick.labelsize": 6,
    "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "legend.frameon": False, "axes.edgecolor": INK, "text.color": INK,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})


def rtsv(path):
    with open(path, encoding="utf-8") as f:
        rows = [r for r in csv.reader(f, delimiter="\t") if r and not r[0].startswith("#")]
    h = [c.strip('"') for c in rows[0]]
    return [dict(zip(h, [c.strip('"') for c in r])) for r in rows[1:]]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def panel(ax, letter, x=-0.18, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def save(fig, folder, name):
    d = os.path.join(OUT, folder)
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, name + ".pdf"))
    fig.savefig(os.path.join(d, name + ".png"), dpi=300)
    fig.savefig(os.path.join(d, name + ".tif"), dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


class Node:
    __slots__ = ("name", "length", "children", "parent", "id", "x", "y", "tips")

    def __init__(self):
        self.name, self.length, self.children, self.parent, self.id = "", 0.0, [], None, None


def parse_newick(s):
    s = s.strip().rstrip(";")
    root = Node(); cur = root; i = 0; tok = re.compile(r"[^(),:;]+")
    while i < len(s):
        c = s[i]
        if c == "(":
            n = Node(); n.parent = cur; cur.children.append(n); cur = n; i += 1
        elif c == ",":
            n = Node(); n.parent = cur.parent; cur.parent.children.append(n); cur = n; i += 1
        elif c == ")":
            cur = cur.parent; i += 1
        elif c == ":":
            m = re.match(r":([-0-9.eE+]+)", s[i:]); cur.length = float(m.group(1)); i += len(m.group(0))
        else:
            m = tok.match(s, i); cur.name = m.group(0); i = m.end()
    return root


def postorder(n):
    st, out = [n], []
    while st:
        x = st.pop(); out.append(x); st.extend(x.children)
    return out[::-1]


def ape_numbering(root):
    """ape node numbers: tips 1..N in Newick order, internal nodes N+1.. in preorder (root = N+1)."""
    tips, internal, st = [], [], [root]
    while st:
        x = st.pop()
        if x.children:
            internal.append(x); st.extend(reversed(x.children))
        else:
            tips.append(x)
    for k, t in enumerate(tips, 1):
        t.id = k
    for k, n in enumerate(internal, len(tips) + 1):
        n.id = k
    return tips, internal
