"""Shared helpers for the v4 figures (denser Figs 2-4). Everything of v3 (and through it v2: palette, statistics, QA,
berg/cell helpers, clade colours) is imported unchanged; outputs are redirected to viz_v4_sty/out, never the package.

New benchmark form added in v4 (panel map: paper/figure_benchmark_map.md, section "Fig. 2-4 v4"):
  berg 3c - Bergeron et al. 2023 Nature, doi:10.1038/s41586-023-05752-y, Fig. 3c: black-outline box plots (median,
            IQR, whiskers) with every observation overlaid as an edgeless dot coloured by group, groups separated by
            thin vertical lines, a centred header per group with the test result, and a bracket with stars or NS
            between the two boxes of each group.
"""
import os, sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "viz_v3_sty"))
from figcommon_v3_sty import *  # noqa: F401,F403  (v3 puts viz_v2_sty on the path)
import figcommon_v2_sty as _v2
V4 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(V4, exist_ok=True)
_v2.V2 = V4


def berg3c_box(ax, x, v, w, color, pts=True, rng=None, s=0.9, alpha=0.35, jit=0.38, z=2):
    """berg 3c: one box (black outline, median line, whiskers to the most extreme value within 1.5 IQR, capped),
    with the observations drawn under it as edgeless dots jittered across `jit` of the box width."""
    v = np.asarray(v, float)
    q1, md, q3 = np.quantile(v, [.25, .5, .75])
    lo = v[v >= q1 - 1.5 * (q3 - q1)].min(); hi = v[v <= q3 + 1.5 * (q3 - q1)].max()
    if pts:
        xs = x + (rng.random(len(v)) - 0.5) * jit * 2 * w
        ax.scatter(xs, v, s=s, c=color, alpha=alpha, lw=0, edgecolors="none", zorder=z, rasterized=False)
    ax.add_patch(Rectangle((x - w / 2, q1), w, q3 - q1, fc="none", ec="k", lw=0.6, zorder=z + 2))
    ax.plot([x - w / 2, x + w / 2], [md, md], color="k", lw=0.9, solid_capstyle="butt", zorder=z + 3)
    for a, b in ((q3, hi), (q1, lo)):
        ax.plot([x, x], [a, b], color="k", lw=0.6, solid_capstyle="butt", zorder=z + 2)
        ax.plot([x - w / 4, x + w / 4], [b, b], color="k", lw=0.6, solid_capstyle="butt", zorder=z + 2)
    return q1, md, q3, lo, hi


def berg3c_bracket(ax, x0, x1, y, txt, dy, fs=6.0):
    """berg 3c: bracket with downward ends over the two boxes of a group, stars (or NS) centred above."""
    ax.plot([x0, x0, x1, x1], [y - dy, y, y, y - dy], color="k", lw=0.5, solid_capstyle="butt", zorder=8,
            clip_on=False)
    ax.text(0.5 * (x0 + x1), y + 0.2 * dy, txt, ha="center", va="bottom", fontsize=fs, color=INK, zorder=10,
            clip_on=False)
