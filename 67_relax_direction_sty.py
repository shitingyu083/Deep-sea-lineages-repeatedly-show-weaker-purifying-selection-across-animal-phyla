"""67 (post hoc, added 2026-10-29 after external critique): direction-specific RELAX cross-clade test.
Re-uses the 25b paired sign-flip permutation (real vs null1, nullhy HOGs, N=10,000, seed 20261022) and the
25c weighted Stouffer (w = sqrt(n_paired), p clipped at 1e-4) on three indicators per HOG:
  any   = RELAX p < 0.05                      (must reproduce 25b 'RELAX p<0.05' paired p exactly: integrity check)
  relax = RELAX p < 0.05 and K < 1            (upper tail: real > null)
  intens= RELAX p < 0.05 and K > 1            (lower tail: real < null, i.e. depletion)
Usage: python 67_relax_direction_sty.py OUTDIR CLADE=codon_dir [...]
"""
import csv, math, os, sys
import numpy as np
from scipy.stats import norm

NPERM, EPS = 10000, 1e-4


def f(x):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return float("nan")


def truthy(x):
    return str(x).strip() in ("1", "1.0", "True", "true")


def ind(r, kind):
    p, k = f(r["relax_p"]), f(r["relax_K"])
    if p != p:
        return float("nan")
    s = p < 0.05
    if kind == "any":
        return float(s)
    if not s:
        return 0.0
    if k != k:
        return float("nan")
    return float(k < 1) if kind == "relax" else float(k > 1)


def perm(d, rng):
    d = np.asarray(d, float); T = d.mean()
    s = rng.choice([-1.0, 1.0], size=(NPERM, len(d)))
    Ts = (s * d).mean(1)
    return T, (1 + np.sum(Ts >= T - 1e-12)) / (1 + NPERM), (1 + np.sum(Ts <= T + 1e-12)) / (1 + NPERM)


def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    rows = []
    for a in sys.argv[2:]:
        clade, d = a.split("=", 1)
        by = {}
        with open(os.path.join(d, "sel_results_v3.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                by.setdefault(r["hog"], {})[r["set"]] = r
        for kind in ("any", "relax", "intens"):
            rng = np.random.default_rng(20261022)   # fresh per clade x kind; 'any' then equals 25b's first test
            diffs, nr, nn = [], [], []
            for s in by.values():
                if "real" in s and "null1" in s and truthy(s["real"].get("in_nullhy")):
                    a_, b_ = ind(s["real"], kind), ind(s["null1"], kind)
                    if a_ == a_ and b_ == b_:
                        diffs.append(a_ - b_); nr.append(a_); nn.append(b_)
            T, pu, pl = perm(diffs, rng)
            rows.append(dict(clade=clade, kind=kind, n_paired=len(diffs), real_rate=np.mean(nr), null1_rate=np.mean(nn),
                             T=T, p_upper=pu, p_lower=pl))
    with open(os.path.join(out, "perclade.tsv"), "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    comb = []
    for kind, tail in (("any", "p_upper"), ("relax", "p_upper"), ("intens", "p_lower")):
        rr = [r for r in rows if r["kind"] == kind]
        ps = [min(max(r[tail], EPS), 1 - EPS) for r in rr]; ws = [math.sqrt(r["n_paired"]) for r in rr]
        z = sum(wi * norm.isf(p) for wi, p in zip(ws, ps)) / math.sqrt(sum(wi * wi for wi in ws))
        loco = []
        for i in range(len(ps)):
            pp, ww = ps[:i] + ps[i + 1:], ws[:i] + ws[i + 1:]
            loco.append(sum(wi * norm.isf(p) for wi, p in zip(ww, pp)) / math.sqrt(sum(wi * wi for wi in ww)))
        comb.append(dict(kind=kind, tail=tail, k=len(ps), n_sig=sum(p < 0.05 for p in ps), z=z, p=norm.sf(z),
                         loco_zmin=min(loco), loco_zmax=max(loco)))
    with open(os.path.join(out, "combined.tsv"), "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(comb[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(comb)
    print(open(os.path.join(out, "perclade.tsv")).read()); print(open(os.path.join(out, "combined.tsv")).read())


if __name__ == "__main__":
    main()
