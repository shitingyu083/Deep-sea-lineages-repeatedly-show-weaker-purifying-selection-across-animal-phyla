#!/usr/bin/env python
"""R5 collect (PREREG_r5.md). D1: real vs null2 (C1), real vs null1 on the same HOGs (C2), null2 vs null1 (C3), with
the frozen 67 machinery (indicators, paired sign-flip N=10,000, seed 20261022 fresh per clade x kind x comparison,
weighted Stouffer). D2a: power and recurrence power of K85/K90/K75h (63c rules). D3: per-row table for the 69e models
plus per-route paired sign-flip/Stouffer. Usage: python 69d_r5_collect_sty.py [D1] [D2] [D3]"""
import os, sys, csv, json, math, itertools, collections, importlib.util
import numpy as np
from scipy.stats import norm
B = os.path.expanduser("~/deepsea_conv"); O = f"{B}/work/addon/R5"
def imp(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
d67 = imp("d67", f"{B}/scripts/67_relax_direction_sty.py"); ps = imp("ps", f"{B}/scripts/63a3_powersim_sty.py")
CL = ["Actinopteri_m", "Annelida", "Arthropoda", "Chondrichthyes", "Cnidaria", "Echinodermata", "Mollusca", "Porifera"]
FEAT = list(csv.DictReader(open(f"{B}/work/addon/A4/sel_features.tsv"), delimiter="\t"))
FE = {(r["clade"], r["hog"], r["set"]): r for r in FEAT}
TAIL = {"any": "p_upper", "relax": "p_upper", "intens": "p_lower", "relax_relaxed": "p_upper"}
WHAT = sys.argv[1:] or ["D1", "D2", "D3"]


def jrelax(p):
    """(relax_p, relax_K) from a RELAX json, as 25b; (nan, nan) if missing/unreadable."""
    try:
        t = json.load(open(p)).get("test results", {})
        return d67.f(t.get("p-value")), d67.f(t.get("relaxation or intensification parameter"))
    except Exception:
        return float("nan"), float("nan")


def stouffer(rows, tail, minn=0):
    rr = [r for r in rows if r["n_paired"] >= minn and r["n_paired"] > 0]
    if not rr: return dict(k=0, z=float("nan"), p=float("nan"), loco_zmin=float("nan"), loco_zmax=float("nan"), pooled_T=float("nan"),
                           n_sig=0, clades="")
    p_ = [min(max(r[tail], d67.EPS), 1 - d67.EPS) for r in rr]; w = [math.sqrt(r["n_paired"]) for r in rr]
    zf = lambda pp, ww: sum(a * norm.isf(b) for a, b in zip(ww, pp)) / math.sqrt(sum(a * a for a in ww))
    z = zf(p_, w); loco = [zf(p_[:i] + p_[i + 1:], w[:i] + w[i + 1:]) for i in range(len(p_))] if len(p_) > 1 else [z]
    pt = sum(r["T"] * r["n_paired"] for r in rr) / sum(r["n_paired"] for r in rr)
    return dict(k=len(rr), z=z, p=norm.sf(z), loco_zmin=min(loco), loco_zmax=max(loco), pooled_T=pt,
                n_sig=sum(x < 0.05 for x in p_), clades=",".join(r["clade"] for r in rr))


def wtsv(fn, rows):
    with open(fn, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)


def ind2(p, k, kind):
    return d67.ind({"relax_p": p, "relax_K": k}, "relax" if kind == "relax_relaxed" else kind)


# ------------------------------------------------------------------ D1
if "D1" in WHAT:
    D = f"{O}/D1"; os.makedirs(D, exist_ok=True); rows, per = [], []
    for C in CL:
        by = {}
        for r in csv.DictReader(open(f"{B}/work/codon/{C}/sel_results_v3.tsv"), delimiter="\t"):
            by.setdefault(r["hog"], {})[r["set"]] = r
        hs = {fn.rsplit(".null2.relax.", 1)[0] for fn in os.listdir(f"{O}/null2/{C}") if ".null2.relax." in fn}  # json or log
        for h in sorted(hs):
            fn = f"{h}.null2.relax.json"; s = by.get(h, {})
            if "real" not in s or "null1" not in s or not d67.truthy(s["real"].get("in_nullhy")): continue
            p2, k2 = jrelax(f"{O}/null2/{C}/{fn}")
            rows.append(dict(clade=C, hog=h, real_p=d67.f(s["real"]["relax_p"]), real_K=d67.f(s["real"]["relax_K"]),
                             null1_p=d67.f(s["null1"]["relax_p"]), null1_K=d67.f(s["null1"]["relax_K"]), null2_p=p2, null2_K=k2))
    wtsv(f"{D}/rows.tsv", rows)
    for C in CL:
        rc = [r for r in rows if r["clade"] == C]
        for kind in ("any", "relax", "intens"):
            v = [(ind2(r["real_p"], r["real_K"], kind), ind2(r["null1_p"], r["null1_K"], kind), ind2(r["null2_p"], r["null2_K"], kind)) for r in rc]
            v = [x for x in v if all(y == y for y in x)]  # same HOG set for C1, C2, C3
            for cmp_, (i, j) in (("C1_real_vs_null2", (0, 2)), ("C2_real_vs_null1", (0, 1)), ("C3_null2_vs_null1", (2, 1))):
                rng = np.random.default_rng(20261022)
                d = [x[i] - x[j] for x in v]
                T, pu, pl = d67.perm(d, rng) if d else (float("nan"), float("nan"), float("nan"))
                per.append(dict(clade=C, kind=kind, comparison=cmp_, n_paired=len(d), rate_a=np.mean([x[i] for x in v]) if v else float("nan"),
                                rate_b=np.mean([x[j] for x in v]) if v else float("nan"), T=T, p_upper=pu, p_lower=pl))
    wtsv(f"{D}/perclade.tsv", per)
    comb = []
    for kind in ("any", "relax", "intens"):
        for cmp_ in ("C1_real_vs_null2", "C2_real_vs_null1", "C3_null2_vs_null1"):
            tail = "p_upper" if cmp_.startswith("C3") else TAIL[kind]
            rr = [r for r in per if r["kind"] == kind and r["comparison"] == cmp_]
            comb.append(dict(kind=kind, comparison=cmp_, tail=tail, **stouffer(rr, tail, minn=30)))
    wtsv(f"{D}/combined.tsv", comb)
    g = {(r["kind"], r["comparison"]): r for r in comb}; V = []
    for kind in ("relax", "intens", "any"):
        c1, c2 = g[(kind, "C1_real_vs_null2")], g[(kind, "C2_real_vs_null1")]
        sgn = -1 if kind == "intens" else 1; t1, t2 = sgn * c1["pooled_T"], sgn * c2["pooled_T"]
        verdict = ("dependent" if not c1["z"] >= 1.645 else "robust" if t1 >= 0.5 * t2 else "robust_sign_attenuated_size")
        V.append(f"D1\t{kind}\tC1_z={c1['z']:.3f}\tC1_T={c1['pooled_T']:.4f}\tC2_z={c2['z']:.3f}\tC2_T={c2['pooled_T']:.4f}\tratio={t1 / t2 if t2 else float('nan'):.3f}\t{verdict}")
    z3 = {k: g[(k, "C3_null2_vs_null1")]["z"] for k in ("any", "relax", "intens")}
    V.append("D1\tcalibration\t" + "\t".join(f"C3_{k}_z={v:.3f}" for k, v in z3.items()) + "\t" +
             ("calibrated" if all(abs(v) < 1.96 for v in z3.values()) else "null_designs_differ"))
    V.append(f"D1\tn_rows={len(rows)}\tnull2_fit_missing={sum(r['null2_p'] != r['null2_p'] for r in rows)}")
    open(f"{D}/verdict.tsv", "w").write("\n".join(V) + "\n"); print("\n".join(V))

# ------------------------------------------------------------------ D2a
if "D2" in WHAT:
    D = f"{O}/D2"; os.makedirs(D, exist_ok=True)
    res = []
    for C in CL:
        for h in sorted(os.listdir(f"{O}/power/{C}")) if os.path.isdir(f"{O}/power/{C}") else []:
            for sc in ("K85", "K90", "K75h"):
                fn = f"{O}/power/{C}/{h}/{sc}/res.tsv"
                if os.path.exists(fn): res += list(csv.DictReader(open(fn), delimiter="\t"))
    wtsv(f"{D}/res_all.tsv", res)
    sel = [r for r in FEAT if r["set"] == "real"]
    fl = lambda x: (lambda v: v if v == v else None)(d67.f(x))
    def qthr(C):  # 63c rule
        v = [(fl(r["relax_p"]), fl(r["relax_q"])) for r in sel if r["clade"] == C and fl(r["relax_p"]) is not None]
        s = [p for p, q in v if q is not None and q < 0.05]
        return (max(s) if s else 0.05 / max(1, len(v)))
    out, power = [], collections.defaultdict(dict)
    for r in csv.DictReader(open(f"{B}/work/addon/R2/main/power_by_clade.tsv"), delimiter="\t"):
        if r["scen"] in ("H0r", "K75", "K50", "K50h"):
            out.append(dict(clade=r["clade"], scen=r["scen"], n=r["n"], rate_p05=r["rate_p05"], q_equiv_p=r["q_equiv_p"],
                            rate_qequiv=r["rate_qequiv"], median_K=r["median_K"], source="R2")); power[r["scen"]][r["clade"]] = float(r["rate_qequiv"])
    for C in sorted(set(r["clade"] for r in res)):
        th = qthr(C)
        for sc in ("K85", "K90", "K75h"):
            rr = [r for r in res if r["clade"] == C and r["scen"] == sc and fl(r["p"]) is not None]
            det = [fl(r["p"]) for r in rr if fl(r["K"]) is not None and fl(r["K"]) < 1]; n = len(rr)
            a05 = sum(p < 0.05 for p in det) / n if n else float("nan"); aq = sum(p < th for p in det) / n if n else float("nan")
            km = sorted(fl(r["K"]) for r in rr if fl(r["K"]) is not None)
            out.append(dict(clade=C, scen=sc, n=n, rate_p05=f"{a05:.3f}", q_equiv_p=f"{th:.3g}", rate_qequiv=f"{aq:.3f}",
                            median_K=f"{km[len(km) // 2]:.3f}" if km else "", source="R5")); power[sc][C] = aq
    wtsv(f"{D}/power_by_clade.tsv", out)
    def pb_ge(p_, m):
        d = [1.0]
        for p in p_: d = [a * (1 - p) + (d[i - 1] * p if i else 0) for i, a in enumerate(d)] + [d[-1] * p]
        return sum(d[m:])
    rec = []
    for sc in ("K85", "K90", "K75h"):
        cl = sorted(c for c in power[sc] if power[sc][c] == power[sc][c])
        for k in range(3, len(cl) + 1):
            v = [pb_ge([power[sc][c] for c in s], 3) for s in itertools.combinations(cl, k)]
            rec.append(dict(scen=sc, k_clades=k, mean_P_ge3=f"{np.mean(v):.3f}", min=f"{min(v):.3f}", max=f"{max(v):.3f}"))
    if rec: wtsv(f"{D}/recurrence_power.tsv", rec)
    print(open(f"{D}/power_by_clade.tsv").read()); print(f"sim results {len(res)} / 480")

# ------------------------------------------------------------------ D3
if "D3" in WHAT:
    D = f"{O}/D3"; os.makedirs(D, exist_ok=True); rows = []
    gapfrac = lambda s: sum(s.count(c) for c in "-Nn?Xx") / len(s)
    for C in CL:
        for rt in ("A", "M"):
            d3 = f"{O}/route/{C}/{rt}"
            if not os.path.exists(f"{d3}/picked.tsv"): continue
            for r in csv.DictReader(open(f"{d3}/picked.tsv"), delimiter="\t"):
                h, s = r["hog"], r["set"]; fe = FE.get((C, h, s), {})
                p, k = jrelax(f"{d3}/{h}.{s}.relax.json")
                if r["reused"] == "1": nt, gd = d67.f(fe.get("n_tips")), d67.f(fe.get("gapdiff"))
                else:
                    t = ps.parse(open(f"{d3}/{h}.{s}.nwk").read()); a = ps.read_fas(f"{d3}/{h}.{s}.fas")
                    tg = {x.name: x.tag for x in ps.nodes(t) if not x.ch}; gp = {q: gapfrac(v) for q, v in a.items() if v}
                    fg = [gp[q] for q in tg if tg[q] == "FG" and q in gp]; bg = [gp[q] for q in tg if tg[q] == "BG" and q in gp]
                    nt = len(tg); gd = (np.mean(fg) - np.mean(bg)) if fg and bg else float("nan")
                rows.append(dict(clade=C, route=rt, hog=h, set=s, frame="pruned", relax_p=p, relax_K=k, n_tips=nt, gapdiff=gd))
                rows.append(dict(clade=C, route=rt, hog=h, set=s, frame="unpruned", relax_p=d67.f(fe.get("relax_p")),
                                 relax_K=d67.f(fe.get("relax_K")), n_tips=d67.f(fe.get("n_tips")), gapdiff=d67.f(fe.get("gapdiff"))))
    wtsv(f"{D}/rows.tsv", rows)
    per, comb = [], []
    for rt in ("A", "M"):
        for fr in ("pruned", "unpruned"):
            for kind in ("relax_relaxed", "any"):
                pr = []
                for C in CL:
                    x = collections.defaultdict(dict)
                    for r in rows:
                        if r["clade"] == C and r["route"] == rt and r["frame"] == fr: x[r["hog"]][r["set"]] = ind2(r["relax_p"], r["relax_K"], kind)
                    d = [v["real"] - v["null1"] for v in x.values() if len(v) == 2 and v["real"] == v["real"] and v["null1"] == v["null1"]]
                    if not d: continue
                    T, pu, pl = d67.perm(d, np.random.default_rng(20261022))
                    pr.append(dict(route=rt, frame=fr, kind=kind, clade=C, n_paired=len(d), T=T, p_upper=pu, p_lower=pl))
                per += pr; comb.append(dict(route=rt, frame=fr, kind=kind, **stouffer(pr, "p_upper")))
    wtsv(f"{D}/perm_perclade.tsv", per); wtsv(f"{D}/perm_combined.tsv", comb)
    miss = sum(1 for r in rows if r["frame"] == "pruned" and r["relax_p"] != r["relax_p"])
    print(f"D3 rows {len(rows)}; pruned fits missing {miss}"); print(open(f"{D}/perm_combined.tsv").read())
