#!/usr/bin/env python3
"""T1 (PREREG_addon.md): depth dose-response of the P6 traits under the three label sets, collated from 44's outputs
(no new fit). Label sets: lab (D > 200 m) -> work/final_v3/pgls_controls, lab_strict -> pgls_controls_strict,
lab_1000 (D > 1000 m) -> pgls_controls_1000 (each CONTROLS.DONE).
Models of 44 (controls_summary.tsv, Animalia; est in trait SD units): all = ref_all (deep, all species), ref (deep, on
the species with Thab), deepT (deep + Thab), ldepthT (log10 depth + Thab), deepdose (log10 depth + Thab within the
deep species; term ldepth) and deepdose_Thab (the Thab term of that model). clade_meta.tsv: random-effects
meta-analysis of the per-clade ref_all (deep) and deepdose (ldepth) effects.
Descriptive flag per label set: dose_same_sign = sign(deepdose) == sign(all deep effect).
Usage: python 60_dose_series.py [--out DIR]   (default work/addon/T1)
Writes OUT/dose_series.tsv (long), OUT/dose_wide.tsv (one row per trait), OUT/report.txt, OUT/T1.DONE."""
import csv, os, sys

R = os.environ.get("DSC_ROOT", os.path.expanduser("~/deepsea_conv"))
A = sys.argv
opt = lambda k, d: A[A.index(k) + 1] if k in A else d
OUT = opt("--out", f"{R}/work/addon/T1")
NA = float("nan")
LS = [("lab", "pgls_controls"), ("lab_strict", "pgls_controls_strict"), ("lab_1000", "pgls_controls_1000")]
P6 = ["busco_f_T", "whole_pI", "whole_R_RK", "whole_acid_base", "whole_net_charge7_per_res", "busco_f_F"]
MODELS = [("all", "deep"), ("ref", "deep"), ("deepT", "deep"), ("ldepthT", "ldepth"), ("deepdose", "ldepth"),
          ("deepdose_Thab", "Thab")]
WIDE = ("all", "ldepthT", "deepdose", "deepdose_Thab")
META = [("ref_all", "deep"), ("deepdose", "ldepth")]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NA


def rows(fn):
    with open(fn) as f:
        return list(csv.DictReader(f, delimiter="\t"))


s = lambda v: "NA" if isinstance(v, float) and v != v else str(v)
g = lambda x, d=3: "NA" if x != x else f"{x:.{d}g}"
gi = lambda x: "NA" if x != x else str(int(round(x)))


def main():
    longr, wide, miss = [], {t: {} for t in P6}, []
    for ls, d in LS:
        D = f"{R}/work/final_v3/{d}"
        if not os.path.exists(f"{D}/CONTROLS.DONE"):
            miss.append(ls); continue
        cs = {r["trait"]: r for r in rows(f"{D}/controls_summary.tsv")}
        cm = {(r["trait"], r["model"]): r for r in rows(f"{D}/clade_meta.tsv")}
        for t in P6:
            r, w = cs.get(t, {}), wide[t]
            for m, term in MODELS:
                v = {k: fnum(r.get(f"{m}_{k}")) for k in ("est", "se", "p", "q", "n")}
                longr.append([t, ls, "Animalia", m, term] + [v[k] for k in ("est", "se", "p", "q", "n")] + [NA, ""])
                if m in WIDE:
                    for k in ("est", "p", "n"):
                        w[f"{ls}_{m}_{k}"] = v[k]
            for m, term in META:
                q = cm.get((t, m), {})
                v = {k: fnum(q.get(k)) for k in ("k", "est", "se", "p", "q", "I2")}
                longr.append([t, ls, "clade_meta", m, term, v["est"], v["se"], v["p"], v["q"], v["k"], v["I2"],
                              q.get("method", "")])
                for k in ("k", "est", "p", "I2"):
                    w[f"{ls}_meta_{m}_{k}"] = v[k]
            a, b = w[f"{ls}_all_est"], w[f"{ls}_deepdose_est"]
            w[f"{ls}_dose_same_sign"] = NA if a != a or b != b else int((a > 0) == (b > 0))
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/dose_series.tsv", "w") as f:
        f.write("trait\tlabelset\tscope\tmodel\tterm\test\tse\tp\tq\tn_or_k\tI2\tmeta_method\n")
        for x in longr:
            f.write("\t".join(s(v) for v in x) + "\n")
    cols = []
    for ls, _ in LS:
        cols += [f"{ls}_{m}_{k}" for m in WIDE for k in ("est", "p", "n")]
        cols += [f"{ls}_meta_{m}_{k}" for m, _ in META for k in ("k", "est", "p", "I2")] + [f"{ls}_dose_same_sign"]
    with open(f"{OUT}/dose_wide.tsv", "w") as f:
        f.write("trait\t" + "\t".join(cols) + "\n")
        for t in P6:
            f.write(t + "\t" + "\t".join(s(wide[t].get(c, NA)) for c in cols) + "\n")
    with open(f"{OUT}/report.txt", "w") as f:
        f.write("T1 depth dose-response (44 outputs, no new fit); est in trait SD units\n")
        if miss:
            f.write(f"missing label sets (no CONTROLS.DONE): {', '.join(miss)}\n")
        for t in P6:
            w = wide[t]; f.write(f"{t}\n")
            for ls, _ in LS:
                if ls in miss:
                    continue
                f.write(f"  {ls:10s} deep {g(w[f'{ls}_all_est'])} (p {g(w[f'{ls}_all_p'], 2)}, n {gi(w[f'{ls}_all_n'])})"
                        f" | log10 depth + Thab {g(w[f'{ls}_ldepthT_est'])} (p {g(w[f'{ls}_ldepthT_p'], 2)})"
                        f" | within deep {g(w[f'{ls}_deepdose_est'])} (p {g(w[f'{ls}_deepdose_p'], 2)},"
                        f" n {gi(w[f'{ls}_deepdose_n'])}; same sign {s(w[f'{ls}_dose_same_sign'])})"
                        f" | clade meta within deep k {gi(w[f'{ls}_meta_deepdose_k'])}"
                        f" est {g(w[f'{ls}_meta_deepdose_est'])} p {g(w[f'{ls}_meta_deepdose_p'], 2)}"
                        f" I2 {g(w[f'{ls}_meta_deepdose_I2'], 3)}\n")
    open(f"{OUT}/T1.DONE", "w").write("ok\n")
    print(open(f"{OUT}/report.txt").read())


if __name__ == "__main__":
    main()
