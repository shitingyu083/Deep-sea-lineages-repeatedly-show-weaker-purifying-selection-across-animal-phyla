#!/usr/bin/env python3
# R2 (PREREG_nature.md): CAFE5 deep-lambda rank test over 40 tip permutations per clade.
# gain = negLnL(base_e) - negLnL(model). Null gains: nullp1-10 (work/cafe/<C>) + nullp11-40 (work/addon/R2/cafe_perm/<C>).
# rank p = (1 + #null >= deep) / (1 + K). Fits with gain < 0 (worse than the 1-lambda base, i.e. failed optimisation)
# are reported and also excluded in a sensitivity column. Output work/addon/R2/cafe_perm/cafe_rank40.tsv.
import os, re, glob
R = os.path.expanduser("~/deepsea_conv"); O = f"{R}/work/addon/R2/cafe_perm"

def nll(d):
    f = f"{d}/Base_results.txt"
    if not os.path.isfile(f) or os.path.getsize(f) == 0: return None, None
    t = open(f).read()
    m = re.search(r"\(-lnL\):\s*([0-9.eE+-]+)", t); l = re.search(r"Lambda:\s*([^\n]+)", t)
    lam = [float(x) for x in l.group(1).replace(",", " ").split()] if l else []
    return (float(m.group(1)) if m else None), lam

rows = ["clade\tsrc\tbase_nll\tdeep_gain\tlam_bg\tlam_deep\tlam_ratio\tK\tn_ge\trank_p\tK_ok\tn_ge_ok\trank_p_ok\tn_failed\tnull_max\tnull_median"]
for P in sorted(d for d in glob.glob(f"{O}/*") if os.path.isdir(d)):
    C = os.path.basename(P); D = f"{R}/work/cafe/{C}"
    if not os.path.isdir(D): continue
    b, _ = nll(f"{D}/base_e"); dn, lam = nll(f"{D}/deep")
    if b is None or dn is None: rows.append(f"{C}\tNA"); continue
    g = b - dn; nulls = []
    for k in range(1, 41):
        src = D if k <= 10 else P
        v, _ = nll(f"{src}/nullp{k}")
        if v is not None: nulls.append(b - v)
    ok = [x for x in nulls if x >= 0]
    K, ge = len(nulls), sum(x >= g for x in nulls); Ko, geo = len(ok), sum(x >= g for x in ok)
    s = sorted(nulls); med = s[len(s) // 2] if s else float("nan")
    lb, ld = (lam + [float("nan")] * 2)[:2]
    rows.append(f"{C}\tperm40\t{b:.1f}\t{g:.1f}\t{lb:.3g}\t{ld:.3g}\t{ld / lb:.3f}\t{K}\t{ge}\t{(1 + ge) / (1 + K):.4f}\t{Ko}\t{geo}\t"
                f"{(1 + geo) / (1 + Ko):.4f}\t{K - Ko}\t{max(nulls):.1f}\t{med:.1f}")
open(f"{O}/cafe_rank40.tsv", "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
open(f"{O}/COLLECT.DONE", "w").close()
