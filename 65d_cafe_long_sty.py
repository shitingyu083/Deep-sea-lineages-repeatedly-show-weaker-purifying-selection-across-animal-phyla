#!/usr/bin/env python3
# Long-format CAFE5 gains for Fig. 4d: gain = negLnL(base_e) - negLnL(model); models deep, null1-5 (matched shallow),
# nullp1-40 (tip permutations; 11-40 from work/addon/R2/cafe_perm). Read-only on inputs; writes cafe_gains_long.tsv.
import os, re, glob
R = os.path.expanduser("~/deepsea_conv"); O = f"{R}/work/addon/R2/cafe_perm"

def nll(d):
    f = f"{d}/Base_results.txt"
    if not os.path.isfile(f) or os.path.getsize(f) == 0: return None
    m = re.search(r"\(-lnL\):\s*([0-9.eE+-]+)", open(f).read())
    return float(m.group(1)) if m else None

out = ["clade\tsrc\tmodel\tgain"]
cl = sorted(os.path.basename(d) for d in glob.glob(f"{O}/*") if os.path.isdir(d) and os.path.isdir(f"{R}/work/cafe/{os.path.basename(d)}"))
for C in cl + ["Porifera"]:
    D = f"{R}/work/cafe/{C}"; P = f"{O}/{C}"; b = nll(f"{D}/base_e")
    src = "perm40" if C in cl else "main10"
    if b is None:
        out.append(f"{C}\t{src}\tbase_e\tNA"); continue
    models = ["deep"] + [f"null{k}" for k in range(1, 6)] + [f"nullp{k}" for k in range(1, 41)]
    for m in models:
        k = int(m[5:]) if m.startswith("nullp") else 0
        d = f"{P}/{m}" if k > 10 else f"{D}/{m}"
        v = nll(d)
        if v is None and not os.path.isdir(d): continue
        out.append(f"{C}\t{src}\t{m}\t" + ("NA" if v is None else f"{b - v:.2f}"))
open(f"{O}/cafe_gains_long.tsv", "w").write("\n".join(out) + "\n")
print(len(out) - 1, "rows")
