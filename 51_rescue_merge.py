#!/usr/bin/env python3
"""RELAX rescue: missingness check and non-destructive merge (2026-10-25).
HyPhy 2.5.96 RELAX wrote empty JSONs for a subset of runs (crash in the reporting-only ancestral 'substitutions'
block after the alternative fit); these were re-run with scripts/hyphy_patch RELAX_noanc.bf into
work/codon/<C>/sel_rescue (same inputs, same options).
Check (per clade, written to OUT/rescue_check_<C>.tsv + summary): empty rate of the original RELAX runs by set
(real vs nullK) and by taxon count; for the rescued runs vs the originally successful runs, p < 0.05 rate, median K,
Mann-Whitney U p on p-values and on log K, separately for real and null sets.
Merge (--apply): a rescued JSON is copied to sel/ only when the sel/ JSON is missing / empty / without a RELAX
p-value and the rescued JSON has one; the replaced sel/ file is moved to sel_pre_rescue/; manifest
sel_rescue/MERGED.tsv (file, md5 old, md5 new, p, K).
Usage: python 51_rescue_merge.py "CLADES" [--apply] [--out work/rescue]"""
import glob, hashlib, json, math, os, re, shutil, sys
from statistics import median
from scipy.stats import mannwhitneyu

R = os.path.expanduser("~/deepsea_conv")
A = sys.argv
CL = A[1].split()
APPLY = "--apply" in A
OUT = A[A.index("--out") + 1] if "--out" in A else f"{R}/work/rescue"
os.makedirs(OUT, exist_ok=True)


def pk(fn):
    try:
        if os.path.getsize(fn) == 0:
            return None
        j = json.load(open(fn))
        tr = j.get("test results", {})
        p = tr.get("p-value")
        K = tr.get("relaxation or intensification parameter")
        return (float(p), float(K)) if p is not None else None
    except Exception:
        return None


def md5(fn):
    return hashlib.md5(open(fn, "rb").read()).hexdigest() if os.path.exists(fn) else "-"


def ntip(W, base):
    fn = f"{W}/trees/{base}.hyphy.nwk"
    if not os.path.exists(fn):
        return -1
    s = re.sub(r"\[[^\]]*\]|\{[^}]*\}", "", open(fn).read())
    return len(re.findall(r"[(,]\s*([^(),:;\s]+)", s))


def stats(v):
    ps = [x[0] for x in v]; ks = [x[1] for x in v]
    return len(v), (sum(p < 0.05 for p in ps) / len(v) if v else float("nan")), (median(ks) if ks else float("nan"))


def mw(a, b):
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    return mannwhitneyu(a, b).pvalue


summ = open(f"{OUT}/rescue_check_summary.tsv", "w")
summ.write("clade\tset\tn_orig\tn_empty\tempty_rate\tn_ok\tok_p05\tok_medK\tn_rescued\tres_p05\tres_medK\t"
           "MW_p_pvalue\tMW_p_logK\tntip_med_ok\tntip_med_empty\tMW_p_ntip\tn_unrescued\n")
for C in CL:
    W = f"{R}/work/codon/{C}"
    fs = sorted(glob.glob(f"{W}/sel/*.relax.json"))
    rows, groups = [], {}
    for fn in fs:
        b = os.path.basename(fn)[:-len(".relax.json")]
        st = "null" if re.search(r"\.null\d+$", b) else "real"
        o = pk(fn)
        rf = f"{W}/sel_rescue/{b}.relax.json"
        r = pk(rf) if o is None and os.path.exists(rf) else None
        nt = ntip(W, b)
        rows.append((b, st, "ok" if o else ("rescued" if r else "empty"), o or r, nt))
    with open(f"{OUT}/rescue_check_{C}.tsv", "w") as f:
        f.write("run\tset\tstatus\tp\tK\tn_tips\n")
        for b, st, s, v, nt in rows:
            f.write(f"{b}\t{st}\t{s}\t{v[0] if v else 'NA'}\t{v[1] if v else 'NA'}\t{nt}\n")
    for st in ("real", "null"):
        x = [r for r in rows if r[1] == st]
        ok = [r[3] for r in x if r[2] == "ok"]; rs = [r[3] for r in x if r[2] == "rescued"]
        ne = sum(r[2] != "ok" for r in x)
        nok, okp, okk = stats(ok); nrs, rsp, rsk = stats(rs)
        lk = lambda v: [math.log(max(k, 1e-12)) for _, k in v]
        tok = [r[4] for r in x if r[2] == "ok" and r[4] > 0]; tem = [r[4] for r in x if r[2] != "ok" and r[4] > 0]
        summ.write(f"{C}\t{st}\t{len(x)}\t{ne}\t{ne / len(x) if x else float('nan'):.4f}\t{nok}\t{okp:.4f}\t{okk:.4g}\t"
                   f"{nrs}\t{rsp:.4f}\t{rsk:.4g}\t{mw([p for p, _ in ok], [p for p, _ in rs]):.4g}\t"
                   f"{mw(lk(ok), lk(rs)):.4g}\t{median(tok) if tok else 'NA'}\t{median(tem) if tem else 'NA'}\t"
                   f"{mw(tok, tem):.4g}\t{sum(r[2] == 'empty' for r in x)}\n")
    if APPLY:
        bk = f"{W}/sel_pre_rescue"; os.makedirs(bk, exist_ok=True)
        man = open(f"{W}/sel_rescue/MERGED.tsv", "a")
        n = 0
        for b, st, s, v, nt in rows:
            if s != "rescued":
                continue
            dst, src = f"{W}/sel/{b}.relax.json", f"{W}/sel_rescue/{b}.relax.json"
            m_old = md5(dst)
            if os.path.exists(dst):
                shutil.move(dst, f"{bk}/{b}.relax.json")
            shutil.copy2(src, dst)
            man.write(f"{b}.relax.json\t{m_old}\t{md5(dst)}\t{v[0]}\t{v[1]}\n"); n += 1
        man.close()
        print(f"{C}: merged {n}", flush=True)
summ.close()
print(open(f"{OUT}/rescue_check_summary.tsv").read())
