#!/usr/bin/env python
"""R1 (PREREG_nature.md): per-tip lineage covariates and per-row FG-vs-BG covariate differences for the A4 HyPhy frame.
Per tip (per clade dir): omega_tip = sum N_sub / sum S_sub over CSUBST terminal branches of all real HOGs (genome-wide
terminal dN/dS, Ne proxy); CDS GC and GC3 over all real codon alignments; genome GC; BUSCO C (prot_C); dated terminal
branch length (LSD2); depth_median, T_hab.  Per row (sel_features real/null1): FG/BG tips from the row's hyphy.nwk tags.
Outputs work/addon/R1/tip_cov.tsv, row_cov.tsv."""
import os, re, csv, glob, math, sys, importlib.util
B = os.path.expanduser("~/deepsea_conv"); OUT = f"{B}/work/addon/R1"; os.makedirs(OUT, exist_ok=True)
sp = importlib.util.spec_from_file_location("ps", f"{B}/scripts/63a_powersim_sty.py"); ps = importlib.util.module_from_spec(sp)
sp.loader.exec_module(ps)
CD = ["Actinopteri", "Annelida", "Arthropoda", "Chondrichthyes", "Cnidaria", "Echinodermata", "Mollusca", "Porifera"]
def rd(p): return list(csv.DictReader(open(p), delimiter="\t"))
def fl(x):
    try: return float(x)
    except Exception: return float("nan")
qc = {r["code"]: r for r in rd(f"{B}/meta/qc_table.tsv")}
gg = {r["code"]: fl(r["genome_gc"]) for r in rd(f"{B}/work/addon/genome_gc.tsv")}
hb = {r["code"]: r for r in rd(f"{B}/meta/habitat_temp.tsv")}
tt = ps.parse(open(f"{B}/work/final_v3/date/lsd.timetree.nwk").read())
tbl = {n.name: n.bl for n in ps.nodes(tt) if not n.ch}
tip = {}
for c in CD:
    W = f"{B}/work/codon/{c}"; NS = {}; nh = 0
    for f in glob.glob(f"{W}/cs/*/csubst_search/csubst_b.tsv"):
        nh += 1
        for r in csv.DictReader(open(f), delimiter="\t"):
            b = r["branch_name"]
            if b in qc:
                a = NS.setdefault(b, [0.0, 0.0, 0]); a[0] += fl(r["N_sub"]); a[1] += fl(r["S_sub"]); a[2] += 1
    G = {}
    for f in glob.glob(f"{W}/aln/*.fas"):
        if ".null" in f: continue
        for k, s in ps.read_fas(f).items():
            a = G.setdefault(k, [0, 0, 0, 0])
            for i in range(0, len(s) - 2, 3):
                cd = s[i:i + 3]
                if cd in ps.CODE:
                    a[0] += sum(ch in "GC" for ch in cd); a[1] += 3; a[2] += cd[2] in "GC"; a[3] += 1
    for k in set(NS) | set(G):
        n, s_, m = NS.get(k, [float("nan")] * 2 + [0]); g = G.get(k, [0, 0, 0, 0])
        tip[(c, k)] = dict(clade_dir=c, code=k, N=n, S=s_, n_branch=m, omega_tip=(n / s_ if s_ and s_ > 0 else float("nan")),
                           gc_cds=(g[0] / g[1] if g[1] else float("nan")), gc3=(g[2] / g[3] if g[3] else float("nan")),
                           genome_gc=gg.get(k, float("nan")), buscoC=fl(qc.get(k, {}).get("prot_C")),
                           tbl=tbl.get(k, float("nan")), depth=fl(hb.get(k, {}).get("depth_median")),
                           T_hab=fl(hb.get(k, {}).get("T_hab")), role=qc.get(k, {}).get("role", ""))
    print(c, "cs HOGs", nh, "tips", len(NS), flush=True)
K = ["clade_dir", "code", "role", "N", "S", "n_branch", "omega_tip", "gc_cds", "gc3", "genome_gc", "buscoC", "tbl", "depth", "T_hab"]
with open(f"{OUT}/tip_cov.tsv", "w") as fo:
    fo.write("\t".join(K) + "\n")
    for v in tip.values(): fo.write("\t".join(str(v[k]) for k in K) + "\n")
def mean(v):
    v = [x for x in v if x == x]; return sum(v) / len(v) if v else float("nan")
def lgd(a, b):
    return math.log(a) - math.log(b) if a == a and b == b and a > 0 and b > 0 else float("nan")
def fgbg(tr):
    fg, bg = set(), set()
    def walk(n, inf):
        f = inf or n.tag == "FG"
        if not n.ch:
            (fg if f else bg if n.tag == "BG" else set()).add(n.name)
        for c in n.ch: walk(c, f)
    walk(tr, False); return fg, bg - fg
RK = ["clade", "hog", "set", "n_fg_tips", "n_bg_tips", "d_omega", "d_gc3", "d_gccds", "d_gc", "d_buscoC", "d_tbl",
      "fg_depth", "d_temp", "min_buscoC"]
nr = 0
with open(f"{OUT}/row_cov.tsv", "w") as fo:
    fo.write("\t".join(RK) + "\n")
    for r in csv.DictReader(open(f"{B}/work/addon/A4/sel_features.tsv"), delimiter="\t"):
        if r["set"] not in ("real", "null1"): continue
        cl = r["clade"]; cdr = "Actinopteri" if cl.startswith("Actinopteri") else cl
        sfx = "" if r["set"] == "real" else ".null1"
        p = f"{B}/work/codon/{cl}/trees/{r['hog']}{sfx}.hyphy.nwk"
        if not os.path.exists(p): continue
        fg, bg = fgbg(ps.parse(open(p).read()))
        T = lambda s, k: [tip.get((cdr, x), {}).get(k, float("nan")) for x in s]
        o = dict(clade=cl, hog=r["hog"], set=r["set"], n_fg_tips=len(fg), n_bg_tips=len(bg),
                 d_omega=lgd(mean(T(fg, "omega_tip")), mean(T(bg, "omega_tip"))),
                 d_gc3=mean(T(fg, "gc3")) - mean(T(bg, "gc3")), d_gccds=mean(T(fg, "gc_cds")) - mean(T(bg, "gc_cds")),
                 d_gc=mean(T(fg, "genome_gc")) - mean(T(bg, "genome_gc")),
                 d_buscoC=mean(T(fg, "buscoC")) - mean(T(bg, "buscoC")),
                 d_tbl=lgd(mean(T(fg, "tbl")), mean(T(bg, "tbl"))),
                 fg_depth=(math.log(mean(T(fg, "depth"))) if mean(T(fg, "depth")) > 0 else float("nan")),
                 d_temp=mean(T(fg, "T_hab")) - mean(T(bg, "T_hab")),
                 min_buscoC=min([x for x in T(fg | bg, "buscoC") if x == x] or [float("nan")]))
        fo.write("\t".join(str(o[k]) for k in RK) + "\n"); nr += 1
print("rows", nr)
open(f"{OUT}/R1A.DONE", "w").write("ok\n")
