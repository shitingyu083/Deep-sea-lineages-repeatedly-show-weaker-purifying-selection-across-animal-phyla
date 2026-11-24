#!/usr/bin/env python
"""R5 prep (PREREG_r5.md): RELAX job list for D1 (null2 on nullhy HOGs) and D3 (route-pure pruned real/null1), and the
D2a simulation job list. RELAX job line = OUT ALN TREE LOG. Existing results are reused (symlinked), never recomputed.
Usage: python 69a_r5_prep_sty.py"""
import os, csv, random, importlib.util
B = os.path.expanduser("~/deepsea_conv"); O = f"{B}/work/addon/R5"; PY = f"{B}/envs/dc/bin/python"
sp = importlib.util.spec_from_file_location("ps", f"{B}/scripts/63a3_powersim_sty.py"); ps = importlib.util.module_from_spec(sp); sp.loader.exec_module(ps)
CL = ["Actinopteri_m", "Annelida", "Arthropoda", "Chondrichthyes", "Cnidaria", "Echinodermata", "Mollusca", "Porifera"]
ROUTE = {"ncbi_annot": "A", "local": "A", "extra_local": "A", "miniprot": "M"}
src = {r["code"]: r["source"] for r in csv.DictReader(open(f"{B}/work/physchem/species_metrics.tsv"), delimiter="\t")}
FEAT = list(csv.DictReader(open(f"{B}/work/addon/A4/sel_features.tsv"), delimiter="\t"))


def oldjson(W, stem):
    return next((p for p in (f"{W}/sel/{stem}.relax.json", f"{W}/sel_rescue/{stem}.relax.json")
                 if os.path.exists(p) and os.path.getsize(p) > 0), None)


def wr(n):
    s = ("(" + ",".join(wr(c) for c in n.ch) + ")") if n.ch else n.name
    if n.tag: s += "{" + n.tag + "}"
    return s + (f":{n.bl:.8g}" if n.up is not None else "")


def prune(n, keep):
    """Pruned copy of subtree n (None if no kept tip); a unary node passes its length to its child, child tag kept."""
    if not n.ch:
        if n.name not in keep: return None
        m = ps.N(); m.name, m.tag, m.bl = n.name, n.tag, n.bl; return m
    kids = [k for k in (prune(c, keep) for c in n.ch) if k is not None]
    if not kids: return None
    if len(kids) == 1:
        k = kids[0]; k.bl += n.bl; return k
    m = ps.N(); m.name, m.tag, m.bl = n.name, n.tag, n.bl; m.ch = kids
    for k in kids: k.up = m
    return m


def tiptags(t):
    return {x.name: x.tag for x in ps.nodes(t) if not x.ch}


jobs, stats = [], []
os.makedirs(O, exist_ok=True)
for C in CL:
    W = f"{B}/work/codon/{C}"; WA = W if os.path.isdir(f"{W}/aln") else f"{B}/work/codon/Actinopteri"
    st = {(r["hog"], r["set"]) for r in FEAT if r["clade"] == C}
    hy = sorted({r["hog"] for r in FEAT if r["clade"] == C and r["set"] == "real" and r["in_nullhy"] == "1"
                 and (r["hog"], "null1") in st})
    # ---- D1: null2 RELAX on the nullhy HOGs
    d1 = f"{O}/null2/{C}"; os.makedirs(d1, exist_ok=True); n1 = nre = 0
    for h in hy:
        tr, fg, aln = f"{W}/trees/{h}.null2.hyphy.nwk", f"{W}/trees/{h}.null2.fg.txt", f"{WA}/aln/{h}.fas"
        if not all(os.path.exists(p) and os.path.getsize(p) > 0 for p in (tr, fg, aln)): continue
        old = oldjson(W, f"{h}.null2"); out = f"{d1}/{h}.null2.relax.json"
        if old and not os.path.exists(out): os.symlink(old, out)
        nre += bool(old); n1 += 1
        if not os.path.exists(out): jobs.append(("D1", out, aln, tr, f"{d1}/{h}.null2.relax.log"))
    stats.append(f"D1\t{C}\tnullhy={len(hy)}\twith_null2={n1}\treused={nre}")
    # ---- D3: route-pure pruning (A = annotation, M = miniprot; TSA dropped)
    for rt in ("A", "M"):
        el = []
        for h in hy:
            ok, ks = True, {}
            for s in ("real", "null1"):
                p = f"{W}/trees/{h}{'' if s == 'real' else '.null1'}.hyphy.nwk"
                if not os.path.exists(p): ok = False; break
                tt = tiptags(ps.parse(open(p).read()))
                keep = {k for k in tt if ROUTE.get(src.get(k)) == rt}
                nf = sum(tt[k] == "FG" for k in keep); nb = sum(tt[k] == "BG" for k in keep)
                if nf < 2 or nb < 2: ok = False; break
                ks[s] = (keep, set(tt))
            if ok: el.append((h, ks))
        pick = sorted(random.Random(20261030).sample(el, min(200, len(el))), key=lambda x: x[0]) if el else []
        d3 = f"{O}/route/{C}/{rt}"; os.makedirs(d3, exist_ok=True); nre = 0
        with open(f"{d3}/picked.tsv", "w") as fo:
            fo.write("hog\tset\tn_keep\tn_all\treused\n")
            for h, kk in pick:
                aln = ps.read_fas(f"{WA}/aln/{h}.fas")
                for s in ("real", "null1"):
                    keep, allt = kk[s]; tag = "" if s == "real" else ".null1"
                    out = f"{d3}/{h}.{s}.relax.json"; reuse = keep == allt
                    if reuse:
                        old = oldjson(W, f"{h}{tag}")
                        if old and not os.path.exists(out): os.symlink(old, out)
                        nre += 1
                    else:
                        t = prune(ps.parse(open(f"{W}/trees/{h}{tag}.hyphy.nwk").read()), keep)
                        open(f"{d3}/{h}.{s}.nwk", "w").write(wr(t) + ";\n")
                        L = max(len(aln[k]) for k in keep)  # drop codon columns that are gap-only among kept tips
                        keepcod = [i for i in range(0, L, 3) if any(aln[k][i:i + 3].strip("-N?") for k in keep)]
                        with open(f"{d3}/{h}.{s}.fas", "w") as fa:
                            for k in sorted(keep):
                                fa.write(f">{k}\n" + "".join(aln[k][i:i + 3].ljust(3, "-") for i in keepcod) + "\n")
                        if not os.path.exists(out):
                            jobs.append(("D3", out, f"{d3}/{h}.{s}.fas", f"{d3}/{h}.{s}.nwk", f"{d3}/{h}.{s}.relax.log"))
                    fo.write(f"{h}\t{s}\t{len(keep)}\t{len(allt)}\t{int(reuse)}\n")
        stats.append(f"D3\t{C}\t{rt}\teligible={len(el)}\tpicked={len(pick)}\treused_rows={nre}")
random.Random(20261030).shuffle(jobs)
with open(f"{O}/relax_jobs.txt", "w") as fo:
    for j in jobs: fo.write(" ".join(j[1:]) + "\n")
# ---- D2a: new scenarios on the 160 R2 HOGs; R2 calibrations reused via _fit symlink
R2 = f"{B}/work/addon/R2/main"; n2 = 0
with open(f"{O}/sim_jobs.txt", "w") as fo:
    for r in csv.DictReader(open(f"{R2}/hogs.tsv"), delimiter="\t"):
        C, h = r["clade"], r["hog"]; d = f"{O}/power/{C}/{h}"; os.makedirs(d, exist_ok=True)
        if not os.path.exists(f"{R2}/{C}/{h}/_fit/scale.txt"): stats.append(f"D2\tno_scale\t{C}\t{h}"); continue
        if not os.path.exists(f"{d}/_fit"): os.symlink(f"{R2}/{C}/{h}/_fit", f"{d}/_fit")
        for sc in ("K85", "K90", "K75h"):
            if not os.path.exists(f"{d}/{sc}/res.tsv"):
                fo.write(f"{PY} {B}/scripts/69c_r5_powersim_sty.py {C} {h} {sc} 0 {d}/{sc}\n"); n2 += 1
stats.append(f"JOBS\trelax_D1={sum(j[0] == 'D1' for j in jobs)}\trelax_D3={sum(j[0] == 'D3' for j in jobs)}\tsim={n2}")
open(f"{O}/prep_stats.tsv", "w").write("\n".join(stats) + "\n"); print("\n".join(stats))
