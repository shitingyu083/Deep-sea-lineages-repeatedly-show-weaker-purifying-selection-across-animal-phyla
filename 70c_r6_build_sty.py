#!/usr/bin/env python3
"""R6-1 step c (PREREG_r6.md): paired original-route (O) and uniform-route (U) codon alignments for the R6 sample.
1. Gene mapping, per re-annotated species (REANNOT.DONE in work/addon/R6/reannot or work/addon/A1_reannot): the
   production HOG protein g is mapped to the miniprot model m by diamond reciprocal best hit (blastp --sensitive;
   best bitscore each way; m's best hit in the full production proteome must be g; pident >= 90, qcovhsp >= 50).
2. Per sampled HOG: species of the production alignment; miniprot-route species enter both arms with their
   production CDS; annotation-route species enter O with the production gene and U with the mapped model (dropped
   from both arms if not re-annotated or not mapped); other routes (tsa ...) are dropped from both arms.
   Each arm then goes through the 22_codon_aln.py steps (CDS must translate to the protein; MAFFT --auto; trimAl
   -automated1 -backtrans; drop > 50 % gap rows). The species kept in both arms are intersected and both alignments
   restricted to them; the HOG is kept if each arm has >= 4 species, >= 2 deep, >= 2 shallow and >= 300 nt, and
   the production real and null1 trees, pruned to these species (ete3, branch lengths kept), each still have
   >= 2 {FG} and >= 2 {BG} tips.
Writes OUT/aln/<clade>/<hog>.{O,U}.fas, OUT/trees/<clade>/<hog>.{real,null1}.nwk, OUT/build_report.tsv,
OUT/map_report.tsv; OUT/BUILD.DONE. Production files are only read.
Usage: envs/dc/bin/python 70c_r6_build_sty.py [NPROC 32]"""
import csv, os, re, subprocess, sys, tempfile
from multiprocessing import Pool
from Bio.Seq import Seq
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv"); E = f"{R}/envs/dc/bin"
OUT = f"{R}/work/addon/R6"
RE_DIRS = [f"{OUT}/reannot", f"{R}/work/addon/A1_reannot"]
ANNOT = ("ncbi_annot", "local", "extra_local")
PID, QCOV = 90.0, 50.0
sys.path.insert(0, f"{R}/scripts")
from importlib import import_module
lab = import_module("21_hog_select").labels()
src = {r["code"]: r["source"] for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}


def read_fa(path, keep=None):
    d, k = {}, None
    for l in open(path):
        if l.startswith(">"):
            k = l[1:].split()[0]; k = k if (keep is None or k in keep) else None
            if k: d[k] = []
        elif k:
            d[k].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}


def reannot_dir(s):
    for d in RE_DIRS:
        if os.path.exists(f"{d}/{s}/REANNOT.DONE"):
            return f"{d}/{s}"
    return None


def best(path):
    b = {}
    for l in open(path):
        q, s, pid, ln, ql, sl, bits, qc = l.rstrip("\n").split("\t")
        if q not in b or float(bits) > b[q][1]:
            b[q] = (s, float(bits), float(pid), float(qc))
    return b


def dmnd(q, db_fa, out, td):
    subprocess.run([f"{E}/diamond", "makedb", "--in", db_fa, "-d", f"{td}/db", "--quiet"], check=True)
    subprocess.run([f"{E}/diamond", "blastp", "-q", q, "-d", f"{td}/db", "-o", out, "--sensitive", "--quiet",
                    "--max-target-seqs", "5", "--threads", "4", "--outfmt", "6", "qseqid", "sseqid", "pident",
                    "length", "qlen", "slen", "bitscore", "qcovhsp"], check=True)
    return best(out)


def map_species(args):
    s, genes = args
    d = reannot_dir(s)
    if not d:
        return s, {}, len(genes), "no_reannot"
    with tempfile.TemporaryDirectory(dir=f"{OUT}/tmp") as td:
        prod = read_fa(f"{R}/data/proteomes/{s}.faa")
        with open(f"{td}/q.faa", "w") as f:
            f.writelines(f">{g}\n{prod[g]}\n" for g in genes if g in prod)
        fw = dmnd(f"{td}/q.faa", f"{d}/prot.faa", f"{td}/fw.tsv", td)
        mp = read_fa(f"{d}/prot.faa", {v[0] for v in fw.values()})
        with open(f"{td}/r.faa", "w") as f:
            f.writelines(f">{m}\n{x}\n" for m, x in mp.items())
        with open(f"{td}/prod.faa", "w") as f:
            f.writelines(f">{g}\n{x}\n" for g, x in prod.items())
        os.makedirs(f"{td}/r", exist_ok=True)
        rv = dmnd(f"{td}/r.faa", f"{td}/prod.faa", f"{td}/rv.tsv", f"{td}/r")
    m = {g: v[0] for g, v in fw.items() if v[2] >= PID and v[3] >= QCOV and rv.get(v[0], ("",))[0] == g}
    return s, m, len(genes), "ok"


SEQ = {}  # (arm, species) -> (prot dict, cds dict) of the genes used, preloaded before the pool forks


def valid(nt, q):  # = 22_codon_aln.py one(): CDS must translate to its protein
    nt = nt.upper()
    if len(nt) % 3:
        return None
    aa = str(Seq(nt).translate())
    if aa.endswith("*"):
        aa, nt = aa[:-1], nt[:-3]
    q = q.rstrip("*")
    if "*" in aa or len(aa) != len(q) or sum(a != b for a, b in zip(aa[1:], q[1:])) > 0.02 * len(q):
        return None
    return aa, nt


def align(prot, cds, td, tag):  # = 22_codon_aln.py: MAFFT --auto, trimAl -automated1 -backtrans, <= 50 % gaps
    with open(f"{td}/{tag}p.fa", "w") as f:
        f.writelines(f">{s}\n{x}\n" for s, x in prot.items())
    with open(f"{td}/{tag}c.fa", "w") as f:
        f.writelines(f">{s}\n{x}\n" for s, x in cds.items())
    with open(f"{td}/{tag}p.aln", "w") as f:
        if subprocess.run([f"{E}/mafft", "--auto", "--thread", "1", "--quiet", f"{td}/{tag}p.fa"], stdout=f).returncode:
            return None
    r = subprocess.run([f"{E}/trimal", "-in", f"{td}/{tag}p.aln", "-backtrans", f"{td}/{tag}c.fa",
                        "-automated1", "-out", f"{td}/{tag}c.aln"], capture_output=True)
    if r.returncode or not os.path.exists(f"{td}/{tag}c.aln"):
        return None
    aln = read_fa(f"{td}/{tag}c.aln")
    L = len(next(iter(aln.values()))) if aln else 0
    return {s: x for s, x in aln.items() if L and x.count("-") / L <= 0.5}


def counts(keys):
    return sum(lab.get(s) == "D" for s in keys), sum(lab.get(s) == "S" for s in keys)


def prune(path, keep):
    t = Tree(open(path).read().strip(), format=1)
    names = {re.sub(r"\{.*\}$", "", l.name): l.name for l in t.get_leaves()}
    t.prune([names[s] for s in keep if s in names], preserve_branch_length=True)
    tips = [l.name for l in t.get_leaves()]
    return t.write(format=1), sum(n.endswith("{FG}") for n in tips), sum(n.endswith("{BG}") for n in tips)


def one(args):
    clade, hog, members = args
    W, A, T = f"{R}/work/codon/{clade}", f"{OUT}/aln/{clade}", f"{OUT}/trees/{clade}"
    sp = [l[1:].split()[0] for l in open(f"{W}/aln/{hog}.fas") if l.startswith(">")]
    gene = dict(m.split(":", 1) for m in reversed(members.split(",")))  # first gene per species, as 22
    arm = {"O": {}, "U": {}}
    drop = {"other": 0, "unmapped": 0, "invalid": 0}
    for s in sp:
        g = gene.get(s)
        if src.get(s) == "miniprot":
            go, gu = g, g
        elif src.get(s) in ANNOT:
            go, gu = g, MAP.get(s, {}).get(g)
            if gu is None:
                drop["unmapped"] += 1; continue
        else:
            drop["other"] += 1; continue
        po, co = SEQ.get(("O", s), ({}, {}))
        pu, cu = SEQ.get(("U", s), ({}, {})) if src.get(s) in ANNOT else (po, co)
        vo = valid(co.get(go, ""), po.get(go, "")) if go in po else None
        vu = valid(cu.get(gu, ""), pu.get(gu, "")) if gu in pu else None
        if not vo or not vu:
            drop["invalid"] += 1; continue
        arm["O"][s], arm["U"][s] = vo, vu
    base = [clade, hog, len(sp), drop["other"], drop["unmapped"], drop["invalid"]]
    nD, nS = counts(arm["O"])
    if len(arm["O"]) < 4 or nD < 2 or nS < 2:
        return base + ["too_few_valid_cds", len(arm["O"]), nD, nS, 0, 0, 0, 0]
    with tempfile.TemporaryDirectory(dir=f"{OUT}/tmp") as td:
        al = {}
        for a in "OU":
            al[a] = align({s: v[0] for s, v in arm[a].items()}, {s: v[1] for s, v in arm[a].items()}, td, a)
            if al[a] is None:
                return base + [f"align_fail_{a}", len(arm[a]), nD, nS, 0, 0, 0, 0]
    keep = sorted(set(al["O"]) & set(al["U"]))
    nD, nS = counts(keep)
    LO, LU = (len(next(iter(al[a].values()))) if al[a] else 0 for a in "OU")
    if len(keep) < 4 or nD < 2 or nS < 2 or min(LO, LU) < 300:
        return base + ["post_trim_filter", len(keep), nD, nS, LO, LU, 0, 0]
    trees, fgbg = {}, []
    for st, f in (("real", f"{W}/trees/{hog}.hyphy.nwk"), ("null1", f"{W}/trees/{hog}.null1.hyphy.nwk")):
        trees[st], nf, nb = prune(f, keep)
        fgbg += [nf, nb]
    if min(fgbg) < 2:
        return base + ["tree_fg_bg", len(keep), nD, nS, LO, LU] + fgbg[:2]
    os.makedirs(A, exist_ok=True); os.makedirs(T, exist_ok=True)
    for a in "OU":
        with open(f"{A}/{hog}.{a}.fas.tmp", "w") as f:
            f.writelines(f">{s}\n{al[a][s]}\n" for s in keep)
        os.replace(f"{A}/{hog}.{a}.fas.tmp", f"{A}/{hog}.{a}.fas")
    for st, t in trees.items():
        with open(f"{T}/{hog}.{st}.nwk", "w") as f:
            f.write(t + "\n")
    return base + ["ok", len(keep), nD, nS, LO, LU] + fgbg[:2]


MAP = {}


def main():
    npr = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    os.makedirs(f"{OUT}/tmp", exist_ok=True)
    samp = list(csv.DictReader(open(f"{OUT}/sample.tsv"), delimiter="\t"))
    mem = {}
    for c in sorted({r["clade"] for r in samp}):
        hc = "Actinopteri" if c == "Actinopteri_m" else c  # Actinopteri_m = relabelled Actinopteri HOGs (DESIGN.md)
        for l in list(open(f"{R}/work/codon/{hc}/hogs.tsv"))[1:]:
            p = l.rstrip("\n").split("\t"); mem[(c, p[0])] = p[4]
    jobs = [(r["clade"], r["hog"], mem[(r["clade"], r["hog"])]) for r in samp]
    need = {}
    for c, h, m in jobs:
        sp = {l[1:].split()[0] for l in open(f"{R}/work/codon/{c}/aln/{h}.fas") if l.startswith(">")}
        g = {}
        for x in reversed(m.split(",")):
            s, gg = x.split(":", 1); g[s] = gg
        for s in sp:
            if s in g:
                need.setdefault(s, set()).add(g[s])
    annot = sorted(s for s in need if src.get(s) in ANNOT)
    with Pool(min(npr, 16)) as pool:
        res = pool.map(map_species, [(s, sorted(need[s])) for s in annot])
    with open(f"{OUT}/map_report.tsv", "w") as f:
        f.write("code\tstatus\tn_query\tn_rbh\n")
        for s, m, n, st in res:
            MAP[s] = m; f.write(f"{s}\t{st}\t{n}\t{len(m)}\n")
    for s, gs in need.items():
        SEQ[("O", s)] = (read_fa(f"{R}/data/proteomes/{s}.faa", gs), read_fa(f"{R}/data/cds/{s}.fna", gs))
        if s in MAP and MAP[s]:
            d, mm = reannot_dir(s), {MAP[s][g] for g in gs if g in MAP[s]}
            SEQ[("U", s)] = (read_fa(f"{d}/prot.faa", mm), read_fa(f"{d}/cds.fna", mm))
    with Pool(npr) as pool, open(f"{OUT}/build_report.tsv", "w") as rep:
        rep.write("clade\thog\tn_prod\tdrop_other\tdrop_unmapped\tdrop_invalid\tstatus\tn\tnD\tnS\tL_O\tL_U"
                  "\tfg_real\tbg_real\n")
        for r in pool.imap_unordered(one, jobs, chunksize=2):
            rep.write("\t".join(map(str, r)) + "\n")
    open(f"{OUT}/BUILD.DONE", "w").close()
    print("build done", len(jobs))


if __name__ == "__main__":
    main()
