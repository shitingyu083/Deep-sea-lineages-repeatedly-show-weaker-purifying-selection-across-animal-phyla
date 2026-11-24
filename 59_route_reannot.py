#!/usr/bin/env python3
"""A1 re-annotation control (PREREG_addon.md, A1). Labelled annot-route tree species are re-annotated by the
production miniprot route (logic of 10_miniprot_annot.py: same reference rule with the own species excluded,
miniprot --outn=3, same model filters); the new proteome goes through compleasm protein (metazoa_odb12, private
lineage view as in 14_protqc_loop.sh) and the 30_physchem.py one() metrics (module imported, its R pointed at a
private view of the new files). delta = annot - miniprot per species and trait. Production files are only read;
everything is written under OUT (default work/addon/A1_reannot). A bad genome writes OUT/<code>/REANNOT.FAIL,
never work/integrity/redownload.txt or data/ncbi/<acc>/.done (the production script's side effects).
Usage
  python 59_route_reannot.py --select  [--out DIR] [--n 40] [--seed 20261025] [--gzmax 1e9] [--force]
  python 59_route_reannot.py --one CODE [--threads 16] [--out DIR] [--force]
  python 59_route_reannot.py --collect [--out DIR]
Exit codes of --one: 0 done, 3 input failure (no genome / reference, md5, gzip, organelle-only), 4 no gene model,
5 compleasm failure or no BUSCO Single hit (REANNOT.FAIL written for 3-5); anything else = crash, retried by
59_run.sh."""
import argparse, csv, glob, gzip, hashlib, importlib.util, json, math, os, platform, random, re, shlex, shutil
import subprocess, sys, time
from collections import Counter
from statistics import mean, stdev

R = os.environ.get("DSC_ROOT") or os.path.expanduser("~/deepsea_conv")
E = f"{R}/envs/dc/bin"
MINCOV, MINLEN = 0.5, 30      # as 10_miniprot_annot.py
MIN_ASM = 1_000_000
ANNOT = ("ncbi_annot", "local", "extra_local")
P6 = ("busco_f_T", "whole_pI", "whole_R_RK", "whole_acid_base", "whole_net_charge7_per_res", "busco_f_F")
GCFAM = ("busco_cds_GC", "busco_GC3", "busco_GC4", "busco_ENC")
LEVELS = ("family", "order", "class", "phylum")
TREE = f"{R}/work/final_v3/date/lsd.timetree.nwk"
q = shlex.quote


def tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def sh(cmd):  # pipefail, as 10_miniprot_annot.py
    subprocess.run(["/bin/bash", "-o", "pipefail", "-c", cmd], check=True)


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def tx(v):  # taxon name or None (empty / NA), for the D-S matching only
    v = (v or "").strip()
    return v if v and v.upper() not in ("NA", "NAN", "NONE", "-") else None


def tips(path=TREE):
    return set(re.findall(r"[(,]\s*([^(),:;\s]+)", open(path).read()))


def physchem():
    """30_physchem.py as a module; one() reads the module global R at call time."""
    p = f"{R}/scripts/30_physchem.py"
    spec = importlib.util.spec_from_file_location("physchem30", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.R = R  # = production root unless DSC_ROOT is set (mock)
    m.SRC_MD5 = md5(p)
    return m


def genome_of(acc):  # = 10_miniprot_annot.py
    d = f"{R}/data/ncbi/{acc}"
    if os.path.exists(f"{d}/genomic.fna.gz"):
        return f"{d}/genomic.fna.gz"
    g = [x for x in glob.glob(f"{d}/ncbi_dataset/data/{acc}/*_genomic.fna.gz") if "cds_from" not in x]
    return g[0] if g else None


def listed_md5(acc, g):  # = 10_miniprot_annot.py
    d = f"{R}/data/ncbi/{acc}"
    if not os.path.exists(f"{d}/md5sum.txt"):
        return None
    rel = os.path.relpath(g, d)
    rel = rel[:-3] if rel.endswith(".gz") else rel
    for l in open(f"{d}/md5sum.txt"):
        x = l.split()
        if len(x) == 2 and os.path.normpath(x[1]) == rel:
            return x[0]
    return None


def level(t, r):  # = the level rule of 10_miniprot_annot.py pick_refs (raw strings compared)
    return 0 if r["family"] == t["family"] else 1 if r["order"] == t["order"] else 2 if r["class"] == t["class"] else 3


def pick_refs(t, rows, k=3):  # = 10_miniprot_annot.py pick_refs, returning (code, level)
    cands = []
    for r in rows:
        if r["code"] == t["code"] or r["source"] == "miniprot":
            continue
        if os.path.exists(f"{R}/work/annotfail/{r['code']}"):
            continue
        if not os.path.exists(f"{R}/data/proteomes/{r['code']}.faa"):
            continue
        if r["phylum"] != t["phylum"]:
            continue
        lvl = level(t, r)
        n = sum(1 for l in open(f"{R}/data/proteomes/{r['code']}.faa") if l.startswith(">"))
        refseq = r["accession"].startswith("GCF_")
        cands.append((lvl, not refseq, -n, r["code"]))
    cands.sort()
    return [(c[3], c[0]) for c in cands[:k]]


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def parse_aln(aln, plen):  # = 10_miniprot_annot.py: mRNA records with their CDS / stop_codon children
    models, cur = [], None
    with gzip.open(aln, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9:
                continue
            if x[2] == "mRNA":
                a = dict(kv.split("=", 1) for kv in x[8].split(";") if "=" in kv)
                tgt = a.get("Target", "").split()
                if not tgt:
                    cur = None
                    continue
                cov = (int(tgt[2]) - int(tgt[1]) + 1) / max(1, plen.get(tgt[0], 1))
                bad = ("StopCodon" in a) or ("Frameshift" in a)
                cur = dict(seq=x[0], s=int(x[3]), e=int(x[4]), strand=x[6],
                           score=float(x[5]) if x[5] != "." else 0.0, cov=cov, bad=bad,
                           ref=tgt[0], id=a["ID"], lines=[x])
                models.append(cur)
            elif cur is not None and x[2] in ("CDS", "stop_codon"):
                cur["lines"].append(x)
    return models


def choose(models):  # = 10_miniprot_annot.py: filters, then greedy by score, <= 20 % overlap on a strand
    keep = [m for m in models if m["cov"] >= MINCOV and not m["bad"]]
    keep.sort(key=lambda m: -m["score"])
    acc_by, chosen = {}, []
    for m in keep:
        key = (m["seq"], m["strand"])
        ok = True
        for o in acc_by.get(key, []):
            ov = overlap((m["s"], m["e"]), (o["s"], o["e"]))
            if ov > 0.2 * min(m["e"] - m["s"], o["e"] - o["s"]):
                ok = False
                break
        if ok:
            acc_by.setdefault(key, []).append(m)
            chosen.append(m)
    return keep, chosen


def write_models(code, chosen, mg):  # = 10_miniprot_annot.py models.gff
    with open(mg, "w") as o:
        o.write("##gff-version 3\n")
        for i, m in enumerate(sorted(chosen, key=lambda m: (m["seq"], m["s"])), 1):
            mid = f"{code}_{i:06d}"
            for x in m["lines"]:
                y = list(x)
                if y[2] == "mRNA":
                    y[8] = f"ID={mid};ref={m['ref']};cov={m['cov']:.2f}"
                else:
                    if y[2] == "stop_codon":
                        y[2] = "CDS"
                    y[8] = f"Parent={mid}"
                o.write("\t".join(y) + "\n")


def final_filter(peps, cdss, refof):  # = 10_miniprot_annot.py final filter (no internal stop, length, frame)
    out = []
    for h, p in peps.items():
        hid = h.split()[0]
        p = p.rstrip(".*")
        c = cdss.get(h, cdss.get(hid, "")).upper()
        if c[-3:] in ("TAA", "TAG", "TGA"):
            c = c[:-3]
        if "*" in p or "." in p or len(p) < MINLEN or len(c) % 3 or len(c) // 3 != len(p):
            continue
        out.append((hid, p, c, refof.get(hid)))
    return out


def compleasm_summary(p):  # summary.txt lines "S:81.25%, 546" ... "N:672"
    s = {}
    for l in open(p):
        m = re.match(r"\s*([A-Z]):\s*([0-9.]+)(%)?(?:,\s*([0-9]+))?", l)
        if not m:
            continue
        if m.group(3):
            s[f"cq_{m.group(1)}_pct"] = float(m.group(2))
            if m.group(4):
                s[f"cq_{m.group(1)}_n"] = int(m.group(4))
        else:
            s[f"cq_{m.group(1)}"] = int(float(m.group(2)))
    return s


def metrics_of(PC, code, prot, cds, ft, root):
    """30_physchem one() on the given files through a private view root/{data/proteomes,data/cds,work/protqc}."""
    for src, dst in ((prot, f"{root}/data/proteomes/{code}.faa"), (cds, f"{root}/data/cds/{code}.fna"),
                     (ft, f"{root}/work/protqc/{code}/full_table.tsv")):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        try:
            os.symlink(os.path.abspath(src), dst)
        except OSError:
            shutil.copy(src, dst)
    r0 = PC.R
    PC.R = root
    try:
        _, rows, sp = PC.one(code)
    finally:
        PC.R = r0
        shutil.rmtree(root, ignore_errors=True)
    return rows, sp


def write_rows(path, rows):
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")


def fail(wd, st, why, detail=""):
    st.update(why=why, detail=str(detail)[:2000], ended=time.strftime("%F %T"))
    with open(f"{wd}/REANNOT.FAIL", "w") as o:
        json.dump(st, o, indent=1)
    print(f"FAIL\t{st['code']}\t{why}\t{detail}", flush=True)
    return 4 if why == "no_model" else 5 if why in ("compleasm", "no_busco") else 3


def run_compleasm(prot, o, th):  # = 14_protqc_loop.sh pq_one: private lineage view, up to 3 tries
    L = f"{R}/db/compleasm"
    shutil.rmtree(o, ignore_errors=True)
    os.makedirs(f"{o}/L")
    for x in ("metazoa_odb12", "eukaryota_odb12", "placement_files"):
        for y in (x, x + ".done"):
            if os.path.exists(f"{L}/{y}"):
                os.symlink(f"{L}/{y}", f"{o}/L/{y}")
    for y in ("file_versions.tsv", "file_versions.tsv.done"):
        shutil.copy(f"{L}/{y}", f"{o}/L/{y}")
    env = dict(os.environ, PATH=f"{E}:{os.environ.get('PATH', '')}")
    ok = 0
    for k in (1, 2, 3):
        time.sleep(random.randint(0, 9) + (k - 1) * 30)
        with open(f"{o}/run.log", "w") as lg:
            rc = subprocess.run(["compleasm", "protein", "-p", prot, "-l", "metazoa_odb12", "-L", f"{o}/L",
                                 "-t", str(th), "-o", o], stdout=lg, stderr=subprocess.STDOUT, env=env).returncode
        if rc == 0:
            ok = k
            break
        os.replace(f"{o}/run.log", f"{o}/run.try{k}.log")
        for p in os.listdir(o):
            if p != "L" and not p.startswith("run.try"):
                pp = f"{o}/{p}"
                shutil.rmtree(pp) if os.path.isdir(pp) and not os.path.islink(pp) else os.remove(pp)
    shutil.rmtree(f"{o}/metazoa_odb12_hmmsearch_output", ignore_errors=True)
    shutil.rmtree(f"{o}/L", ignore_errors=True)
    return ok


def run_one(code, th, out, force=False):
    wd = f"{out}/{code}"
    if os.path.exists(f"{wd}/REANNOT.DONE") and not force:
        print(f"{code}\talready done", flush=True)
        return 0
    os.makedirs(wd, exist_ok=True)
    for p in os.listdir(wd):  # fresh start; a complete alignment (written atomically) is re-used if refs agree
        if p not in ("aln.gff.gz", "aln.refs"):
            pp = f"{wd}/{p}"
            shutil.rmtree(pp) if os.path.isdir(pp) and not os.path.islink(pp) else os.remove(pp)
    T0 = time.time()
    rows = tsv(f"{R}/meta/species_codes.tsv")
    t = next(r for r in rows if r["code"] == code)
    acc = t["accession"]
    st = dict(code=code, accession=acc, phylum=t["phylum"], started=time.strftime("%F %T"), host=platform.node(),
              threads=th, script_md5=md5(os.path.abspath(__file__)))
    g = genome_of(acc)
    if not g:
        return fail(wd, st, "no_genome", acc)
    st.update(genome=g, genome_gz_bytes=os.path.getsize(g))
    refs = pick_refs(t, rows)
    if not refs:
        return fail(wd, st, "no_ref", t["phylum"])
    st.update(refs=[c for c, _ in refs], ref_levels=[LEVELS[l] for _, l in refs])
    fa = f"{wd}/genome.fa"
    t1 = time.time()
    try:
        sh(f"zcat {q(g)} > {q(fa)}.tmp")
        exp, got = listed_md5(acc, g), md5(f"{fa}.tmp")
        why = "md5_mismatch" if exp and got != exp else None
    except subprocess.CalledProcessError:
        why, exp, got = "gzip_error", None, None
    st.update(md5_listed=exp, md5_genome=got)
    if why:
        if os.path.exists(f"{fa}.tmp"):
            os.remove(f"{fa}.tmp")
        return fail(wd, st, why, g)
    os.replace(f"{fa}.tmp", fa)
    sh(f"{E}/samtools faidx {q(fa)} || {E}/seqkit faidx {q(fa)}")
    asm = sum(int(l.split("\t")[1]) for l in open(fa + ".fai"))
    st.update(asm_bp=asm, t_genome_s=round(time.time() - t1, 1))
    if asm < MIN_ASM:
        for p in (fa, fa + ".fai"):
            os.remove(p)
        return fail(wd, st, "organelle_only", asm)
    with open(f"{wd}/refs.txt", "w") as o:
        o.write("".join(f"{c}\t{LEVELS[l]}\n" for c, l in refs))
    ref, plen, multi = f"{wd}/ref.faa", {}, 0
    with open(ref, "w") as o:
        for c, _ in refs:
            name = None
            for line in open(f"{R}/data/proteomes/{c}.faa"):
                o.write(line)
                if line.startswith(">"):
                    name = line[1:].strip()
                    plen[name] = 0
                    multi += len(name.split()) > 1
                else:
                    plen[name] += len(line.strip())
    st.update(ref_proteins=len(plen), ref_multitoken_headers=multi)
    aln, key = f"{wd}/aln.gff.gz", ",".join(st["refs"])
    if os.path.exists(aln) and not (os.path.exists(f"{wd}/aln.refs") and open(f"{wd}/aln.refs").read().strip() == key):
        os.remove(aln)
    st["aln_reused"] = os.path.exists(aln)
    t1 = time.time()
    if not st["aln_reused"]:
        sh(f"{E}/miniprot -t {th} --gff -I --outn=3 {q(fa)} {q(ref)} | gzip > {q(aln)}.tmp && mv {q(aln)}.tmp {q(aln)}")
        with open(f"{wd}/aln.refs", "w") as o:
            o.write(key + "\n")
    st["t_miniprot_s"] = round(time.time() - t1, 1)
    models = parse_aln(aln, plen)
    keep, chosen = choose(models)
    mg = f"{wd}/models.gff"
    write_models(code, chosen, mg)
    t1 = time.time()
    sh(f"{E}/gffread {q(mg)} -g {q(fa)} -x {q(wd)}/cds.fa -y {q(wd)}/pep.fa -S")
    PC = physchem()
    peps, cdss = PC.read_fa(f"{wd}/pep.fa"), PC.read_fa(f"{wd}/cds.fa")
    st["pep_multitoken_headers"] = sum(1 for l in open(f"{wd}/pep.fa") if l.startswith(">") and len(l[1:].split()) > 1)
    refof = {}
    for line in open(mg):
        if "\tmRNA\t" in line:
            a = dict(kv.split("=", 1) for kv in line.rstrip().split("\t")[8].split(";"))
            refof[a["ID"]] = a["ref"]
    fin = final_filter(peps, cdss, refof)
    with open(f"{wd}/prot.faa", "w") as fp, open(f"{wd}/cds.fna", "w") as fc, open(f"{wd}/idmap.tsv", "w") as mp:
        for hid, p, c, rf in fin:
            fp.write(f">{hid}\n{p}\n")
            fc.write(f">{hid}\n{c}\n")
            mp.write(f"{hid}\tminiprot\t{rf}\n")
    for p in (fa, fa + ".fai", ref, f"{wd}/pep.fa", f"{wd}/cds.fa"):
        if os.path.exists(p):
            os.remove(p)
    sh(f"gzip -f {q(mg)}")
    st.update({"aln_models": len(models), "pass": len(keep), "chosen": len(chosen), "final": len(fin),
               "t_gffread_s": round(time.time() - t1, 1)})
    print(f"{code}\trefs={key}\taln_models={len(models)}\tpass={len(keep)}\tchosen={len(chosen)}\tfinal={len(fin)}",
          flush=True)
    if not fin:
        return fail(wd, st, "no_model", "no gene model passed the filters")
    t1 = time.time()
    ok = run_compleasm(f"{wd}/prot.faa", f"{wd}/protqc", th)
    st.update(compleasm_try=ok, t_compleasm_s=round(time.time() - t1, 1))
    ft = f"{wd}/protqc/full_table.tsv"
    if not ok or not (os.path.exists(ft) and os.path.exists(f"{wd}/protqc/summary.txt")):
        return fail(wd, st, "compleasm", f"see {wd}/protqc/run.try*.log")
    cq = compleasm_summary(f"{wd}/protqc/summary.txt")
    t1 = time.time()
    rows_g, sp = metrics_of(PC, code, f"{wd}/prot.faa", f"{wd}/cds.fna", ft, f"{wd}/_pcroot")
    if sp is None or not rows_g:
        return fail(wd, st, "no_busco", "no compleasm Single hit with protein metrics")
    write_rows(f"{wd}/metrics.tsv", [dict(sp, **cq)])
    write_rows(f"{wd}/gene_metrics.tsv", rows_g)
    _, _, sp0 = PC.one(code)  # production proteome re-measured (read only): reproduces work/physchem?
    if sp0:
        write_rows(f"{wd}/prod_metrics.tsv", [sp0])
    st.update(physchem_md5=PC.SRC_MD5, n_busco_genes=sp["n_busco_genes"], n_proteins=sp["n_proteins"],
              t_metrics_s=round(time.time() - t1, 1), t_total_s=round(time.time() - T0, 1),
              ended=time.strftime("%F %T"), **cq)
    with open(f"{wd}/REANNOT.DONE", "w") as o:
        json.dump(st, o, indent=1)
    print(f"{code}\tdone\tS+D={cq.get('cq_S_pct', 0) + cq.get('cq_D_pct', 0):.2f}\tbusco_genes={sp['n_busco_genes']}"
          f"\t{st['t_total_s']} s", flush=True)
    return 0


def select(out, n, seed, gzmax, force=False):
    if os.path.exists(f"{out}/SELECT.DONE") and not force:
        print(f"{out}/SELECT.DONE exists: selection kept (--force to redo)", flush=True)
        return 0
    os.makedirs(out, exist_ok=True)
    inputs = (TREE, f"{R}/meta/labels_v2.tsv", f"{R}/work/physchem/species_metrics.tsv", f"{R}/meta/species_codes.tsv")
    tip = tips()
    lab = {r["code"]: r["lab"] for r in tsv(inputs[1])}
    sm = {r["code"]: r for r in tsv(inputs[2])}
    meta = {r["code"]: r for r in tsv(inputs[3])}
    cand = []
    for c in sorted(tip):
        src = (sm.get(c) or meta.get(c) or {}).get("source", "")
        if lab.get(c) not in ("D", "S") or src not in ANNOT or c not in meta:
            continue
        m = meta[c]
        g = genome_of(m["accession"])
        gz = os.path.getsize(g) if g else ""
        why = ("no_physchem" if c not in sm else
               "no_proteome" if not os.path.exists(f"{R}/data/proteomes/{c}.faa") else
               "no_full_table" if not os.path.exists(f"{R}/work/protqc/{c}/full_table.tsv") else
               "annotfail" if os.path.exists(f"{R}/work/annotfail/{c}") else
               "no_genome" if not g else "genome_gz_too_large" if gz > gzmax else "")
        cand.append({"code": c, "species": m.get("species", ""), "lab": lab[c], "source": src, "phylum": m["phylum"],
                     "class": m["class"], "order": m["order"], "family": m["family"], "accession": m["accession"],
                     "genome_gz_bytes": gz, "eligible": int(not why), "reason": why or "eligible"})
    write_rows(f"{out}/candidates.tsv", cand)
    by = {}
    for r in cand:
        if r["eligible"]:
            by.setdefault(r["phylum"], {"D": [], "S": []})[r["lab"]].append(r)
    m = {p: min(len(v["D"]), len(v["S"])) for p, v in by.items()}
    cap = 0  # largest c with sum_p min(m_p, c) <= n / 2
    while cap < max(m.values(), default=0) and sum(min(x, cap + 1) for x in m.values()) <= n // 2:
        cap += 1
    quota = {p: min(x, cap) for p, x in m.items()}
    sel = []
    for p in sorted(quota):
        if not quota[p]:
            continue
        rng = random.Random(f"{seed}|{p}")
        D = sorted(by[p]["D"], key=lambda r: r["code"])
        S = sorted(by[p]["S"], key=lambda r: r["code"])
        used = set()
        for k, d in enumerate(rng.sample(D, quota[p]), 1):
            pool = [s for s in S if s["code"] not in used]
            best, lvl = pool, "phylum"
            for key in ("order", "class"):
                v = tx(d[key])
                b = [s for s in pool if v and tx(s[key]) == v]
                if b:
                    best, lvl = b, key
                    break
            s = rng.choice(best)
            used.add(s["code"])
            for r in (d, s):
                sel.append(dict({"pair": f"{p}_{k:02d}", "draw": k, "match": lvl}, **r))
    write_rows(f"{out}/selection.tsv", sel)
    rep = [f"59_route_reannot.py --select {time.strftime('%F %T')}  python {platform.python_version()}",
           f"script md5 {md5(os.path.abspath(__file__))}; seed {seed}; n {n}; gzmax {gzmax:g}"]
    rep += [f"input {f} md5 {md5(f)}" for f in inputs]
    rep += [f"population (labelled annot-route tree species) {len(cand)}; eligible {sum(r['eligible'] for r in cand)}",
            "ineligible: " + ", ".join(f"{k} {v}" for k, v in sorted(Counter(r["reason"] for r in cand
                                                                            if not r["eligible"]).items())),
            f"cap c = {cap}; pairs {len(sel) // 2}; species {len(sel)}", "phylum\tnD\tnS\tquota"]
    rep += [f"{p}\t{len(by[p]['D'])}\t{len(by[p]['S'])}\t{quota[p]}" for p in sorted(by)]
    rep += ["pair\tD\tS\tmatch"] + [f"{a['pair']}\t{a['code']}\t{b['code']}\t{a['match']}" for a, b in zip(sel[::2], sel[1::2])]
    with open(f"{out}/select_report.txt", "w") as o:
        o.write("\n".join(rep) + "\n")
    with open(f"{out}/SELECT.DONE", "w") as o:
        o.write(time.strftime("%F %T") + f"\t{len(sel)} species\n")
    print("\n".join(rep), flush=True)
    return 0


def _fin(v):
    return [x for x in v if math.isfinite(x)]


def collect(out):
    from scipy import stats
    nan = float("nan")

    def tci(v):
        if len(v) < 2:
            return (v[0] if v else nan), nan, nan
        h = float(stats.t.ppf(0.975, len(v) - 1)) * stdev(v) / math.sqrt(len(v))
        return mean(v), mean(v) - h, mean(v) + h

    def safe(f, *a):  # descriptive tests; degenerate input (n < 2, all zero, constant) -> NaN
        try:  # tuple unpacking: old scipy namedtuples (SpearmanrResult.correlation) and new result objects
            s, p = f(*a)
            return float(s), float(p)
        except Exception:
            return nan, nan

    sel = tsv(f"{out}/selection.tsv")
    sm = {r["code"]: r for r in tsv(f"{R}/work/physchem/species_metrics.tsv")}
    meta = {r["code"]: r for r in tsv(f"{R}/meta/species_codes.tsv")}
    qc = {}
    for r in tsv(f"{R}/meta/qc_table.tsv"):
        qc.setdefault(r["code"], r)
    lab = {r["code"]: r["lab"] for r in tsv(f"{R}/meta/labels_v2.tsv")}
    tip = tips()
    traits = [k for k in next(iter(sm.values())) if k.startswith(("busco_", "whole_"))] + ["n_busco_genes", "n_proteins"]
    pop = [c for c in sm if c in tip and lab.get(c) in ("D", "S")]
    sd = {}
    for k in traits:
        v = _fin([fnum(sm[c].get(k)) for c in pop])
        sd[k] = stdev(v) if len(v) > 2 else nan
    cs = f"{R}/work/final_v3/pgls_controls/controls_summary.tsv"
    deep = {r["trait"]: fnum(r.get("all_est")) for r in tsv(cs)} if os.path.exists(cs) else {}
    sp_rows, mp, lvl = [], {}, {}
    for r in sel:
        c, wd = r["code"], f"{out}/{r['code']}"
        st, prod, status = {}, {}, "missing"
        if os.path.exists(f"{wd}/REANNOT.DONE"):
            st, status = json.load(open(f"{wd}/REANNOT.DONE")), "done"
            mp[c] = tsv(f"{wd}/metrics.tsv")[0]
            if os.path.exists(f"{wd}/prod_metrics.tsv"):
                prod = tsv(f"{wd}/prod_metrics.tsv")[0]
            lvl[c] = LEVELS.index(st["ref_levels"][0])
        elif os.path.exists(f"{wd}/REANNOT.FAIL"):
            st, status = json.load(open(f"{wd}/REANNOT.FAIL")), "fail"
        rep = _fin([abs(fnum(prod[k]) - fnum(sm[c][k])) for k in traits if k in prod and k in sm.get(c, {})])
        m, qa = mp.get(c, {}), qc.get(c, {})
        sp_rows.append({"pair": r["pair"], "lab": r["lab"], "code": c, "phylum": r["phylum"], "status": status,
                        "why": st.get("why", ""), "refs": ",".join(st.get("refs", [])),
                        "ref_levels": ",".join(st.get("ref_levels", [])), "asm_bp": st.get("asm_bp", ""),
                        "md5_listed": st.get("md5_listed", ""), "md5_genome": st.get("md5_genome", ""),
                        **{k: st.get(k, "") for k in ("aln_models", "pass", "chosen", "final")},
                        "n_prot_annot": sm.get(c, {}).get("n_proteins", ""), "n_prot_mp": m.get("n_proteins", ""),
                        "n_busco_annot": sm.get(c, {}).get("n_busco_genes", ""), "n_busco_mp": m.get("n_busco_genes", ""),
                        "protC_annot": qa.get("prot_C", ""),
                        "protC_mp": fnum(m.get("cq_S_pct")) + fnum(m.get("cq_D_pct")) if m else "",
                        "repro_maxabs": max(rep) if rep else "",
                        **{k: st.get(k, "") for k in ("t_genome_s", "t_miniprot_s", "t_compleasm_s", "t_total_s")}})
    write_rows(f"{out}/species.tsv", sp_rows)
    D = []
    for r in sel:
        c = r["code"]
        if c not in mp:
            continue
        for k in traits:
            a, b = fnum(sm[c].get(k)), fnum(mp[c].get(k))
            D.append({"code": c, "lab": r["lab"], "pair": r["pair"], "phylum": r["phylum"], "trait": k, "annot": a,
                      "miniprot": b, "delta_raw": a - b, "sd": sd[k],
                      "delta_sd": (a - b) / sd[k] if sd[k] > 0 else nan})
    write_rows(f"{out}/deltas.tsv", D or [{"code": "", "trait": ""}])
    S = []
    for k in traits:
        rows = [d for d in D if d["trait"] == k and math.isfinite(d["delta_sd"])]
        v = [d["delta_sd"] for d in rows]
        vD, vS = [d["delta_sd"] for d in rows if d["lab"] == "D"], [d["delta_sd"] for d in rows if d["lab"] == "S"]
        pr = {}
        for d in rows:
            pr.setdefault(d["pair"], {})[d["lab"]] = d["delta_sd"]
        pv = [x["D"] - x["S"] for x in pr.values() if "D" in x and "S" in x]
        mu, lo, hi = tci(v)
        dDS = (mean(vD) - mean(vS)) if vD and vS else nan
        de = deep.get(k, nan)
        rat = (lambda x: x / de if math.isfinite(de) and de != 0 else nan)
        a, b = [d["annot"] for d in rows], [d["miniprot"] for d in rows]
        S.append({"trait": k, "class": "P6" if k in P6 else "gc" if k in GCFAM else "other", "n": len(v),
                  "n_D": len(vD), "n_S": len(vS), "mean_sd": mu, "ci_lo": lo, "ci_hi": hi,
                  "wilcoxon_p": safe(stats.wilcoxon, v)[1], "mean_D": mean(vD) if vD else nan,
                  "mean_S": mean(vS) if vS else nan, "diff_DS": dDS,
                  "mw_p": safe(lambda x, y: stats.mannwhitneyu(x, y, alternative="two-sided"), vD, vS)[1],
                  "n_pairs": len(pv), "pair_mean": mean(pv) if pv else nan, "pair_wilcoxon_p": safe(stats.wilcoxon, pv)[1],
                  "pearson_r": safe(stats.pearsonr, a, b)[0], "spearman_rho": safe(stats.spearmanr, a, b)[0],
                  "deep_est": de, "ratio_mean": rat(mu), "ratio_DS": rat(dDS), "ratio_pair": rat(mean(pv) if pv else nan)})
    write_rows(f"{out}/summary.tsv", S)
    RL = []  # first (best) reference: re-annotated species vs the production miniprot tree species
    for r in sel:
        c, wd = r["code"], f"{out}/{r['code']}"
        if os.path.exists(f"{wd}/refs.txt"):
            x = open(f"{wd}/refs.txt").readline().split()
            RL.append({"group": "reannot", "code": c, "lab": r["lab"], "phylum": r["phylum"],
                       "ref1": x[0] if x else "NA", "level": x[1] if len(x) > 1 else "NA",
                       "status": "done" if c in mp else "fail" if os.path.exists(f"{wd}/REANNOT.FAIL") else "missing"})
    for c in sorted(tip):
        if lab.get(c) not in ("D", "S") or (sm.get(c) or meta.get(c) or {}).get("source", "") != "miniprot":
            continue
        f = f"{R}/work/miniprot/{c}/refs.txt"
        x = open(f).readline().split() if os.path.exists(f) else []
        lv = LEVELS[level(meta[c], meta[x[0]])] if x and c in meta and x[0] in meta else "NA"
        RL.append({"group": "production_miniprot", "code": c, "lab": lab[c], "phylum": meta.get(c, {}).get("phylum", "NA"),
                   "ref1": x[0] if x else "NA", "level": lv, "status": "production"})
    write_rows(f"{out}/reflevel.tsv", RL or [{"group": ""}])

    def f3(v):
        return f"{v:.3g}" if isinstance(v, (int, float)) and math.isfinite(v) else "NA"

    cnt = Counter(x["status"] for x in sp_rows)
    rx = _fin([fnum(x["repro_maxabs"]) for x in sp_rows])
    rep = [f"59_route_reannot.py --collect {time.strftime('%F %T')}  python {platform.python_version()}",
           f"script md5 {md5(os.path.abspath(__file__))}; out {out}",
           f"selected {len(sel)}: done {cnt['done']}, fail {cnt['fail']}, missing {cnt['missing']}"]
    rep += [f"  FAIL {x['code']} ({x['lab']}): {x['why']}" for x in sp_rows if x["status"] == "fail"]
    rep += [f"production proteome re-measured by 30_physchem one(): max |diff| vs species_metrics {f3(max(rx))}"
            if rx else "production proteome re-measure: none",
            "delta = annot - miniprot in SD units of the labelled tree species (positive = miniprot route lower);"
            " diff_DS = mean delta D - mean delta S; ratio_DS = diff_DS / 44 ref_all deep estimate",
            "trait\tn(D/S)\tmean_sd [95% CI]\twilcoxon_p\tdiff_DS\tmw_p\tpairs\tpair_mean\tpair_p\tdeep_est\tratio_DS"]
    Sk = {s["trait"]: s for s in S}
    for name, ks in (("P6", P6), ("GC set", GCFAM)):
        rep.append(f"-- {name}")
        for k in ks:
            if k in Sk:
                s = Sk[k]
                rep.append(f"{k}\t{s['n']}({s['n_D']}/{s['n_S']})\t{f3(s['mean_sd'])} [{f3(s['ci_lo'])}, {f3(s['ci_hi'])}]"
                           f"\t{f3(s['wilcoxon_p'])}\t{f3(s['diff_DS'])}\t{f3(s['mw_p'])}\t{s['n_pairs']}"
                           f"\t{f3(s['pair_mean'])}\t{f3(s['pair_wilcoxon_p'])}\t{f3(s['deep_est'])}\t{f3(s['ratio_DS'])}")
    rep.append("first-reference level by group and label (" + " ".join(LEVELS) + " NA)")
    for gname in ("reannot", "production_miniprot"):
        for lb in ("D", "S"):
            cc = Counter(x["level"] for x in RL if x["group"] == gname and x["lab"] == lb)
            rep.append(f"{gname}\t{lb}\t" + "\t".join(f"{l} {cc.get(l, 0)}" for l in LEVELS + ("NA",)))
    rep.append("Spearman |delta_sd| vs first-reference level (0 family .. 3 phylum), done species")
    for k in P6:
        z = [(lvl[d["code"]], abs(d["delta_sd"])) for d in D
             if d["trait"] == k and d["code"] in lvl and math.isfinite(d["delta_sd"])]
        rho, p = safe(stats.spearmanr, [a for a, _ in z], [b for _, b in z]) if len(z) > 2 else (nan, nan)
        rep.append(f"{k}\tn {len(z)}\trho {f3(rho)}\tp {f3(p)}")
    with open(f"{out}/report.txt", "w") as o:
        o.write("\n".join(rep) + "\n")
    with open(f"{out}/COLLECT.DONE", "w") as o:
        o.write(time.strftime("%F %T") + f"\tdone {cnt['done']} fail {cnt['fail']} missing {cnt['missing']}\n")
    print("\n".join(rep), flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser(description="A1 re-annotation control (PREREG_addon.md, A1)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--select", action="store_true")
    g.add_argument("--one", metavar="CODE")
    g.add_argument("--collect", action="store_true")
    ap.add_argument("--out", default=f"{R}/work/addon/A1_reannot")
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20261025)
    ap.add_argument("--gzmax", type=float, default=1e9)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    if a.select:
        return select(out, a.n, a.seed, a.gzmax, a.force)
    if a.one:
        return run_one(a.one, a.threads, out, a.force)
    return collect(out)


if __name__ == "__main__":
    sys.exit(main())
