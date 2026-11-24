#!/usr/bin/env python3
"""Homology-projection annotation of an unannotated assembly with miniprot.

Reference = proteomes of the up-to-3 taxonomically closest annotated species
(same family > order > class > phylum; RefSeq preferred, then gene count),
concatenated. miniprot spliced alignments -> greedy non-overlapping best-score
gene models (coverage >= 0.5 of reference protein, no in-frame stop or
frameshift) -> CDS/protein via gffread.

Input checks (2026-10-19): the decompressed genome is verified against the NCBI
md5sum.txt when listed (mismatch -> .done removed and accession queued in
work/integrity/redownload.txt, exit 0 so the orthology gates keep waiting);
assemblies < 1 Mb are organelle-only records and fail; a run that yields no gene
model fails instead of leaving an empty proteome (the annotation loop re-ran
such species forever). Shell steps run with pipefail.

Usage: python 10_miniprot_annot.py CODE [threads]
Writes data/proteomes/CODE.faa, data/cds/CODE.fna, meta/idmap/CODE.tsv,
work/miniprot/CODE/{refs.txt,aln.gff.gz,models.gff}
"""
import csv, glob, gzip, hashlib, os, re, subprocess, sys

R = os.path.expanduser("~/deepsea_conv")
E = f"{R}/envs/dc/bin"
MINCOV, MINLEN = 0.5, 30
MIN_ASM = 1_000_000  # bp; smaller records are mitogenome/organelle-only deposits


def sh(cmd):  # pipefail: a crashed producer must not leave a truncated cache behind
    subprocess.run(["/bin/bash", "-o", "pipefail", "-c", cmd], check=True)


def table():
    return list(csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"))


def genome_of(acc):
    d = f"{R}/data/ncbi/{acc}"
    if os.path.exists(f"{d}/genomic.fna.gz"):
        return f"{d}/genomic.fna.gz"
    g = [x for x in glob.glob(f"{d}/ncbi_dataset/data/{acc}/*_genomic.fna.gz") if "cds_from" not in x]
    return g[0] if g else None


def listed_md5(acc, g):
    """md5 of the uncompressed genome as listed in the datasets md5sum.txt, else None."""
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


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def pick_refs(t, rows, k=3):
    cands = []
    for r in rows:
        if r["code"] == t["code"] or r["source"] == "miniprot":
            continue
        if os.path.exists(f"{R}/work/annotfail/{r['code']}"):  # projected, not a real annotation
            continue
        if not os.path.exists(f"{R}/data/proteomes/{r['code']}.faa"):
            continue
        if r["phylum"] != t["phylum"]:
            continue
        lvl = 0 if r["family"] == t["family"] else 1 if r["order"] == t["order"] else \
            2 if r["class"] == t["class"] else 3
        n = sum(1 for l in open(f"{R}/data/proteomes/{r['code']}.faa") if l.startswith(">"))
        refseq = r["accession"].startswith("GCF_")
        cands.append((lvl, not refseq, -n, r["code"]))
    cands.sort()
    return [c[3] for c in cands[:k]]


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def main():
    code = sys.argv[1]; th = sys.argv[2] if len(sys.argv) > 2 else "16"
    rows = table(); t = next(r for r in rows if r["code"] == code)
    acc = t["accession"]
    wd = f"{R}/work/miniprot/{code}"; os.makedirs(wd, exist_ok=True)
    g = genome_of(acc)
    if not g:
        sys.exit(f"no genome for {code}")
    refs = pick_refs(t, rows)
    if not refs:
        sys.exit(f"no annotated reference in phylum {t['phylum']} for {code}")
    # decompress genome once (gffread needs random access), verified against NCBI md5
    fa = f"{wd}/genome.fa"
    if not (os.path.exists(fa + ".fai") and os.path.exists(fa + ".ok")):
        try:
            sh(f"zcat {g} > {fa}.tmp")
            exp = listed_md5(acc, g)
            why = "md5_mismatch" if exp and md5(f"{fa}.tmp") != exp else None
        except subprocess.CalledProcessError:
            why = "gzip_error"
        if why:
            if os.path.exists(f"{fa}.tmp"):
                os.remove(f"{fa}.tmp")
            os.makedirs(f"{R}/work/integrity", exist_ok=True)
            with open(f"{R}/work/integrity/redownload.txt", "a") as o:
                o.write(f"{code}\t{acc}\t{why}\t{g}\n")
            done = f"{R}/data/ncbi/{acc}/.done"
            if os.path.exists(done):
                os.remove(done)  # leaves the annotation queue; no .fail, so orthology gates keep waiting
            print(f"{code}\t{why} for {g}: queued for re-download", flush=True)
            return
        os.replace(f"{fa}.tmp", fa)
        sh(f"{E}/samtools faidx {fa} || {E}/seqkit faidx {fa}")
        open(fa + ".ok", "w").close()
    asm = sum(int(l.split("\t")[1]) for l in open(fa + ".fai"))
    if asm < MIN_ASM:
        for p in (fa, fa + ".fai", fa + ".ok"):
            os.remove(p)
        sys.exit(f"{code}: assembly is {asm} bp (< {MIN_ASM}): organelle-only record, not a nuclear genome")
    open(f"{wd}/refs.txt", "w").write("\n".join(refs) + "\n")
    ref = f"{wd}/ref.faa"
    plen = {}
    with open(ref, "w") as o:
        for c in refs:
            name = None
            for line in open(f"{R}/data/proteomes/{c}.faa"):
                o.write(line)
                if line.startswith(">"):
                    name = line[1:].strip(); plen[name] = 0
                else:
                    plen[name] += len(line.strip())
    aln = f"{wd}/aln.gff.gz"
    if not os.path.exists(aln):
        sh(f"{E}/miniprot -t {th} --gff -I --outn=3 {fa} {ref} | gzip > {aln}.tmp && mv {aln}.tmp {aln}")
    # parse mRNA + children
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
                    cur = None; continue
                cov = (int(tgt[2]) - int(tgt[1]) + 1) / max(1, plen.get(tgt[0], 1))
                bad = ("StopCodon" in a) or ("Frameshift" in a)
                cur = dict(seq=x[0], s=int(x[3]), e=int(x[4]), strand=x[6],
                           score=float(x[5]) if x[5] != "." else 0.0, cov=cov, bad=bad,
                           ref=tgt[0], id=a["ID"], lines=[x])
                models.append(cur)
            elif cur is not None and x[2] in ("CDS", "stop_codon"):
                cur["lines"].append(x)
    keep = [m for m in models if m["cov"] >= MINCOV and not m["bad"]]
    keep.sort(key=lambda m: -m["score"])
    acc_by = {}
    chosen = []
    for m in keep:
        key = (m["seq"], m["strand"])
        ok = True
        for o in acc_by.get(key, []):
            ov = overlap((m["s"], m["e"]), (o["s"], o["e"]))
            if ov > 0.2 * min(m["e"] - m["s"], o["e"] - o["s"]):
                ok = False; break
        if ok:
            acc_by.setdefault(key, []).append(m); chosen.append(m)
    mg = f"{wd}/models.gff"
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
    sh(f"{E}/gffread {mg} -g {fa} -x {wd}/cds.fa -y {wd}/pep.fa -S")
    # final filter: no internal stops, length, frame
    import importlib; loc = importlib.import_module("06_prep_local")
    peps = dict(loc.read_fa(f"{wd}/pep.fa")); cdss = dict(loc.read_fa(f"{wd}/cds.fa"))
    refof = {m_id: None for m_id in peps}
    for line in open(mg):
        if "\tmRNA\t" in line:
            a = dict(kv.split("=", 1) for kv in line.rstrip().split("\t")[8].split(";"))
            refof[a["ID"]] = a["ref"]
    n = 0
    out = (f"{R}/data/proteomes/{code}.faa", f"{R}/data/cds/{code}.fna", f"{R}/meta/idmap/{code}.tsv")
    with open(out[0], "w") as fa_o, open(out[1], "w") as fn_o, open(out[2], "w") as mp:
        for h, p in peps.items():
            hid = h.split()[0]
            p = p.rstrip(".*")
            c = cdss.get(h, cdss.get(hid, "")).upper()
            if c[-3:] in ("TAA", "TAG", "TGA"):
                c = c[:-3]
            if "*" in p or "." in p or len(p) < MINLEN or len(c) % 3 or len(c) // 3 != len(p):
                continue
            fa_o.write(f">{hid}\n{p}\n"); fn_o.write(f">{hid}\n{c}\n")
            mp.write(f"{hid}\tminiprot\t{refof.get(hid)}\n"); n += 1
    for p in (fa, fa + ".fai", fa + ".ok", ref):
        os.remove(p)
    print(f"{code}\trefs={','.join(refs)}\taln_models={len(models)}\tpass={len(keep)}\tchosen={len(chosen)}\tfinal={n}", flush=True)
    if n == 0:
        for p in out:
            os.remove(p)
        sys.exit(f"{code}: no gene model passed the filters")


if __name__ == "__main__":
    main()
