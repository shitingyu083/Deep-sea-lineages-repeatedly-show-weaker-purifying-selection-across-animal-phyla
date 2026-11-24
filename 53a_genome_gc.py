#!/usr/bin/env python3
"""Genome GC of every non-TSA species (add-on A1, PREREG_addon.md).

gc_fasta = (G + C) / (A + C + G + T) over the assembly FASTA on disk (genome_of = 10_miniprot_annot.py;
case-insensitive; N and other IUPAC codes excluded; header lines skipped, any line length); ncbi_gc =
assemblyStats.gcPercent / 100 of the datasets assembly_data_report.jsonl. genome_gc = gc_fasta, else ncbi_gc
(gc_src fasta / ncbi / none); flag_diff = 1 if |gc_fasta - ncbi_gc| > 0.01. Genomes of local species: optional
work/addon/local_genomes.tsv (code, path). Per-species cache work/addon/gc_cache/CODE.json (resumable; failed reads
are not cached). Read-only on data/; run under nice / ionice.
Usage: python 53a_genome_gc.py [workers=8]      Writes work/addon/genome_gc.tsv
"""
import csv, glob, json, os, shutil, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed

R = os.path.expanduser("~/deepsea_conv")
OUT = f"{R}/work/addon"
CACHE = f"{OUT}/gc_cache"
BLOCK = 1 << 24


def genome_of(acc):  # = 10_miniprot_annot.genome_of
    d = f"{R}/data/ncbi/{acc}"
    if os.path.exists(f"{d}/genomic.fna.gz"):
        return f"{d}/genomic.fna.gz"
    g = [x for x in glob.glob(f"{d}/ncbi_dataset/data/{acc}/*_genomic.fna.gz") if "cds_from" not in x]
    return g[0] if g else None


def ncbi_gc(acc):
    for rep in sorted(glob.glob(f"{R}/data/ncbi/{acc}/**/assembly_data_report.jsonl", recursive=True)):
        for l in open(rep):
            try:
                v = json.loads(l).get("assemblyStats", {}).get("gcPercent")
            except ValueError:
                continue
            if v is not None:
                return float(v) / 100
    return None


def count(path):
    """base counts of a (gzipped) FASTA, streaming; '>' only occurs in header lines"""
    if path.endswith(".gz"):
        cmd = ["pigz", "-dc", path] if shutil.which("pigz") else ["zcat", path]
    else:
        cmd = ["cat", path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=BLOCK)
    n = dict.fromkeys(["A", "C", "G", "T", "N", "letters"], 0)
    in_hdr = False
    for block in iter(lambda: p.stdout.read(BLOCK), b""):
        i, parts = 0, []
        if in_hdr:
            j = block.find(b"\n")
            if j < 0:
                continue
            i, in_hdr = j + 1, False
        while True:
            j = block.find(b">", i)
            if j < 0:
                parts.append(block[i:]); break
            parts.append(block[i:j])
            k = block.find(b"\n", j)
            if k < 0:
                in_hdr = True; break
            i = k + 1
        s = b"".join(parts)
        for b in "ACGTN":
            n[b] += s.count(b.encode()) + s.count(b.lower().encode())
        n["letters"] += len(s) - s.count(b"\n") - s.count(b"\r") - s.count(b" ")
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError(f"{cmd[0]} exit {p.returncode} on {path}")
    return n


def one(code, acc, src, path):
    cf = f"{CACHE}/{code}.json"
    if os.path.exists(cf):
        return json.load(open(cf))
    r = dict(code=code, accession=acc, source=src, file=path or "", ncbi_gc=ncbi_gc(acc) if acc else None,
             gc_fasta=None, n_acgt=0, n_N=0, n_other=0, error="")
    if path and os.path.exists(path):
        try:
            n = count(path)
            acgt = n["A"] + n["C"] + n["G"] + n["T"]
            r.update(n_acgt=acgt, n_N=n["N"], n_other=n["letters"] - acgt - n["N"],
                     gc_fasta=(n["G"] + n["C"]) / acgt if acgt else None)
        except Exception as e:  # not cached: retried on the next run
            r["error"] = str(e)[:200]
            return r
    json.dump(r, open(cf + ".tmp", "w")); os.replace(cf + ".tmp", cf)
    return r


def main():
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    os.makedirs(CACHE, exist_ok=True)
    loc = {}
    if os.path.exists(f"{OUT}/local_genomes.tsv"):
        loc = {r["code"]: r["path"] for r in csv.DictReader(open(f"{OUT}/local_genomes.tsv"), delimiter="\t")}
    jobs = []
    for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"):
        if r["source"] == "tsa":
            continue
        acc = r["accession"]
        jobs.append((r["code"], acc, r["source"], loc.get(r["code"]) or (genome_of(acc) if acc else None)))
    jobs.sort(key=lambda j: -(os.path.getsize(j[3]) if j[3] and os.path.exists(j[3]) else 0))  # largest first
    res = []
    with ProcessPoolExecutor(nw) as ex:
        fs = [ex.submit(one, *j) for j in jobs]
        for f in as_completed(fs):
            res.append(f.result())
    cols = ["code", "accession", "source", "file", "n_acgt", "n_N", "n_other", "gc_fasta", "ncbi_gc", "genome_gc",
            "gc_src", "flag_diff", "error"]
    with open(f"{OUT}/genome_gc.tsv.tmp", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in sorted(res, key=lambda x: x["code"]):
            g, c = r["gc_fasta"], r["ncbi_gc"]
            r["genome_gc"] = g if g is not None else c
            r["gc_src"] = "fasta" if g is not None else ("ncbi" if c is not None else "none")
            r["flag_diff"] = int(g is not None and c is not None and abs(g - c) > 0.01)
            fo.write("\t".join("" if r.get(k) is None else (f"{r[k]:.6f}" if isinstance(r[k], float) else str(r[k]))
                               for k in cols) + "\n")
    os.replace(f"{OUT}/genome_gc.tsv.tmp", f"{OUT}/genome_gc.tsv")
    src = [r["gc_src"] for r in res]
    print(f"species {len(res)}: fasta {src.count('fasta')}, ncbi {src.count('ncbi')}, none {src.count('none')}, "
          f"errors {sum(bool(r['error']) for r in res)}, flag_diff {sum(r['flag_diff'] for r in res)}")


if __name__ == "__main__":
    main()
