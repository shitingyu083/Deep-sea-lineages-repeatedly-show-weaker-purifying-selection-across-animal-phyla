#!/usr/bin/env python3
"""Standardise the lab's local mollusc proteomes/CDS into data/proteomes/.

For each species: pair pep<->cds records (by order, verified by translation;
fallback by normalised ID), drop pairs whose CDS translation disagrees with the
peptide, collapse isoforms where the ID pattern encodes the gene (.tN), keep
the longest isoform, and rename to CODE_000001 with a mapping table.
Outputs: data/proteomes/<CODE>.faa, data/cds/<CODE>.fna, meta/idmap/<CODE>.tsv,
meta/local_prep_report.tsv
"""
import gzip, os, re, sys

HOME = os.path.expanduser("~")
SRC = f"{HOME}/deepsea"
OUT = f"{HOME}/deepsea_conv"

# code, species, role, pep, cds
# Adu (pred2 exon models, avg 74 aa, BUSCO 1%) -> miniprot re-annotation;
# Gae/Mme/Hru/Lgi -> NCBI RefSeq via 08_prep_ncbi_annot.py (isoform-aware)
LOCAL = [
    ("Ama", "Archivesica marissinica", "deep", "allpep/Ama_pep.fasta", "allcds/Ama_cds.fasta"),
    ("Bla", "Bathyacmaea lactea", "deep", "allpep/Bla_pep.fasta", "allcds/Bla_cds.fasta"),
    ("Bpl", "Gigantidas platifrons", "deep", "allpep/Bpl_pep.fasta", "allcds/Bpl_cds.fasta"),
    ("Cbi", "Conchocele bisecta", "deep", "allpep/Cbi_pep.fasta", "allcds/Cbi_cds.fasta"),
    ("Csq", "Chrysomallon squamiferum", "deep", "allpep/Csq_pep.fasta", "allcds/Csq_cds.fasta"),
    ("Xdo", "Xylophaga dorsalis", "deep", "D_ Xylophaga dorsalis/Xte.pep.fasta", "allcds/Xte_cds.fasta"),
    ("Hha", "Hanleya hanleyi", "shallow", "allpep/Hha.pep.fasta", "wq_Hanleya_hanleyi/Hanleya_hanleyi.cds.fa"),
    ("Mco", "Mytilus coruscus", "shallow", "allpep/Mco_pep.fasta", "allcds/Mco_cds.fasta"),
    ("Rve", "Rapana venosa", "shallow", "allpep/Rve_pep.fasta", "allcds/Rve_cds.fasta"),
    ("Sco", "Sinonovacula constricta", "shallow", "allpep/Sco_pep.fasta", "allcds/Sco_cds.fasta"),
    ("Sph", "Sepia pharaonis", "shallow", "allpep/Sph_pep.fasta", "allcds/Sph_cds.fasta"),
    ("Tna", "Teredo navalis", "shallow", "allpep/Tna.pep.fasta", "allcds/Tna_cds.fasta"),
]

CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aas[i]; i += 1


def translate(s):
    s = s.upper()
    return "".join(CODON.get(s[j:j + 3], "X") for j in range(0, len(s) - 2, 3))


def read_fa(path):
    op = gzip.open if path.endswith(".gz") else open
    recs, h, buf = [], None, []
    with op(path, "rt") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(buf)))
                h, buf = line[1:].strip(), []
            elif line:
                buf.append(line.strip())
    if h is not None:
        recs.append((h, "".join(buf)))
    return recs


def norm_id(h, code):
    t = h.split()[0]
    if t.startswith(code + "_"):
        t = t[len(code) + 1:]
    t = re.sub(r"_cds$", "", t)
    return t


def agree(pep, cds):
    p = pep.replace("*", "").upper()
    for off in (0,):  # CDS should start in frame
        q = translate(cds[off:]).rstrip("*")
        if q.replace("*", "") == p or q == p:
            return True
        # tolerate start-codon/terminal X and ambiguous bases
        if len(q) == len(p) and sum(x != y for x, y in zip(q, p)) <= max(2, len(p) // 100):
            return True
    return False


def gene_key(pid, code):
    # only Hha (BRAKER jgN.tN) encodes isoforms; Mco ContigN.tN are distinct genes
    m = re.match(r"^(jg\d+)\.t\d+$", pid) if code == "Hha" else None
    return m.group(1) if m else pid


def main():
    for d in ("data/proteomes", "data/cds", "meta/idmap"):
        os.makedirs(f"{OUT}/{d}", exist_ok=True)
    rep = open(f"{OUT}/meta/local_prep_report.tsv", "w")
    rep.write("code\tspecies\trole\tn_pep\tn_cds\tpairing\tn_agree\tn_genes_kept\n")
    for code, sp, role, pf, cf in LOCAL:
        pep = read_fa(f"{SRC}/{pf}")
        cds = read_fa(f"{SRC}/{cf}")
        pairs, mode = [], "order"
        if len(pep) == len(cds):
            pairs = list(zip(pep, cds))
            ok = sum(agree(p[1], c[1]) for p, c in pairs[:500])
            if ok < 0.9 * min(500, len(pairs)):
                pairs = []
        if not pairs:
            mode = "id"
            cmap = {norm_id(h, code): (h, s) for h, s in cds}
            for h, s in pep:
                k = norm_id(h, code)
                if k not in cmap:
                    k2 = re.sub(r"(\.\d+)+$", "", k)
                    k = k2 if k2 in cmap else None
                    if k is None:
                        for kk in (norm_id(h, code) + ".1.1",):
                            if kk in cmap:
                                k = kk
                if k is not None and k in cmap:
                    pairs.append(((h, s), cmap[k]))
        good = [(p, c) for p, c in pairs if agree(p[1], c[1])]
        best = {}
        for (ph, ps), (ch, cs) in good:
            pid = norm_id(ph, code)
            g = gene_key(pid, code)
            L = len(ps.replace("*", ""))
            if L < 30:
                continue
            if g not in best or L > best[g][0]:
                best[g] = (L, pid, ps, cs)
        fa = open(f"{OUT}/data/proteomes/{code}.faa", "w")
        fn = open(f"{OUT}/data/cds/{code}.fna", "w")
        mp = open(f"{OUT}/meta/idmap/{code}.tsv", "w")
        for n, g in enumerate(sorted(best), 1):
            L, pid, ps, cs = best[g]
            nid = f"{code}_{n:06d}"
            p = ps.replace("*", "")
            c = cs.upper()
            if translate(c[-3:]) == "*":
                c = c[:-3]
            fa.write(f">{nid}\n{p}\n")
            fn.write(f">{nid}\n{c}\n")
            mp.write(f"{nid}\t{pid}\t{g}\n")
        for h in (fa, fn, mp):
            h.close()
        rep.write(f"{code}\t{sp}\t{role}\t{len(pep)}\t{len(cds)}\t{mode}\t{len(good)}\t{len(best)}\n")
        rep.flush()
        print(code, len(pep), len(cds), mode, len(good), len(best), flush=True)
    rep.close()


if __name__ == "__main__":
    main()
