#!/usr/bin/env python3
"""Drop the species that are not tips of the labelling tree from a clade's codon alignments (2026-10-20).
Needed when a taxon left the species tree after the alignments were made, e.g. Euncla (Cnidaria): removed by
identity QC (its assembly is near-identical to Sibsp; work/draft_v2/fix.log.tsv, work/tree_v3/log/loci.log).
23_label_trees.py prunes the tree to the alignment's species, so such a species would leave HyPhy and CSUBST
with an alignment that does not match the tree.
Per work/codon/<clade>/aln/<hog>.fas with >= 1 such sequence: the original is copied to aln_pre43f/<hog>.fas and
the filtered alignment replaces it (codons that are gaps in every remaining sequence removed, frame kept); if
the remainder fails the 22_codon_aln.py retention rule (>= 4 species, >= 2 deep, >= 2 shallow, >= 300 nt) the
HOG leaves aln/ instead (original moved to aln_pre43f/). Any change removes trees/SOURCE_TREE.txt so that
43b_sel_launch.sh relabels. Idempotent. Run it before 23 / 24, never while 24 runs on the clade.
Log: aln_pre43f/sync43f.tsv (hog, dropped, n_left, nD, nS, len_nt, action); summary on stdout.
Usage: python 43f_sync_aln.py CLADE TREE [--dry]"""
import os
import shutil
import sys
import time

from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv")
sys.path.insert(0, f"{R}/scripts")
from importlib import import_module  # noqa: E402

LAB = import_module("21_hog_select").labels()
EMPTY = set("-?Nn")


def read_fa(p):
    seqs, name = {}, None
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def drop_empty_codons(keep):
    L = {len(x) for x in keep.values()}
    assert len(L) == 1 and L.pop() % 3 == 0, "ragged or out-of-frame alignment"
    n = len(next(iter(keep.values())))
    cod = [i for i in range(0, n, 3) if any(set(x[i:i + 3]) - EMPTY for x in keep.values())]
    return {s: "".join(x[i:i + 3] for i in cod) for s, x in keep.items()}


def main():
    clade, tree = sys.argv[1], sys.argv[2]
    dry = "--dry" in sys.argv[3:]
    W = f"{R}/work/codon/{clade}"
    A, B = f"{W}/aln", f"{W}/aln_pre43f"
    tips = set(Tree(tree, format=1).get_leaf_names())
    rows, n_seen = [], 0
    for fn in sorted(os.listdir(A)):
        if not fn.endswith(".fas"):
            continue
        n_seen += 1
        seqs = read_fa(f"{A}/{fn}")
        drop = sorted(set(seqs) - tips)
        if not drop:
            continue
        keep = drop_empty_codons({s: x for s, x in seqs.items() if s in tips})
        nD = sum(LAB.get(s) == "D" for s in keep)
        nS = sum(LAB.get(s) == "S" for s in keep)
        nt = len(next(iter(keep.values()))) if keep else 0
        ok = len(keep) >= 4 and nD >= 2 and nS >= 2 and nt >= 300
        rows.append((fn[:-4], ",".join(drop), len(keep), nD, nS, nt, "filtered" if ok else "removed"))
        if dry:
            continue
        os.makedirs(B, exist_ok=True)
        if os.path.exists(f"{B}/{fn}"):
            sys.exit(f"{B}/{fn} exists: refusing to overwrite a saved original")
        if ok:
            shutil.copy2(f"{A}/{fn}", f"{B}/{fn}")
            with open(f"{A}/{fn}.tmp43f", "w") as f:
                f.writelines(f">{s}\n{x}\n" for s, x in keep.items())
            os.replace(f"{A}/{fn}.tmp43f", f"{A}/{fn}")
        else:
            os.replace(f"{A}/{fn}", f"{B}/{fn}")
    if rows and not dry:
        log = f"{B}/sync43f.tsv"
        new = not os.path.exists(log)
        with open(log, "a") as f:
            if new:
                f.write("hog\tdropped\tn_left\tnD\tnS\tlen_nt\taction\ttree\ttime\n")
            stamp = time.strftime("%F %T")
            f.writelines("\t".join(map(str, r)) + f"\t{tree}\t{stamp}\n" for r in rows)
        src = f"{W}/trees/SOURCE_TREE.txt"
        if os.path.exists(src):
            os.remove(src)
            print(f"removed {src} (labels must be redone)")
    nf = sum(r[6] == "filtered" for r in rows)
    codes = sorted({c for r in rows for c in r[1].split(",")})
    print(f"{time.strftime('%F_%T')} 43f {clade}{' DRY' if dry else ''}: {n_seen} alignments; species not in "
          f"{tree}: {codes or 'none'}; filtered {nf}, removed {len(rows) - nf}")
    if dry and rows:
        nd = sorted(r[5] for r in rows)
        print(f"  len_nt after filter: min {nd[0]} median {nd[len(nd) // 2]}; removed HOGs: "
              f"{[r[0] for r in rows if r[6] == 'removed'][:20]}")


if __name__ == "__main__":
    main()
