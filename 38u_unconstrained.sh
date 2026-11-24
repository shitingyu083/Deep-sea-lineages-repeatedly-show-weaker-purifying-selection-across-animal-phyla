#!/bin/bash
# Topology sensitivity of the draft analyses: the unconstrained wASTRAL draft (work/tree_v2_fast, no topology fixes)
# with the taxa dropped in meta/topo_fix.tsv removed as in the fixed tree, analysed exactly like the fixed draft:
# gene-tree branch lengths re-estimated on this topology (36_rer_trees.sh), then 38c_chain.sh with the main labels
# (dating, ASR, PGLS, RERconverge).
# Usage: bash scripts/38u_unconstrained.sh [W=work/draft_v2_unc] [R=work/rer_v2_unc]           W/CHAIN.DONE
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
U=work/tree_v2_fast; W=${1:-work/draft_v2_unc}; R=${2:-work/rer_v2_unc}
mkdir -p $W/log
log() { echo "$(date '+%F %T') [38u] $*"; }
[ -e $W/CHAIN.DONE ] && exit 0
log "waiting for $U/DRAFT_BL.DONE"
while [ ! -e $U/DRAFT_BL.DONE ]; do sleep 60; done
if [ ! -s $W/tree_og.nwk ]; then
  python - $U/species_tree_astral_og.nwk meta/topo_fix.tsv $W/tree_og.tmp <<'EOF' || { log "pruning failed"; exit 1; }
import csv, sys
from ete3 import Tree
t = Tree(open(sys.argv[1]).read().strip(), format=1, quoted_node_names=True)
drop = [r["subject"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t") if r["op"] == "drop"]
drop = [x for x in drop if x in set(t.get_leaf_names())]
t.prune([x for x in t.get_leaf_names() if x not in drop], preserve_branch_length=True)
t.write(outfile=sys.argv[3], format=5)
print("dropped:", ",".join(drop) or "-", "; tips", len(t))
EOF
  mv $W/tree_og.tmp $W/tree_og.nwk
fi
if [ ! -s $W/concat.faa ]; then                  # the concatenated alignment on the same taxa (LSD2 site count)
  python - $U/concat.faa $W/tree_og.nwk $W/concat.tmp <<'EOF' || { log "alignment subset failed"; exit 1; }
import sys
from ete3 import Tree
keep = set(Tree(open(sys.argv[2]).read().strip(), format=1, quoted_node_names=True).get_leaf_names())
n = k = 0; w = False
with open(sys.argv[1]) as f, open(sys.argv[3], "w") as o:
    for line in f:
        if line.startswith(">"):
            n += 1; w = line[1:].split()[0] in keep; k += w
        if w: o.write(line)
print("alignment: kept", k, "of", n, "sequences")
sys.exit(0 if k == len(keep) else "tree taxa missing from the alignment")
EOF
  mv $W/concat.tmp $W/concat.faa
fi
bash scripts/36_rer_trees.sh $U $W/tree_og.nwk $R 48 > $W/log/36.log 2>&1 || { log "RER gene trees failed"; exit 1; }
bash scripts/38c_chain.sh $W $W/tree_og.nwk $W/concat.faa $R $R/rer main || { log "chain failed"; exit 1; }
log "done"
