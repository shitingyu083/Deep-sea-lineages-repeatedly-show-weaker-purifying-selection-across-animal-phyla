#!/bin/bash
# Stages 2-3 of the fast DRAFT species tree (exec'd by 19f_fast_draft.sh); resumable via *.DONE.
#   2 TreeShrink per-species q = 0.05 on the FastTree gene trees; removals applied to the trees and the
#     filtered alignments -> W/final_aln, W/genetrees_shrunk.tre, W/partition.nex (LG+G4 per locus)  SHRINK.DONE
#   3 wASTRAL hybrid weighting, SH-like support range (-L); default search rounds (the -R two-step
#     search has a serial NNI phase: >15 min on 12 genes in the smoke test)                          WASTRAL.DONE
#   then exec scripts/19f_fast_bl.sh (branch lengths, rooting, gCF, RF; read fresh when it starts)
# Usage: bash scripts/19f_fast_post.sh W
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
P=~/anaconda3/envs/py310env/bin/python                 # TreeShrink lives here
W=$(readlink -f $1)
TOOLS="$PWD/envs/dc/bin/python scripts/19c_tools.py"   # ete3 3.1.3 lives here
FTOOLS="$PWD/envs/dc/bin/python scripts/19f_tools.py"
log() { echo "$(date '+%F %T') [post] $*"; }

if [ ! -e $W/SHRINK.DONE ]; then
  log "TreeShrink"
  ls $W/gt | grep '\.tre$' | sed 's/\.tre$//' > $W/genetrees.list
  while read b; do cat $W/gt/$b.tre; done < $W/genetrees.list > $W/genetrees.tre
  rm -rf $W/treeshrink
  $P ~/anaconda3/envs/py310env/bin/run_treeshrink.py -t $W/genetrees.tre -o $W/treeshrink -m per-species \
    -q 0.05 > $W/log/treeshrink.log 2>&1 || { log "treeshrink failed"; exit 1; }
  $FTOOLS shrink $W LG+G4 || { log "shrink failed"; exit 1; }
  touch $W/SHRINK.DONE
fi

if [ ! -e $W/WASTRAL.DONE ]; then
  log "wASTRAL"
  wastral -L --mode 1 -u 1 -t 32 --seed 20261019 -i $W/genetrees_shrunk.tre -o $W/wastral.tre \
    > $W/log/wastral.log 2>&1 || { log "wastral failed"; exit 1; }
  $TOOLS topo $W/wastral.tre $W/wastral_topo.tre
  touch $W/WASTRAL.DONE
fi

exec bash scripts/19f_fast_bl.sh $W
