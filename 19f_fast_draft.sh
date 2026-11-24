#!/bin/bash
# Fast DRAFT species tree on the 396 QC-pass species of work/tree_v2 (631 BUSCO loci), so that dating,
# ancestral depth, PGLS and RERconverge can run while the panel is being completed. Everything downstream is
# re-run on the final tree (19c_species_tree.sh on the completed panel: L-INS-i, ModelFinder + aBayes gene
# trees, wASTRAL + concatenated ML, gCF/sCF).
#   1 gene trees (19f_fast_locus.sh, NJ in parallel)                                        GENETREES.DONE
#   2-4 19f_fast_post.sh: TreeShrink, wASTRAL, partitioned branch lengths, rooting, gCF      DRAFT_BL.DONE
# Usage: setsid nohup bash scripts/19f_fast_draft.sh work/tree_v2_fast 32 >> logs/tree_v2_fast.log 2>&1 < /dev/null &
#   extra session sharing the locus list (smallest loci first, exits when the list is done):
#        setsid nohup bash scripts/19f_fast_draft.sh work/tree_v2_fast 32 worker >> logs/tree_v2_fast_worker.log 2>&1 < /dev/null &
set -u
shopt -s nullglob
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
W=$(readlink -f $1); NJ=${2:-32}; ROLE=${3:-main}
LOCI=work/tree_v2/loci
log() { echo "$(date '+%F %T') [$ROLE] $*"; }
mkdir -p $W/aln $W/trim $W/filt $W/gt $W/skip $W/log $W/claim
run_pass() {  # $1 = ls order: -S largest first (main), -Sr smallest first (worker)
  ls $1 $LOCI | sed 's/\.faa$//' | xargs -P $NJ -I{} bash scripts/19f_fast_locus.sh $W {}
}
if [ "$ROLE" = worker ]; then
  log "worker start NJ=$NJ"; run_pass -Sr; log "worker finished"; exit 0
fi
if [ ! -e $W/GENETREES.DONE ]; then
  log "gene trees, NJ=$NJ"
  run_pass -S
  while :; do              # loci still held by a worker session
    held=0
    for f in $W/claim/*.lock; do flock -n $f true || held=$((held+1)); done
    [ $held -eq 0 ] && break
    log "waiting for $held loci held by other sessions"; sleep 60
  done
  run_pass -S              # retry anything that failed once
  a=$(ls $LOCI | wc -l); n=$(ls $W/gt | grep -c '\.tre$'); s=$(ls $W/skip | wc -l)
  ls $LOCI | sed 's/\.faa$//' | while read b; do [ -s $W/gt/$b.tre ] || [ -e $W/skip/$b ] || echo $b; done > $W/failed_loci.txt
  log "gene trees $n, skipped $s, failed $(wc -l < $W/failed_loci.txt), of $a loci"
  [ $n -ge $((a / 2)) ] || { log "too few gene trees"; exit 1; }
  touch $W/GENETREES.DONE
fi
exec bash scripts/19f_fast_post.sh $W
