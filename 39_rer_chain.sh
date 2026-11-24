#!/bin/bash
# RERconverge for one species tree: waits for the gene trees (36_rer_trees.sh, G/TREES.DONE) and the time tree
# (38_downstream.sh, written with date/input.tsv), then 37_rerconverge.R per label set, all sets in parallel,
# NCORES permulation workers each (memguard 48G). Rerun unchanged on the final tree.
# Usage: bash scripts/39_rer_chain.sh G TIMETREE OUT [NPERM=1000] [NCORES=32] [SETS=all|main]
#   G = 36_rer_trees.sh output dir; TIMETREE = W/date/lsd.timetree.nwk of 38_downstream.sh
#   OUT/{lab,lab_strict,lab_1000}/ (37 outputs), OUT/log/                                      OUT/RER.DONE
set -u
cd ~/deepsea_conv
G=$1; TT=$2; OUT=$3; NP=${4:-1000}; NC=${5:-32}; SETS=${6:-all}
RS=$PWD/envs/rphy/bin/Rscript
# single-threaded BLAS/OpenMP: RERconverge's many small per-gene fits stall in R's default 128-thread OpenBLAS pool on
# the loaded server; the parallelism is the NCORES permulation workers
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p $OUT/log
log() { echo "$(date '+%F %T') [39 $(basename $OUT)] $*"; }
[ -e $OUT/RER.DONE ] && exit 0
log "waiting for $G/TREES.DONE and $(dirname $TT)/input.tsv"
while [ ! -e $G/TREES.DONE ] || [ ! -s $(dirname $TT)/input.tsv ]; do sleep 60; done
printf "gene_trees\t%s\t%s\ntime_tree\t%s\t%s\n" $(readlink -f $G/rer_trees.txt) $(md5sum < $G/rer_trees.txt | cut -c1-32) \
  $(readlink -f $TT) $(md5sum < $TT | cut -c1-32) > $OUT/inputs.tsv
[ $SETS = main ] && sets="lab" || sets="lab lab_strict lab_1000"
pids=""
for lc in $sets; do
  [ -e $OUT/$lc/DONE ] && continue
  log "$lc: $NP permulations on $NC cores"
  ( bash scripts/memguard.sh 48 $RS scripts/37_rerconverge.R $G/rer_trees.txt $TT $OUT/$lc $lc $NP $NC \
      > $OUT/log/$lc.log 2>&1
    rc=$?
    [ $rc = 0 ] && [ -e $OUT/$lc/DONE ] && log "$lc done" || { log "$lc failed (rc $rc)"; tail -5 $OUT/log/$lc.log; exit 1; }
  ) &
  pids="$pids $!"
done
fail=0
for p in $pids; do wait $p || fail=1; done
[ $fail = 0 ] || exit 1
touch $OUT/RER.DONE
log "RER.DONE"
