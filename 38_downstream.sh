#!/bin/bash
# Downstream comparative analyses on one species tree: LSD2 dating (33_date_lsd2.py; the tree's own branch lengths
# kept, -blfix), then on the time tree ancestral depth habitat (31_asr_depth.R) and physicochemical PGLS
# (32_pgls_physchem.R) for each label set, all label sets in parallel. Rerun unchanged on the final tree.
# Usage: bash scripts/38_downstream.sh W TREE_OG ALN [SETS=all|main] [WAIT_FOR=file]
#   W        output dir (W/date, W/asr{,_strict,_1000}, W/pgls{,_strict,_1000}, W/log/38_*.log)
#   TREE_OG  species tree with branch lengths, outgroups Monbre,Salros included (rooted on them for dating)
#   ALN      the alignment the lengths were estimated from (site count for LSD2)
#   SETS     all = lab/asr, lab_strict/asr_strict, lab_1000/asr_1000; main = lab/asr only
#   WAIT_FOR start when this file exists (e.g. the BL.DONE marker of the tree builder)
# Calibrations whose clade is not monophyletic in TREE_OG are skipped (LSD_INTRUDERS=skip) and listed in
# W/date/calibrations_used.tsv.                                                               W/DOWNSTREAM.DONE
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
W=$1; T=$2; A=$3; SETS=${4:-all}; WF=${5:-}
RS=$PWD/envs/rphy/bin/Rscript
# single-threaded BLAS/OpenMP for the R steps (R's default 128-thread OpenBLAS pool stalls small repeated fits on the
# loaded server; their parallelism is NCORES workers); IQ-TREE keeps its -T threads
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p $W/log
log() { echo "$(date '+%F %T') [38 $(basename $W)] $*"; }
[ -e $W/DOWNSTREAM.DONE ] && exit 0
if [ -n "$WF" ]; then
  log "waiting for $WF"
  while [ ! -e $WF ]; do sleep 60; done
fi
[ -s $T ] && [ -s $A ] || { log "missing $T or $A"; exit 1; }

TT=$W/date/lsd.timetree.nwk
if [ ! -s $TT ]; then
  log "LSD2 dating of $T"
  LSD_BLFIX=1 LSD_MEM=64G LSD_INTRUDERS=skip bash scripts/memguard.sh 64 $PWD/envs/dc/bin/python \
    scripts/33_date_lsd2.py $T $A $W/date Monbre,Salros 16 > $W/log/38_date.log 2>&1
  rc=$?
  [ $rc = 0 ] && [ -s $TT ] || { log "dating failed (rc $rc)"; tail -5 $W/log/38_date.log; exit 1; }
  printf "tree\t%s\ttree_md5\t%s\n" $(readlink -f $T) $(md5sum < $T | cut -c1-32) > $W/date/input.tsv
fi
log "time tree $TT: $(grep -c . $W/date/dates.txt) calibrations; $(grep -i "root date\|tMRCA" $W/date/lsd.timetree.lsd 2>/dev/null | head -1)"

[ $SETS = main ] && sfx="main" || sfx="main strict 1000"
pids=""
for s in $sfx; do
  [ $s = main ] && { ac=asr; lc=lab; d=""; } || { ac=asr_$s; lc=lab_$s; d=_$s; }
  if [ ! -s $W/asr$d/transitions.tsv ]; then
    ASR_COL=$ac $RS scripts/31_asr_depth.R $TT $W/asr$d 1000 > $W/log/38_asr$d.log 2>&1 &
    pids="$pids $!"
  fi
  if [ ! -s $W/pgls$d/pgls_gene_summary.tsv ]; then
    LABEL_COL=$lc DATED=1 NCORES=16 bash scripts/memguard.sh 48 $RS scripts/32_pgls_physchem.R $TT $W/pgls$d \
      > $W/log/38_pgls$d.log 2>&1 &
    pids="$pids $!"
  fi
done
log "ASR and PGLS running ($SETS):$pids"
fail=0
for p in $pids; do wait $p || fail=$((fail + 1)); done
for s in $sfx; do
  [ $s = main ] && d="" || d=_$s
  [ -s $W/asr$d/transitions.tsv ] || { log "asr$d incomplete"; fail=$((fail + 1)); }
  [ -s $W/pgls$d/pgls_gene_summary.tsv ] || { log "pgls$d incomplete"; fail=$((fail + 1)); }
done
[ $fail = 0 ] || { log "$fail step(s) failed; see $W/log/38_*.log"; exit 1; }
touch $W/DOWNSTREAM.DONE
log "DOWNSTREAM.DONE"
