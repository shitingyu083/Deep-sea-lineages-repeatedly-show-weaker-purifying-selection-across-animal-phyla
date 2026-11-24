#!/bin/bash
# FINAL species tree (2026-10-19) on the frozen panel: the QC-pass species of meta/qc_table.tsv (rule of
# 15_qc_table.py: protein-mode compleasm C >= 60 %, M <= 25 %) minus the taxa dropped in meta/topo_fix.tsv, plus the
# choanoflagellate outgroups Monbre and Salros.
#   0 loci: 18c_busco_loci_v2.py (metazoa_odb12 single-copy BUSCOs present in >= 50 % of the ingroup); expected
#     MAFFT cost per locus (19h_tools.py cost) -> W/cost.tsv = run order and MAFFT threads              LOCI.DONE
#   1 gene trees (19h_locus.sh): this driver and NW worker sessions (own sessions = own CPU shares under autogroup
#     scheduling), NJ loci each; all go most expensive first except the last worker (cheapest first); one retry of
#     unfinished loci; stops without the marker if > 5 % of the loci fail or W/CUT exists         GENETREES.DONE
#   then the autogroup nice of my other sessions is restored (19h_ag.py, W/autogroup.tsv) and 19h_post.sh runs
#   stages 2-7 (TreeShrink, wASTRAL, fixed topology, concatenated ML, branch lengths, downstream chains).
# One driver per W (flock on W/driver.lock, kept by 19h_post.sh). Env: NOPOST=1 stops after the gene trees (smoke
# test); F and RV are passed on to 19h_post.sh.
# Usage: cd ~/deepsea_conv; (setsid nohup bash scripts/19h_final_tree.sh W [NJ=48] [NW=3] > W.log 2>&1 < /dev/null &)
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
mkdir -p $1; W=$(readlink -f $1); NJ=${2:-48}; NW=${3:-3}
PY=$PWD/envs/dc/bin/python
log() { echo "$(date '+%F %T') [19h] $*"; }
exec 8> $W/driver.lock
flock -n 8 || { log "another driver holds $W/driver.lock"; exit 1; }
mkdir -p $W/aln $W/trim $W/filt $W/gt $W/skip $W/log $W/workers $W/claim
log "start W=$W NJ=$NJ NW=$NW (pid $$)"

# 0 loci
if [ ! -e $W/LOCI.DONE ]; then
  EX=$($PY scripts/19g_tools.py dropped meta/topo_fix.tsv) || { log "drop list failed"; exit 1; }
  log "loci: STATUSES=pass EXCLUDE=$EX"
  STATUSES=pass EXCLUDE=$EX $PY scripts/18c_busco_loci_v2.py $W 0.5 > $W/log/loci.log 2>&1 \
    || { log "18c failed"; tail -n 20 $W/log/loci.log; exit 1; }
  $PY scripts/19h_tools.py cost $W >> $W/log/loci.log 2>&1 || { log "cost failed"; tail -n 5 $W/log/loci.log; exit 1; }
  tail -n 4 $W/log/loci.log
  touch $W/LOCI.DONE
fi

# 1 gene trees
unfinished() { tail -n +2 $W/cost.tsv | cut -f1 | while read b; do
  [ -s $W/gt/$b.treefile ] || [ -e $W/skip/$b ] || echo $b; done; }
alive() { local f p n=0; for f in $W/workers/*.pid; do [ -e "$f" ] || continue; p=$(cat $f)
    if tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -q 19h_worker.sh; then n=$((n + 1)); else rm -f $f; fi
  done; echo $n; }
held() { unfinished | while read b; do flock -n $W/claim/$b.lock true 2>/dev/null || echo $b; done | wc -l; }
if [ ! -e $W/GENETREES.DONE ]; then
  [ -e $W/CUT ] && { log "W/CUT present: not starting"; exit 0; }
  nl=$(($(wc -l < $W/cost.tsv) - 1))
  log "gene trees: $nl loci, $(unfinished | wc -l) to do; driver + $NW worker sessions x $NJ loci"
  if [ $(alive) = 0 ]; then
    for i in $(seq 1 $NW); do
      o=desc; [ $i = $NW ] && o=asc
      (setsid nohup bash scripts/19h_worker.sh $W $NJ $o > $W/log/worker_$i.log 2>&1 < /dev/null 8>&- &)
      sleep 1
    done
  else
    log "$(alive) workers already alive, none started"
  fi
  tail -n +2 $W/cost.tsv | cut -f1,5 | xargs -P $NJ -L1 bash scripts/19h_locus.sh $W 8>&-
  log "driver pass done; waiting for the workers and for running loci"
  while [ $(alive) -gt 0 ] || [ $(held) -gt 0 ]; do sleep 120; done
  [ -e $W/CUT ] && { log "W/CUT present: stopping (no GENETREES.DONE)"; exit 0; }
  unfinished > $W/retry.txt
  if [ -s $W/retry.txt ]; then
    log "retrying $(wc -l < $W/retry.txt) unfinished loci"
    awk -F'\t' 'NR == FNR {r[$1]; next} FNR > 1 && ($1 in r) {print $1 "\t" $5}' $W/retry.txt $W/cost.tsv \
      | xargs -P $NJ -L1 bash scripts/19h_locus.sh $W 8>&-
    while [ $(held) -gt 0 ]; do sleep 120; done
  fi
  unfinished > $W/failed_loci.txt
  nf=$(wc -l < $W/failed_loci.txt)
  $PY scripts/19h_tools.py times $W
  log "gene trees: $(ls $W/gt | grep -c '\.treefile$') trees, $(ls $W/skip | wc -l) skipped, $nf failed"
  [ $((nf * 20)) -gt $nl ] && { log "more than 5 % of the loci failed: stopping (no GENETREES.DONE)"; exit 1; }
  touch $W/GENETREES.DONE
fi
[ -s $W/autogroup.tsv ] && $PY scripts/19h_ag.py restore $W
[ -n "${NOPOST:-}" ] && { log "NOPOST: stopping after the gene trees"; exit 0; }
while [ ! -s scripts/19h_post.sh ]; do log "waiting for scripts/19h_post.sh"; sleep 600; done
exec bash scripts/19h_post.sh $W ${F:-work/final_v3} ${RV:-work/rer_v3}
