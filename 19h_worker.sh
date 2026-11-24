#!/bin/bash
# Gene-tree worker session of the FINAL species tree (19h_final_tree.sh, 2026-10-19). Started in its own session,
# so it gets its own CPU share under autogroup scheduling. Runs the loci of W/cost.tsv in cost order (desc = most
# expensive first, asc = cheapest first) with the MAFFT threads of cost.tsv, NJ loci at a time; the per-locus flock
# of 19h_locus.sh keeps every locus to one run. The pid file W/workers/<pid>.pid tells the driver it is alive.
# Usage: cd ~/deepsea_conv; (setsid nohup bash scripts/19h_worker.sh W NJ desc|asc > W/log/worker_N.log 2>&1 < /dev/null &)
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
W=$(readlink -f $1); NJ=${2:-48}; ORDER=${3:-desc}
mkdir -p $W/workers
echo $$ > $W/workers/$$.pid
echo "$(date '+%F %T') worker $$ start W=$W NJ=$NJ order=$ORDER"
if [ $ORDER = asc ]; then tail -n +2 $W/cost.tsv | tac; else tail -n +2 $W/cost.tsv; fi \
  | cut -f1,5 | xargs -P $NJ -L1 bash scripts/19h_locus.sh $W
echo "$(date '+%F %T') worker $$ finished (xargs rc $?)"
rm -f $W/workers/$$.pid
