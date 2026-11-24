#!/bin/bash
# All downstream analyses of one species tree: 38_downstream.sh (LSD2 dating -> ASR + PGLS) and, as soon as the
# time tree is written, 39_rer_chain.sh (RERconverge on the gene trees in G) alongside ASR/PGLS.
# Rerun unchanged on the final tree.
# Usage: [NPERM=1000] [NCORES=32] bash scripts/38c_chain.sh W TREE_OG ALN G RER_OUT [SETS=all|main] [WAIT_FOR=file]
#   arguments as in 38_downstream.sh (W TREE_OG ALN SETS WAIT_FOR) and 39_rer_chain.sh (G, RER_OUT)   W/CHAIN.DONE
set -u
cd ~/deepsea_conv
W=$1; T=$2; A=$3; G=$4; RO=$5; SETS=${6:-all}; WF=${7:-}
mkdir -p $W/log
log() { echo "$(date '+%F %T') [38c $(basename $W)] $*"; }
[ -e $W/CHAIN.DONE ] && exit 0
bash scripts/38_downstream.sh $W $T $A $SETS $WF >> $W/log/38.log 2>&1 &
p=$!
while [ ! -s $W/date/input.tsv ] && kill -0 $p 2>/dev/null; do sleep 60; done
[ -s $W/date/input.tsv ] || { log "dating failed, see $W/log/38.log"; exit 1; }
log "time tree ready; RERconverge on $G"
bash scripts/39_rer_chain.sh $G $W/date/lsd.timetree.nwk $RO ${NPERM:-1000} ${NCORES:-32} $SETS >> $W/log/39.log 2>&1
r39=$?
wait $p; r38=$?
[ $r38 = 0 ] && [ $r39 = 0 ] || { log "failed: 38 rc $r38, 39 rc $r39"; exit 1; }
touch $W/CHAIN.DONE
log "CHAIN.DONE"
