#!/bin/bash
# Tiered selection run for one clade (2026-10-20): 43f_sync_aln.py on TREE (drops alignment species absent
# from TREE, e.g. Euncla) -> 23_label_trees.py on TREE (5 shallow-shallow null sets)
# -> 43c_sel_lists.py (hylist = 1000 random eligible HOGs, nullhy = its first 400) ->
# 24_selection_run.sh v3 (NJ slots, K 5 CSUBST null sets on hylist, HyPhy on null set 1 of nullhy).
# Start it in its own session (setsid): it sets that session's autogroup nice to AGN and records its
# PID (= session id) in work/chain43/sel_<clade>.pid. Log logs/sel_<clade>.log.
# Usage: bash scripts/43b_sel_launch.sh CLADE TREE [NJ=48] [AGN=10]
R=~/deepsea_conv; cd $R
C=$1; T=$2; NJ=${3:-48}; AGN=${4:-10}
echo $AGN > /proc/self/autogroup 2>/dev/null
mkdir -p work/chain43; echo $$ > work/chain43/sel_$C.pid
W=work/codon/$C; L=logs/sel_$C.log
log() { echo "$(date +%F_%T) $*" >> $L; }
[ -e $W/PREP.DONE ] || { log "no $W/PREP.DONE"; exit 1; }
[ -s $T ] || { log "no tree $T"; exit 1; }
envs/dc/bin/python scripts/43f_sync_aln.py $C $T >> $L 2>&1 || { log "43f failed (rc $?)"; exit 1; }
if [ -s $W/trees/SOURCE_TREE.txt ] && [ "$(cut -f2 $W/trees/SOURCE_TREE.txt)" = "$(md5sum < $T | cut -c1-32)" ]; then
  log "labels already from $T"
else
  log "23_label_trees.py on $T"
  envs/dc/bin/python scripts/23_label_trees.py $C $T --null-from-shallow 5 > $W/label23.log 2>&1 \
    || { log "23 failed (rc $?, $W/label23.log)"; exit 1; }
  printf "%s\t%s\n" $T $(md5sum < $T | cut -c1-32) > $W/trees/SOURCE_TREE.txt
  log "23 done: $(($(wc -l < $W/lineages.tsv) - 1)) HOGs; missing-in-tree warnings: $(grep -c 'missing in species tree' $W/label23.log)"
fi
envs/dc/bin/python scripts/43c_sel_lists.py $C >> $L 2>&1 || { log "43c failed"; exit 1; }
log "24 start (NJ $NJ, autogroup nice $AGN)"
HYLIST=$R/$W/tiers/hylist.txt NULLLIST=$R/$W/tiers/hylist.txt NULLHY=$R/$W/tiers/nullhy.txt \
  bash scripts/24_selection_run.sh $C $NJ 5 1 >> $L 2>&1
log "24 exit $?"
