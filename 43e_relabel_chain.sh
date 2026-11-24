#!/bin/bash
# 43e_relabel_chain.sh (2026-10-20): move a clade whose selection run (43b -> 24) started on the DRAFT species tree
# onto the FINAL tree without stopping that run.
#   1. wait for $(dirname FINAL)/BL.DONE (the rooted FINAL tree of 19h_post.sh)
#   2. 43d_relabel.py --apply: HOGs whose alignment-induced rooted topology changed get 23's labels on FINAL and
#      their draft-tree outputs go to relabel43d/stale/ (in-flight jobs of those HOGs are drained first; jobs the
#      first pass dispatches later run on the FINAL labels). If the drain fails: retry once after SEL.DONE.
#   3. 43g_topup.py: HyPhy tiers on the FINAL-eligible frame (tiers/hylist_final.txt, nullhy_final.txt)
#   4. after the first pass (SEL.DONE): 24 again with TAG re43d on the final tiers, at autogroup nice 0; it runs
#      only what is missing (relabelled HOGs, jobs skipped while their label files were away, top-up HOGs).
# Restartable (each step is skipped when its marker exists). Start it in its own session (setsid).
# Usage: bash scripts/43e_relabel_chain.sh CLADE DRAFT FINAL [NJ=48]   (log logs/relabel43_<clade>.log)
R=~/deepsea_conv; cd $R || exit 1
C=$1; D=$2; F=$3; NJ=${4:-48}; PY=envs/dc/bin/python
W=work/codon/$C; O=$W/relabel43d; L=logs/relabel43_$C.log; B=$(dirname $F)/BL.DONE
log() { echo "$(date +%F_%T) [43e $C] $*" >> $L; }
mkdir -p work/chain43; echo $$ > work/chain43/relabel_$C.pid
log "start (pid $$, autogroup $(tr -s ' ' < /proc/self/autogroup)): DRAFT $D, FINAL $F, NJ $NJ"
[ -s $D ] && [ -s $W/trees/SOURCE_TREE.txt ] || { log "missing $D or $W/trees/SOURCE_TREE.txt"; exit 1; }
while [ ! -e $B ]; do sleep 300; done
[ -s $F ] || { log "$B present but $F missing or empty"; exit 1; }
log "$B present; FINAL md5 $(md5sum < $F | cut -c1-32)"
if [ ! -e $O/APPLIED ]; then
  $PY scripts/43d_relabel.py $C $D $F --apply >> $L 2>&1; rc=$?
  if [ $rc = 1 ] && [ -e $O/APPLY_FAILED ] && [ ! -e $O/APPLIED ] && [ ! -e $O/APPLYING ] && [ ! -e $O/APPLY_ERROR ]
  then
    log "apply failed (jobs of changed HOGs did not drain); retry once after $W/SEL.DONE"
    while [ ! -e $W/SEL.DONE ]; do sleep 300; done
    $PY scripts/43d_relabel.py $C $D $F --apply >> $L 2>&1; rc=$?
  fi
  [ $rc = 0 ] && [ -e $O/APPLIED ] || { log "43d --apply rc $rc: chain stopped"; exit 1; }
fi
$PY scripts/43g_topup.py $C >> $L 2>&1 || { log "43g failed: chain stopped"; exit 1; }
[ -s $O/changed.txt ] || { log "no HOG changed: the first pass is final"; exit 0; }
[ -e $W/SEL_re43d.DONE ] && { log "SEL_re43d.DONE exists"; exit 0; }
log "waiting for $W/SEL.DONE (first pass)"
while [ ! -e $W/SEL.DONE ]; do sleep 300; done
for i in $(seq 10); do echo 0 > /proc/self/autogroup 2>/dev/null; grep -q 'nice 0$' /proc/self/autogroup && break; sleep 0.3; done
log "second pass: 24 TAG re43d on tiers/hylist_final.txt (autogroup $(tr -s ' ' < /proc/self/autogroup))"
TAG=re43d HYLIST=$R/$W/tiers/hylist_final.txt NULLLIST=$R/$W/tiers/hylist_final.txt \
  NULLHY=$R/$W/tiers/nullhy_final.txt bash scripts/24_selection_run.sh $C $NJ 5 1 >> logs/sel_$C.log 2>&1
rc=$?
log "second pass exit $rc$([ -e $W/SEL_re43d.DONE ] && echo '; SEL_re43d.DONE')"
exit $rc
