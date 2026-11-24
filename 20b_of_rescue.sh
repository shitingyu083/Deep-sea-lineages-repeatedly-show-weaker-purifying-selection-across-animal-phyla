#!/bin/bash
# Rescue OrthoFinder runs that failed in the orthologue stage (numpy-2 scalar bug in
# scripts_of/orthologues.py:263, patched 2026-10-17; runs started before the patch still fail).
# If orthogroups were finished, resume from them with -fg; on success mark .done.
# Usage: nohup bash scripts/20b_of_rescue.sh 48 8 > logs/of_rescue.log 2>&1 &
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH
T=${1:-48}; A=${2:-8}
CLADES="Actinopteri Chondrichthyes Mollusca Annelida Arthropoda Cnidaria Porifera Echinodermata"
while true; do
  left=0
  for c in $CLADES; do
    d=work/orthofinder/$c
    [ -f $d/.done ] || [ -f $d/.rescue_fail ] && continue
    left=$((left+1))
    [ -f $d/.fail ] || continue
    res=$(ls -d $d/of/Results_* 2>/dev/null | while read x; do [ -s $x/Orthogroups/Orthogroups.tsv ] && echo $x; done | tail -1)
    if [ -z "$res" ]; then
      echo "$(date +%F_%T) $c failed before orthogroups; not rescuable here"; touch $d/.rescue_fail; continue
    fi
    echo "$(date +%F_%T) rescue $c from $res"
    if nice -n 19 orthofinder -fg $res -t $T -a $A -X > $d/of_resume.log 2>&1; then
      mv $d/.fail $d/.fail_resumed; touch $d/.done; echo "$(date +%F_%T) rescued $c"
    else
      touch $d/.rescue_fail; echo "$(date +%F_%T) rescue FAILED $c"
    fi
  done
  [ $left = 0 ] && break
  sleep 300
done
echo RESCUE_DONE
