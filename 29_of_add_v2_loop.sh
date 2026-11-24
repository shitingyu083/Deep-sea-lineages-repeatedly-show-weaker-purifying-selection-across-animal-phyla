#!/bin/bash
# Add v2 panel species to the 4 clades whose OrthoFinder already ran on v1 (reuses v1 BLASTs via -b).
# Waits for v1 .done and for all v2 clade members to be QC'd (qc_table refreshed by 20_orthofinder_loop).
# Usage: setsid nohup bash scripts/29_of_add_v2_loop.sh 48 8 > logs/of_add_v2.log 2>&1 < /dev/null &
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH
P=~/anaconda3/envs/py310env/bin/python
T=${1:-48}; A=${2:-8}
CLADES="Annelida Arthropoda Cnidaria Echinodermata"
while true; do
  # refresh QC table ourselves only once the main OrthoFinder loop (which does it) has exited
  [ -f work/orthofinder/ORTHOFINDER.DONE ] && $P scripts/15_qc_table.py > /dev/null
  left=0
  for c in $CLADES; do
    d=work/orthofinder/$c
    [ -f $d/.v2_started ] && continue
    left=$((left+1))
    [ -f $d/.done ] || continue
    $P scripts/29_of_add_v2.py $c; rc=$?
    if [ $rc = 3 ]; then touch $d/.v2_started $d/.v2_done; continue; fi
    [ $rc = 0 ] || continue
    wd=$(ls -d $d/of/Results_*/WorkingDirectory 2>/dev/null | while read w; do
           ls $w/Blast0_0.txt* >/dev/null 2>&1 && echo $w; done | head -1)
    [ -n "$wd" ] || { echo "$(date +%F_%T) $c no v1 WorkingDirectory with BLASTs"; touch $d/.v2_started $d/.v2_fail; continue; }
    touch $d/.v2_started
    echo "$(date +%F_%T) add_v2 $c ($(ls $d/add_v2 | wc -l) species) on $wd"
    ( orthofinder -b $wd -f $d/add_v2 -t $T -a $A -X > $d/of_add_v2.log 2>&1 \
        && touch $d/.v2_done || touch $d/.v2_fail ) &
  done
  [ $left = 0 ] && break
  sleep 600
done
wait
echo DONE > work/orthofinder/ADD_V2.DONE
