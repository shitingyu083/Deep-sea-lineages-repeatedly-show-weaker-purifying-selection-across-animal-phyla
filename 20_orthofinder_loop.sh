#!/bin/bash
# Launch OrthoFinder 2.5.5 for each clade as soon as all its members are QC'd.
# (Run 15_qc_table.py each cycle so readiness is current.) Each clade: -t T -a A.
# Usage: nohup bash scripts/20_orthofinder_loop.sh 48 8 > logs/orthofinder_loop.log 2>&1 &
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH
P=~/anaconda3/envs/py310env/bin/python
T=${1:-48}; A=${2:-8}
CLADES="Actinopteri Chondrichthyes Mollusca Annelida Arthropoda Cnidaria Porifera Echinodermata"
while true; do
  $P scripts/15_qc_table.py > /dev/null
  left=0
  for c in $CLADES; do
    d=work/orthofinder/$c
    [ -f $d/.started ] && continue
    left=$((left+1))
    if $P scripts/16_orthofinder_setup.py --clade=$c >> logs/orthofinder_loop.log 2>&1; then
      touch $d/.started
      echo "$(date +%F_%T) start $c ($(ls $d/input | wc -l) species)"
      ( orthofinder -f $d/input -t $T -a $A -X -o $d/of > $d/of.log 2>&1 \
          && touch $d/.done || touch $d/.fail ) &
    fi
  done
  [ $left = 0 ] && break
  sleep 600
done
wait
echo DONE > work/orthofinder/ORTHOFINDER.DONE
