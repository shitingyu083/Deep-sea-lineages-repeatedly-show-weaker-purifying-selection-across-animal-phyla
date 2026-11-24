#!/bin/bash
# R3 fish genome-scale RER chain: 64b prep -> 36_rer_trees.sh (LG+G4, fixed fish topology, NJ 32)
# -> 37_rerconverge.R (lab, 1000 permulations, 24 cores) under memguard.        work/addon/R3/fishrer/RER.DONE
cd ~/deepsea_conv || exit 0
O=work/addon/R3/fishrer; mkdir -p $O/log
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
[ -s $O/src/locus_stats.tsv ] || envs/dc/bin/python scripts/64b_fishrer_prep_sty.py > $O/log/prep64b.txt 2>&1
# share the machine: start the gene trees only after the RTF1 HyPhy runs and the R2 pilot have finished
while [ ! -e work/addon/R3/rtf1/RTF1.DONE ] || [ ! -e work/addon/R2/pilot/PILOT.DONE ]; do sleep 120; done
bash scripts/memguard.sh 48 bash scripts/36_rer_trees.sh $O/src $O/topo_fish.tre $O/trees 32 > $O/log/trees.log 2>&1
[ -e $O/trees/TREES.DONE ] || { echo "trees failed" >> $O/log/chain.txt; exit 0; }
bash scripts/memguard.sh 48 envs/rphy/bin/Rscript scripts/37_rerconverge.R $O/trees/rer_trees.txt work/final_v3/date/lsd.timetree.nwk $O/lab lab 1000 24 > $O/log/rer.log 2>&1
[ -e $O/lab/DONE ] && touch $O/RER.DONE
exit 0
