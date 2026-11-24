#!/bin/bash
# One locus of 19c_species_tree.sh: MAFFT L-INS-i -> trimAl -automated1 -> fragment filter -> IQ-TREE
# (ModelFinder, aBayes); every step resumable.
# 2026-10-19: extra worker sessions (19c_worker.sh) may run loci of the same W. Each locus runs under flock on
#   W/claim/<locus>.lock (released when the whole process tree ends, also after a crash); a held locus is
#   skipped. W/CUT stops new loci from starting (time-boxed draft); loci already running finish.
W=$1; b=$2; T=${3:-3}; P=~/anaconda3/envs/py310env/bin/python; S=~/deepsea_conv/scripts
[ -s $W/gt/$b.treefile ] && exit 0
[ -e $W/skip/$b ] && exit 0
[ -e $W/CUT ] && exit 0
mkdir -p $W/claim
exec 9> $W/claim/$b.lock
flock -n 9 || exit 0
[ -s $W/gt/$b.treefile ] && exit 0   # finished by another run between the first check and the lock
[ -e $W/skip/$b ] && exit 0
if [ ! -s $W/aln/$b.faa ]; then
  mafft --localpair --maxiterate 1000 --thread $T --quiet $W/loci/$b.faa > $W/aln/$b.tmp 2> $W/log/$b.mafft \
    && [ -s $W/aln/$b.tmp ] && mv $W/aln/$b.tmp $W/aln/$b.faa || { echo "$b mafft failed"; exit 1; }
fi
if [ ! -s $W/filt/$b.faa ]; then
  trimal -in $W/aln/$b.faa -out $W/trim/$b.faa -automated1 > $W/log/$b.trimal 2>&1 || { echo "$b trimal failed"; exit 1; }
  $P $S/19c_tools.py frag $W/trim/$b.faa $W/filt/$b.faa 0.5 100 4; rc=$?
  if [ $rc = 2 ]; then echo "too short/few taxa after trimming" > $W/skip/$b; exit 0; fi
  [ $rc = 0 ] || { echo "$b frag failed"; exit 1; }
fi
iq() { iqtree3 -s $W/filt/$b.faa -m MFP -mset LG,WAG,JTT,Q.pfam -mfreq FU,F -mrate G,R -cmax 6 --abayes \
         -T $T -seed 20261019 --prefix $W/gt/$b $1 --quiet > $W/log/$b.iqtree 2>&1; }
if [ -s $W/gt/$b.ckp.gz ]; then iq "" || iq -redo; else iq -redo; fi || { echo "$b iqtree failed"; exit 1; }
[ -s $W/gt/$b.treefile ] || { echo "$b no treefile"; exit 1; }
rm -f $W/gt/$b.ckp.gz $W/gt/$b.mldist $W/gt/$b.model.gz $W/gt/$b.uniqueseq.phy
