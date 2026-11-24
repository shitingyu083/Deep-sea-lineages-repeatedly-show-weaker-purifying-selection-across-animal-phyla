#!/bin/bash
# One locus of the fast DRAFT species tree (19f_fast_draft.sh): alignment = the work/tree_v2 MAFFT L-INS-i
# alignment when that run had finished it, else MAFFT --auto -> trimAl -automated1 -> fragment filter (as in
# 19c) -> FastTree LG+Gamma (SH-like local supports). Every step resumable. Runs under flock on
# W/claim/<locus>.lock so a second worker session can share the list; a held locus is skipped.
# The draft is only a scaffold for the downstream pipelines; everything is re-run on the final 19c tree.
W=$1; b=$2; P=~/anaconda3/envs/py310env/bin/python; S=~/deepsea_conv/scripts
SRC=~/deepsea_conv/work/tree_v2
[ -s $W/gt/$b.tre ] && exit 0
[ -e $W/skip/$b ] && exit 0
mkdir -p $W/claim
exec 9> $W/claim/$b.lock
flock -n 9 || exit 0
[ -s $W/gt/$b.tre ] && exit 0
[ -e $W/skip/$b ] && exit 0
if [ ! -s $W/aln/$b.faa ]; then
  if [ -s $SRC/aln/$b.faa ]; then
    cp $SRC/aln/$b.faa $W/aln/$b.tmp && mv $W/aln/$b.tmp $W/aln/$b.faa && echo linsi > $W/aln/$b.method \
      || { echo "$b copy failed"; exit 1; }
  else
    mafft --auto --thread 1 --quiet $SRC/loci/$b.faa > $W/aln/$b.tmp 2> $W/log/$b.mafft \
      && [ -s $W/aln/$b.tmp ] && mv $W/aln/$b.tmp $W/aln/$b.faa && echo auto > $W/aln/$b.method \
      || { echo "$b mafft failed"; exit 1; }
  fi
fi
if [ ! -s $W/filt/$b.faa ]; then
  trimal -in $W/aln/$b.faa -out $W/trim/$b.faa -automated1 > $W/log/$b.trimal 2>&1 || { echo "$b trimal failed"; exit 1; }
  $P $S/19c_tools.py frag $W/trim/$b.faa $W/filt/$b.faa 0.5 100 4; rc=$?
  if [ $rc = 2 ]; then echo "too short/few taxa after trimming" > $W/skip/$b; exit 0; fi
  [ $rc = 0 ] || { echo "$b frag failed"; exit 1; }
fi
OMP_NUM_THREADS=1 FastTree -lg -gamma -quiet -nopr $W/filt/$b.faa > $W/gt/$b.tmp 2> $W/log/$b.fasttree \
  && [ -s $W/gt/$b.tmp ] && mv $W/gt/$b.tmp $W/gt/$b.tre || { echo "$b fasttree failed"; exit 1; }
