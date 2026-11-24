#!/bin/bash
# One locus of the FINAL species tree (19h_final_tree.sh, 2026-10-19): MAFFT L-INS-i -> trimAl -automated1 ->
# fragment filter -> IQ-TREE (ModelFinder over LG/WAG/JTT/Q.pfam x FU/F x G4, default tree search, aBayes);
# every step resumable. As 19c_locus.sh, plus: PATH exported (the workers start from a bare shell); a memory gate
# before each heavy step (waits until MemAvailable >= MINMEM GB, default 80, never fails); wall time and peak
# memory of MAFFT and IQ-TREE in W/log/<locus>.{mafft,iqtree}.time (/usr/bin/time "%e %M"); the aligner and its
# options in W/aln/<locus>.method; MAFFT threads = argument 3 (W/cost.tsv); IQ-TREE threads from the trimmed
# alignment, ceil(h / 3) with h = 1.5 * n * L / 60040 the expected single-thread hours (the default search on a
# 316 x 190 locus took ~1.5 h in the timing test), 1..8 (W/log/<locus>.iqtree.threads).
# Each locus runs under flock on W/claim/<locus>.lock (released when the whole process tree ends, also after a
# crash); a held locus is skipped. W/CUT stops new loci from starting; loci already running finish. A locus
# counts as finished when its IQ-TREE log has the closing "Date and Time:" line, not when its treefile exists
# (IQ-TREE writes intermediate treefiles); an unfinished run resumes from its checkpoint.
# Usage: bash scripts/19h_locus.sh W LOCUS [MAFFT_THREADS=1]
export PATH=~/deepsea_conv/envs/dc/bin:$PATH
W=$1; b=$2; T=${3:-1}; P=~/anaconda3/envs/py310env/bin/python; S=~/deepsea_conv/scripts
# IQ-TREE writes intermediate treefiles: a locus is finished once its log ends the run (2026-10-20 fix)
done_() { [ -s $W/gt/$b.treefile ] && grep -q "^Date and Time:" $W/gt/$b.log 2>/dev/null; }
done_ && exit 0
[ -e $W/skip/$b ] && exit 0
[ -e $W/CUT ] && exit 0
mkdir -p $W/claim
exec 9> $W/claim/$b.lock
flock -n 9 || exit 0
done_ && exit 0   # finished by another run between the first check and the lock
[ -e $W/skip/$b ] && exit 0
memwait() { while [ $(awk '/^MemAvailable/{print int($2/1048576)}' /proc/meminfo) -lt ${MINMEM:-80} ]; do sleep 60; done; }
if [ ! -s $W/aln/$b.faa ]; then
  memwait
  /usr/bin/time -f "%e %M" -o $W/log/$b.mafft.time mafft --localpair --maxiterate 1000 --thread $T --quiet \
      $W/loci/$b.faa > $W/aln/$b.tmp 2> $W/log/$b.mafft \
    && [ -s $W/aln/$b.tmp ] && mv $W/aln/$b.tmp $W/aln/$b.faa || { echo "$b mafft failed"; exit 1; }
  echo "MAFFT $(mafft --version 2>&1 | head -1 | awk '{print $1}') L-INS-i (--localpair --maxiterate 1000), $T threads" \
    > $W/aln/$b.method
fi
if [ ! -s $W/filt/$b.faa ]; then
  trimal -in $W/aln/$b.faa -out $W/trim/$b.faa -automated1 > $W/log/$b.trimal 2>&1 || { echo "$b trimal failed"; exit 1; }
  $P $S/19c_tools.py frag $W/trim/$b.faa $W/filt/$b.faa 0.5 100 4; rc=$?
  if [ $rc = 2 ]; then echo "too short/few taxa after trimming" > $W/skip/$b; exit 0; fi
  [ $rc = 0 ] || { echo "$b frag failed"; exit 1; }
fi
TIQ=$(awk '/^>/{next} {c += length($0)} END{h = 1.5 * c / 60040; t = int(h / 3); if (t < h / 3) t++
           if (t < 1) t = 1; if (t > 8) t = 8; print t}' $W/filt/$b.faa)
iq() { /usr/bin/time -f "%e %M" -o $W/log/$b.iqtree.time iqtree3 -s $W/filt/$b.faa -m MFP -mset LG,WAG,JTT,Q.pfam \
         -mfreq FU,F -mrate G --abayes -T $TIQ -seed 20261019 --prefix $W/gt/$b $1 --quiet > $W/log/$b.iqtree 2>&1; }
memwait
echo $TIQ > $W/log/$b.iqtree.threads
if [ -s $W/gt/$b.ckp.gz ]; then iq "" || iq -redo; else iq -redo; fi || { echo "$b iqtree failed"; exit 1; }
done_ || { echo "$b IQ-TREE did not finish"; exit 1; }
rm -f $W/gt/$b.ckp.gz $W/gt/$b.mldist $W/gt/$b.model.gz $W/gt/$b.uniqueseq.phy
