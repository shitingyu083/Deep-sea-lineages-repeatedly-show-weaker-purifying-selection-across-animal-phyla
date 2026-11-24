#!/bin/bash
# Gene trees for RERconverge: per-locus branch lengths on the fixed species topology. IQ-TREE 3 with the
# locus model of locus_stats.tsv (LG+G4), topology fixed (-te), one thread per locus, GNU parallel.
# Inputs from 36_rer_prep.py; output OUT/rer_trees.txt ("locus<TAB>newick", the readTrees format).  TREES.DONE
# Usage: bash scripts/36_rer_trees.sh SRC TOPO OUT [NJ=48]
#   SRC = species-tree work dir (final_aln/, locus_stats.tsv); TOPO = rooted species topology
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
SRC=$(readlink -f $1); TOPO=$(readlink -f $2); mkdir -p $3/log $3/bl; OUT=$(readlink -f $3); NJ=${4:-48}
log() { echo "$(date '+%F %T') [36] $*"; }
[ -e $OUT/TREES.DONE ] && exit 0
if [ ! -s $OUT/jobs.tsv ]; then
  python scripts/36_rer_prep.py $SRC $TOPO $OUT > $OUT/log/prep.txt 2>&1 || { log "prep failed"; cat $OUT/log/prep.txt; exit 1; }
  printf "SRC\t%s\nTOPO\t%s\nTOPO_md5\t%s\n" $SRC $TOPO $(md5sum < $TOPO | cut -c1-32) > $OUT/inputs.tsv
fi
log "$(tail -1 $OUT/log/prep.txt); branch lengths with $NJ parallel IQ-TREE jobs"
tail -n +2 $OUT/jobs.tsv | parallel -j $NJ --memfree 40G --colsep '\t' --joblog $OUT/log/parallel.joblog \
  "[ -s $OUT/bl/{1}.treefile ] || iqtree3 -s $OUT/aln/{1}.faa -te $OUT/topo/{1}.tre -m {3} -T 1 --keep-ident \
     -seed 20261019 --prefix $OUT/bl/{1} -redo -quiet > /dev/null 2>&1; rm -f $OUT/bl/{1}.ckp.gz"
n=0; miss=""
tail -n +2 $OUT/jobs.tsv | cut -f1 > $OUT/loci.tmp
while read l; do
  if [ -s $OUT/bl/$l.treefile ]; then printf "%s\t%s\n" $l "$(head -1 $OUT/bl/$l.treefile)"; n=$((n + 1))
  else miss="$miss $l"; fi
done < $OUT/loci.tmp > $OUT/rer_trees.tmp
rm -f $OUT/loci.tmp
[ -n "$miss" ] && { echo "$miss" | tr ' ' '\n' | sed '/^$/d' > $OUT/failed_loci.txt; log "no tree for$miss"; }
nj=$(($(wc -l < $OUT/jobs.tsv) - 1))
[ $n -ge $((nj * 95 / 100)) ] || { log "only $n of $nj loci have trees"; exit 1; }
mv $OUT/rer_trees.tmp $OUT/rer_trees.txt
echo "IQ-TREE $(iqtree3 --version 2>&1 | head -1 | awk '{for (i = 1; i < NF; i++) if ($i == "version") {print $(i + 1); exit}}') -te (fixed species topology) -m locus model, --keep-ident, seed 20261019" > $OUT/method.txt
touch $OUT/TREES.DONE
log "TREES.DONE: $n of $nj loci"
