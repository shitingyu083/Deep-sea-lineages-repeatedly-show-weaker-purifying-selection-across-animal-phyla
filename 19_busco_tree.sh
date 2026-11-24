#!/bin/bash
# Per-locus MAFFT -> trimAl -> IQ-TREE gene trees, then ASTRAL species tree and
# a partitioned concatenated IQ-TREE ML tree (UFBoot 1000).
# Usage: nohup bash scripts/19_busco_tree.sh NJ > logs/tree.log 2>&1 &
set -u
R=~/deepsea_conv; export PATH=$R/envs/dc/bin:$PATH
NJ=${1:-32}; W=${TREE_DIR:-$R/work/tree}; mkdir -p $W/aln $W/trim $W/gt
export W
one() {
  b=$1
  [ -s $W/gt/$b.treefile ] && return 0
  mafft --auto --thread 2 --quiet $W/loci/$b.faa > $W/aln/$b.faa
  trimal -in $W/aln/$b.faa -out $W/trim/$b.faa -automated1 > /dev/null 2>&1
  n=$(grep -c ">" $W/trim/$b.faa); [ "$n" -lt 4 ] && return 0
  iqtree -s $W/trim/$b.faa -m MFP -mset LG,WAG,JTT -T 2 -B 1000 --prefix $W/gt/$b -quiet > /dev/null 2>&1
}
export -f one
ls $W/loci | sed 's/.faa$//' | xargs -P $NJ -I{} bash -c 'one {}'
# ASTRAL species tree (topology + local posterior support)
cat $W/gt/*.treefile > $W/genetrees.tre
astral -i $W/genetrees.tre -o $W/astral.tre > $W/astral.log 2>&1
# branch lengths: concatenated alignment on the fixed ASTRAL topology, each locus with
# the model ModelFinder chose for its gene tree
{ echo "#nexus"; echo "begin sets;"
  for t in $W/gt/*.iqtree; do
    b=$(basename $t .iqtree); m=$(grep -m1 "^Best-fit model" $t | awk '{print $NF}')
    echo "  charset $b = $W/trim/$b.faa: *;"
  done
  printf "  charpartition mine = "
  for t in $W/gt/*.iqtree; do
    b=$(basename $t .iqtree); m=$(grep -m1 "^Best-fit model" $t | awk '{print $NF}')
    printf "%s:%s, " "$m" "$b"
  done | sed 's/, $/;\n/'
  echo "end;"; } > $W/partition.nex
iqtree -p $W/partition.nex -te $W/astral.tre -T 64 --prefix $W/concat_bl > $W/concat_bl.stdout 2>&1
echo DONE > $W/TREE.DONE
