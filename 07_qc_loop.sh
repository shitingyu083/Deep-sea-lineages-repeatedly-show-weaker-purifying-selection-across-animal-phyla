#!/bin/bash
# compleasm (metazoa_odb12) QC daemon: processes genomes as downloads finish.
# Also runs protein-mode QC on the local mollusc proteomes once.
# Usage (from ~/deepsea_conv): nohup bash scripts/07_qc_loop.sh 4 16 > logs/qc.log 2>&1 &
set -u
R=~/deepsea_conv; E=$R/envs/dc/bin; export PATH=$E:$PATH
NJ=${1:-4}; T=${2:-16}
L=$R/db/compleasm; Q=$R/work/qc; mkdir -p $Q

genome_of() {  # flat (05b) or datasets layout (05)
  local d=$R/data/ncbi/$1
  [ -s $d/genomic.fna.gz ] && { echo $d/genomic.fna.gz; return; }
  ls $d/ncbi_dataset/data/$1/*_genomic.fna.gz 2>/dev/null | grep -v cds_from | head -1
}
export -f genome_of; export R E L Q T

qc_one() {
  acc=$1; o=$Q/$acc
  [ -s $o/summary.txt ] && return 0
  [ -f $o/.fail ] && return 0
  g=$(genome_of $acc); [ -z "$g" ] && return 0
  rm -rf $o; mkdir -p $o
  compleasm run -a $g -o $o -l metazoa_odb12 -L $L -t $T > $o/run.log 2>&1 || { echo "QCFAIL $acc" >> $Q/failed.txt; touch $o/.fail; }
  rm -rf $o/metazoa_odb12/hmmer_output 2>/dev/null   # keep it small
}
export -f qc_one

# local proteomes (protein mode), once
for f in $R/data/proteomes/*.faa; do
  c=$(basename $f .faa); o=$Q/local_$c
  [ -s $o/summary.txt ] && continue
  mkdir -p $o
  compleasm protein -p $f -l metazoa_odb12 -L $L -t $((NJ*T)) -o $o > $o/run.log 2>&1 || echo "QCFAIL local_$c" >> $Q/failed.txt
done

while true; do
  fin=0; [ -f $R/data/ncbi/DOWNLOAD.DONE ] && fin=1
  ls $R/data/ncbi | grep ^GC | while read a; do [ -f $R/data/ncbi/$a/.done ] && [ ! -s $Q/$a/summary.txt ] && [ ! -f $Q/$a/.fail ] && echo $a; done > $Q/todo.txt
  if [ -s $Q/todo.txt ]; then
    xargs -P $NJ -I{} bash -c 'qc_one {}' < $Q/todo.txt
  elif [ $fin = 1 ]; then
    break
  else
    sleep 300
  fi
done
echo QCDONE > $Q/QC.DONE
