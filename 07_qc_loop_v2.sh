#!/bin/bash
# compleasm genome-mode QC daemon, v2.
# Changes vs v1: queue only accessions whose genome is on disk (05d phase 1 marks
# annotated accessions .done before the genome arrives); skip accessions with a
# compleasm already running (inherited from v1); finish only after DOWNLOAD.DONE
# and no .nogenome markers remain.
# Usage: nohup bash scripts/07_qc_loop_v2.sh 6 16 > logs/qc.log 2>&1 &
set -u
R=~/deepsea_conv; E=$R/envs/dc/bin; export PATH=$E:$PATH
NJ=${1:-4}; T=${2:-16}
L=$R/db/compleasm; Q=$R/work/qc; mkdir -p $Q

genome_of() {
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
  rm -rf $o/metazoa_odb12/hmmer_output 2>/dev/null
}
export -f qc_one

while true; do
  fin=0; [ -f $R/data/ncbi/DOWNLOAD.DONE ] && [ -f $R/data/ncbi/DOWNLOAD_B.DONE ] && ! ls $R/data/ncbi/*/.nogenome >/dev/null 2>&1 && fin=1
  running=$(pgrep -af "compleasm run" | grep -o "$Q/GC[^ ]*" | xargs -n1 basename 2>/dev/null)
  ls $R/data/ncbi | grep ^GC | while read a; do
    [ -f $R/data/ncbi/$a/.done ] || continue
    [ -s $Q/$a/summary.txt ] || [ -f $Q/$a/.fail ] && continue
    echo "$running" | grep -qx "$a" && continue
    [ -n "$(genome_of $a)" ] && echo $a
  done > $Q/todo.txt
  if [ -s $Q/todo.txt ]; then
    xargs -P $NJ -I{} bash -c 'qc_one {}' < $Q/todo.txt
  elif [ $fin = 1 ] && [ -z "$running" ]; then
    break
  else
    sleep 300
  fi
done
echo QCDONE > $Q/QC.DONE
