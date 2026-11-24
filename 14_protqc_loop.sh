#!/bin/bash
# Proteome QC daemon: compleasm protein (metazoa_odb12) on every final proteome,
# re-run whenever the proteome is newer than its summary; ends after ANNOT.DONE
# once nothing is left.
# 2026-10-19: each job runs on a private lineage view ($Q/<code>/L = symlinks to
# the shared lineage data + a private copy of file_versions.tsv). compleasm
# re-fetches file_versions.tsv.hash into -L at every start, and jobs sharing one
# -L read each other's half-written hash file (IndexError at compleasm.py:124):
# 12 false fails on 2026-10-17/18. Up to 3 tries per job; failed-try logs kept.
# Usage: nohup bash scripts/14_protqc_loop.sh 4 8 > logs/protqc_loop.log 2>&1 &
set -u
R=~/deepsea_conv; export PATH=$R/envs/dc/bin:$PATH
NJ=${1:-4}; T=${2:-8}
L=$R/db/compleasm; Q=$R/work/protqc; mkdir -p $Q
export R L Q T
pq_one() {
  c=$1; f=$R/data/proteomes/$c.faa; o=$Q/$c
  rm -rf $o; mkdir -p $o/L
  for x in metazoa_odb12 eukaryota_odb12 placement_files; do
    [ -e $L/$x ] && ln -s $L/$x $o/L/$x
    [ -e $L/$x.done ] && ln -s $L/$x.done $o/L/$x.done
  done
  cp $L/file_versions.tsv $L/file_versions.tsv.done $o/L/
  ok=0
  for try in 1 2 3; do
    sleep $(( RANDOM % 10 + (try - 1) * 30 ))
    if compleasm protein -p $f -l metazoa_odb12 -L $o/L -t $T -o $o > $o/run.log 2>&1; then
      ok=1; break
    fi
    mv $o/run.log $o/run.try$try.log
    find $o -mindepth 1 -maxdepth 1 ! -name L ! -name 'run.try*.log' -exec rm -rf {} +
  done
  [ $ok = 1 ] || touch $o/.fail
  rm -rf $o/metazoa_odb12_hmmsearch_output $o/L
}
export -f pq_one
while true; do
  done_flag=0; [ -f $R/work/ANNOT.DONE ] && done_flag=1
  for f in $R/data/proteomes/*.faa; do
    c=$(basename $f .faa); s=$Q/$c/summary.txt
    [ -n "$(find $f -mmin -3)" ] && continue   # may still be being written
    [ -f $Q/$c/.fail ] && [ $Q/$c/.fail -nt $f ] && continue
    { [ ! -s $s ] || [ $f -nt $s ]; } && echo $c
  done > $Q/todo.txt
  if [ -s $Q/todo.txt ]; then
    xargs -P $NJ -I{} bash -c 'pq_one {}' < $Q/todo.txt
  elif [ $done_flag = 1 ]; then
    break
  else
    sleep 300
  fi
done
echo DONE > $Q/PROTQC.DONE
