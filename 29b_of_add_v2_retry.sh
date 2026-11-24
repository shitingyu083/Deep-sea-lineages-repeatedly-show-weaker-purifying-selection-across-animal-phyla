#!/bin/bash
# Re-run a clade's v2 OrthoFinder -b pass after a failure (29_of_add_v2_loop.sh never retries a clade
# once .v2_started exists). Adds a watchdog for hung `diamond makedb`: on 2026-10-18 Annelida's
# makedb deadlocked in futex wait (0.09 s CPU over 33 h). makedb on a proteome takes seconds, so any
# makedb of this clade alive > 15 min is killed; OrthoFinder then exits non-zero and we retry (<= 3x).
# Usage: bash scripts/29b_of_add_v2_retry.sh CLADE [THREADS=48] [ALIGN_THREADS=8]
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
c=$1; T=${2:-48}; A=${3:-8}
d=work/orthofinder/$c
wd=$(ls -d $d/of/Results_*/WorkingDirectory 2>/dev/null | while read w; do
       ls $w/Blast0_0.txt* >/dev/null 2>&1 && echo $w; done | head -1)
[ -n "$wd" ] || { echo "$c: no v1 WorkingDirectory with BLASTs"; exit 1; }
touch $d/.v2_started
for try in 1 2 3; do
  rm -f $d/.v2_fail $d/.v2_done
  touch $d/.v2_try_stamp
  echo "$(date '+%F %T') $c try $try: orthofinder -b $wd -f $d/add_v2 -t $T -a $A -X"
  nice -n 10 orthofinder -b $wd -f $d/add_v2 -t $T -a $A -X > $d/of_add_v2.log 2>&1 &
  of=$!
  while kill -0 $of 2>/dev/null; do
    sleep 120
    for p in $(pgrep -f "^diamond makedb .*/work/orthofinder/$c/"); do
      et=$(ps -o etimes= -p $p 2>/dev/null | tr -d ' ')
      if [ -n "$et" ] && [ "$et" -gt 900 ]; then
        echo "$(date '+%F %T') watchdog: makedb $p alive ${et}s ($(awk '{print $14+$15}' /proc/$p/stat) ticks) -> kill"
        kill $p; sleep 5; kill -0 $p 2>/dev/null && kill -9 $p
      fi
    done
  done
  wait $of; rc=$?
  cp $d/of_add_v2.log $d/of_add_v2.try$try.log
  n0=$(find $wd/OrthoFinder -path "*Phylogenetic_Hierarchical_Orthogroups/N0.tsv" -newer $d/.v2_try_stamp 2>/dev/null | head -1)
  if [ $rc = 0 ] && [ -n "$n0" ]; then
    touch $d/.v2_done; echo "$(date '+%F %T') $c v2 done"; exit 0
  fi
  echo "$(date '+%F %T') $c try $try failed rc=$rc: $(grep -m1 -i error $d/of_add_v2.log)"
  sleep 60
done
touch $d/.v2_fail; exit 1
