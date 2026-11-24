#!/bin/bash
# Panel-v3 additions for one clade (2026-10-20). The meta/v3_hold.txt species that passed QC were held out of the
# v1/v2 OrthoFinder runs (16_orthofinder_setup.py, 29_of_add_v2.py) but are in the final species tree (tree_v3).
#   1 wait for the clade's base run (DONE_MARKER: work/orthofinder/<clade>/.done or .v2_done; .fail aborts)
#   2 base = WorkingDirectory of the newest HOG run (21_hog_select.latest_results); stage the QC-pass HOLD members
#     missing from it (29_of_add_v3.py -> add_v3/, list in add_v3.txt)
#   3 move the old selection prep work/codon/<clade> to work/superseded/ (43_autochain.sh can then not start the
#     selection on it), OrthoFinder -b <base Results dir> -f add_v3 -n v3 (only the new searches; the base Log.txt
#     and previous_wd.txt chain the earlier WorkingDirectories; results in <base>/../Results_v3, found by
#     21_hog_select.latest_results). -b must be the Results dir (with Log.txt): given its WorkingDirectory, OrthoFinder
#     2.5.5 falls back to the old-structure locator and sees only that directory's Species*.fa (failed 2026-10-20 16:31)
#   4 27_sel_prep_v2.sh (21 -> 22) on the v3 HOGs -> work/codon/<clade>/PREP.DONE; 43_autochain.sh then launches
#     the final-tree selection (sel_<clade>: BL.DONE + PREP.DONE)
# Markers work/orthofinder/<clade>/.v3_{started,done,fail}; resumable (rerun after a failure).
# Usage: bash scripts/29d_of_add_v3.sh CLADE DONE_MARKER [THREADS=48] [ALIGN_THREADS=8] [NPROC=16]
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH; PY=$R/envs/dc/bin/python
c=$1; DM=$2; T=${3:-48}; A=${4:-8}; NP=${5:-16}
d=work/orthofinder/$c; L=$d/of_add_v3.log
log() { echo "$(date +%F_%T) [29d $c] $*"; }
latest() { $PY -c "import sys; sys.path.insert(0, 'scripts'); from importlib import import_module as m; print(m('21_hog_select').latest_results('$c') or '')"; }
ulimit -n 65536 || { log "cannot raise the nofile limit"; exit 1; }
log "start (pid $$)"
if [ ! -e $d/.v3_done ]; then
  while [ ! -e $DM ]; do [ -e $d/.fail ] && { log "base run failed ($d/.fail)"; exit 1; }; sleep 300; done
  res=$(latest)
  { [ -n "$res" ] && [ -s $res/WorkingDirectory/SpeciesIDs.txt ]; } || { log "no base results"; exit 1; }
  wd=$(readlink -f $res/WorkingDirectory)
  case "$res" in *Results_v3*) log "newest results are already v3 ($res) but no .v3_done: check by hand"; exit 1;; esac
  $PY scripts/29_of_add_v3.py $c $wd > $d/add_v3.txt 2>&1; rc=$?
  cat $d/add_v3.txt
  [ $rc = 3 ] && { touch $d/.v3_done; log "nothing to add"; exit 0; }
  [ $rc = 0 ] || { log "staging failed (rc $rc)"; exit 1; }
  if [ -d work/codon/$c ]; then
    mkdir -p work/superseded; s=work/superseded/codon_${c}_v2prep_$(date +%Y%m%d_%H%M%S)
    mv work/codon/$c $s || { log "cannot move work/codon/$c"; exit 1; }
    log "old prep moved to $s"
  fi
  rd=$(readlink -f $res); rm -f $d/.v3_fail; touch $d/.v3_started $d/.v3_stamp
  log "orthofinder -b $rd -f $d/add_v3 -t $T -a $A -X -n v3 (nofile $(ulimit -n))"
  orthofinder -b $rd -f $d/add_v3 -t $T -a $A -X -n v3 > $L 2>&1; rc=$?
  n0=$(find $(dirname $rd) -maxdepth 3 -path "*/Results_v3*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" -newer $d/.v3_stamp 2>/dev/null | head -1)
  { [ $rc = 0 ] && [ -n "$n0" ]; } || { touch $d/.v3_fail; log "orthofinder failed rc=$rc: $(grep -m1 -i error $L)"; exit 1; }
  touch $d/.v3_done; log "v3 HOGs: $n0"
fi
res=$(latest)
case "$res" in *Results_v3*) ;; *) log "latest_results is not the v3 run ($res): 21_hog_select.py not patched?"; exit 1;; esac
if [ -e work/codon/$c/PREP.DONE ] && [ work/codon/$c/PREP.DONE -nt $d/.v3_done ]; then log "PREP.DONE exists (v3)"; exit 0; fi
if [ -e work/codon/$c/PREP.DONE ]; then log "work/codon/$c/PREP.DONE predates the v3 run: move it away by hand"; exit 1; fi
log "27_sel_prep_v2.sh $NP $c on $res"
bash scripts/27_sel_prep_v2.sh $NP $c & p=$!
sleep 20; echo 0 > /proc/self/autogroup 2>/dev/null   # 27 sets its session to nice 10; the prep is on the critical path
wait $p
[ -e work/codon/$c/PREP.DONE ] || { log "prep failed (work/codon/$c/prep.log)"; exit 1; }
log "PREP.DONE: $(tail -n 1 work/codon/$c/prep.log)"
exit 0
