#!/bin/bash
# v2 selection prep, the tree-independent part: 21_hog_select.py (near-single-copy HOGs from the
# v2 per-clade OrthoFinder run) + 22_codon_aln.py (codon alignments) per clade, clades in parallel.
# A clade is skipped unless its newest N0.tsv dates from the v2 runs (>= 2026-10-19); Cnidaria
# first waits (max 72 h) for its -fg rerun (work/orthofinder/Cnidaria/.v2_done or .v2_fail, 29c_of_fg_retry.sh).
# 23_label_trees.py + 24_selection_run.sh need the final species tree and run later.
# Usage: bash scripts/27_sel_prep_v2.sh NPROC CLADE...   (detached; logs in work/codon/<clade>/)
R=~/deepsea_conv; cd $R; PY=$R/envs/dc/bin/python
NP=${1:-16}; shift
echo 10 > /proc/self/autogroup 2>/dev/null   # yield CPU to the gene-tree workers
T0=$(date -d '2026-10-19 00:00' +%s)
prep() {
  c=$1; W=work/codon/$c; mkdir -p $W
  [ -e $W/PREP.DONE ] && return 0
  if [ "$c" = Cnidaria ]; then
    n=0; while [ ! -e work/orthofinder/Cnidaria/.v2_done ] && [ ! -e work/orthofinder/Cnidaria/.v2_fail ] && [ $n -lt 864 ]; do
      sleep 300; n=$((n+1)); done
  fi
  read r mt <<< "$($PY -c "import os,sys; sys.path.insert(0,'scripts'); from importlib import import_module as m; r=m('21_hog_select').latest_results('$c'); print(r, int(os.path.getmtime(r+'/Phylogenetic_Hierarchical_Orthogroups/N0.tsv')) if r else 0)")"
  echo "$(date +%F_%T) $c results=$r N0_mtime=$(date -d @${mt:-0} +%F_%T)" >> $W/prep.log
  if [ "${mt:-0}" -lt "$T0" ]; then echo "$(date +%F_%T) $c newest N0 predates v2: skipped" >> $W/prep.log; return 1; fi
  $PY scripts/21_hog_select.py $c > $W/prep21.log 2>&1 || { echo "$(date +%F_%T) $c 21 failed" >> $W/prep.log; return 1; }
  echo "$(date +%F_%T) $c hogs=$(($(wc -l < $W/hogs.tsv)-1))" >> $W/prep.log
  $PY scripts/22_codon_aln.py $c $NP > $W/prep22.log 2>&1 || { echo "$(date +%F_%T) $c 22 failed" >> $W/prep.log; return 1; }
  echo "$(date +%F_%T) $c aln=$(ls $W/aln | wc -l) $(cut -f2 $W/aln_report.tsv | sort | uniq -c | tr -s ' ' | tr '\n' ' ')" >> $W/prep.log
  touch $W/PREP.DONE
}
for c in "$@"; do prep $c & done
wait
echo "$(date +%F_%T) all done: $*"
