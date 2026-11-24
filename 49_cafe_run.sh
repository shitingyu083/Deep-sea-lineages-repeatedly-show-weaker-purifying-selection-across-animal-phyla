#!/bin/bash
# CAFE5 per clade (2026-10-25): prep (49_cafe_prep.py) -> one-rate base model with error model (-e) -> deep two-rate
# model (-y lambda_deep.nwk) and the K null shallow two-rate models with that error model -> 49_cafe_collect.py.
# A run that fails with the error model is repeated without it (logged). Each clade runs its fits serially with -c TH.
# Usage: bash 49_cafe_run.sh "CLADES" [TH]   ; markers work/cafe/<C>/CAFE.DONE, work/cafe/CAFE_ALL.DONE
cd ~/deepsea_conv
CL=${1:-"Chondrichthyes Porifera Arthropoda Echinodermata Cnidaria Annelida Mollusca Actinopteri"}
TH=${2:-4}
PY=envs/dc/bin/python; CAFE=envs/dc/bin/cafe5
fit() {  # dir name extra-args...
  local D=$1 N=$2; shift 2
  [ -s $D/$N/Base_results.txt ] && return 0
  rm -rf $D/$N
  nice -n 10 $CAFE -i $D/counts.tsv -t $D/tree.nwk -c $TH -o $D/$N "$@" > $D/logs/$N.log 2>&1
  local rc=$?; echo "$(date +%F_%T) $N rc=$rc" >> $D/logs/runs.log; return $rc
}
one() {
  local C=$1 D=work/cafe/$1
  mkdir -p $D/logs
  $PY scripts/49_cafe_prep.py $C > $D/logs/prep.log 2>&1 || { echo "$(date +%F_%T) prep FAIL" >> $D/logs/runs.log; return; }
  local E=""
  if fit $D base_e -e && [ -s $D/base_e/Base_error_model.txt ]; then
    E="-e$D/base_e/Base_error_model.txt"
  else
    echo "$(date +%F_%T) error model failed -> no error model" >> $D/logs/runs.log
    fit $D base_noe
  fi
  fit $D deep -y $D/lambda_deep.nwk $E
  for y in $D/lambda_null*.nwk; do n=$(basename $y .nwk); fit $D ${n#lambda_} -y $y $E; done
  $PY scripts/49_cafe_collect.py $C > $D/logs/collect.log 2>&1
  touch $D/CAFE.DONE
}
for C in $CL; do one $C & done
wait
touch work/cafe/CAFE_ALL.DONE
exit 0
