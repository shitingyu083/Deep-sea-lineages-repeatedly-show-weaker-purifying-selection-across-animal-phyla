#!/bin/bash
# CAFE5 rescue for clades whose base fit failed to initialise (2026-10-28). Filters high-range families
# (49b_cafe_filter.py tier 1, then tier 2), same model sequence as 49_cafe_run.sh, two-rate fits run in parallel.
# Output work/cafe/<C>_rescue/, collected by the unchanged 49_cafe_collect.py; marker work/cafe/CAFE_RESCUE.DONE
cd ~/deepsea_conv
PY=envs/dc/bin/python; CAFE=envs/dc/bin/cafe5; TH=4
fit() { local D=$1 N=$2; shift 2; rm -rf $D/$N
  nice -n 10 $CAFE -i $D/counts.tsv -t $D/tree.nwk -c $TH -o $D/$N "$@" > $D/logs/$N.log 2>&1
  echo "$(date +%F_%T) $N rc=$?" >> $D/logs/runs.log; [ -s $D/$N/Base_results.txt ]; }
one() {
  local C=$1 D=work/cafe/$1_rescue E="" ok=0
  for tier in 1 2; do
    rm -rf $D; $PY scripts/49b_cafe_filter.py work/cafe/$C $D $tier > /tmp/cafe_filter_${C}_sty.log 2>&1
    mkdir -p $D/logs; cat /tmp/cafe_filter_${C}_sty.log >> $D/logs/runs.log
    if fit $D base_e -e && [ -s $D/base_e/Base_error_model.txt ]; then E="-e$D/base_e/Base_error_model.txt"; ok=1; break; fi
    if fit $D base_noe; then ok=1; break; fi
  done
  [ $ok = 1 ] || { echo "$(date +%F_%T) rescue FAIL" >> $D/logs/runs.log; touch $D/CAFE.FAIL; return; }
  fit $D deep -y $D/lambda_deep.nwk $E &
  for y in $D/lambda_null*.nwk; do n=$(basename $y .nwk); fit $D ${n#lambda_} -y $y $E & done
  wait
  $PY scripts/49_cafe_collect.py ${C}_rescue > $D/logs/collect.log 2>&1
  touch $D/CAFE.DONE
}
for C in Mollusca Actinopteri_m; do one $C & done
wait
touch work/cafe/CAFE_RESCUE.DONE
exit 0
