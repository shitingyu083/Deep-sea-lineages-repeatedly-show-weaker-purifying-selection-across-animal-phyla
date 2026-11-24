#!/bin/bash
# D-R2-2: after the R2 main run, redo with 63a3 (AliSim --skip-bl-check) every HOG whose calibration failed and every
# job without res.tsv (AliSim rejected mean branch length > 10 for saturated high-omega / episodic classes).
B=$HOME/deepsea_conv; PY=$B/envs/dc/bin/python; R=$B/work/addon/R2/main; S=$B/scripts/63a3_powersim_sty.py
cd $B || exit 0
while [ ! -e $R/R2MAIN.DONE ]; do sleep 120; done
tail -n +2 $R/hogs.tsv | while read C h; do
  [ -s $R/$C/$h/_fit/scale.txt ] || { rm -rf $R/$C/$h/_fit/calwd; echo "$PY $S calib $C $h $R/$C/$h/_fit"; }
done > $R/calib_fix_jobs.txt
cat $R/calib_fix_jobs.txt | xargs -P 48 -I{} bash -c "{} < /dev/null" > $R/calib_fix.log 2>&1
: > $R/jobs_fix.txt
while read l; do od=${l##* }; [ -s $od/res.tsv ] || echo "${l/63a2_powersim_sty.py/63a3_powersim_sty.py}" >> $R/jobs_fix.txt; done < $R/jobs.txt
cat $R/jobs_fix.txt | xargs -P 48 -I{} bash -c "{} < /dev/null" > $R/run_fix.log 2>&1
find $R -name res.tsv | sort | xargs cat | awk 'NR==1||$1!="clade"' > $R/res_all.tsv
echo "calib_fix $(wc -l < $R/calib_fix_jobs.txt) jobs_fix $(wc -l < $R/jobs_fix.txt) res $(($(wc -l < $R/res_all.tsv)-1))" > $R/fix_summary.txt
$PY scripts/63c_r2collect_sty.py < /dev/null > $R/collect.log 2>&1
touch $R/R2FIX.DONE
exit 0
