#!/bin/bash
# R2 main power run: wait for pilot2, install the Actinopteri_m path fix of 63a2, gate check, calib, jobs.
B=$HOME/deepsea_conv; PY=$B/envs/dc/bin/python; R=$B/work/addon/R2/main
cd $B || exit 0
while [ ! -e work/addon/R2/pilot2/PILOT2.DONE ]; do sleep 120; done
if [ -e /tmp/addon_sty/63a2_powersim_sty.v2b ]; then
  cp scripts/63a2_powersim_sty.py backups/63a2_powersim_sty.py.$(date +%Y%m%d_%H%M%S)
  cp /tmp/addon_sty/63a2_powersim_sty.v2b scripts/63a2_powersim_sty.py
fi
mkdir -p $R
$PY scripts/63b_sample_sty.py < /dev/null > $R/sample.log 2>&1
if [ $? -ne 0 ]; then echo "gate fail or sampler error" >> $R/sample.log; touch $R/GATE_FAIL; exit 0; fi
cat $R/calib_jobs.txt | xargs -P 40 -I{} bash -c "{} < /dev/null" > $R/calib.log 2>&1
touch $R/CALIB.DONE
cat $R/jobs.txt | xargs -P 48 -I{} bash -c "{} < /dev/null" > $R/run.log 2>&1
find $R -name res.tsv | sort | xargs cat | awk 'NR==1||$1!="clade"' > $R/res_all.tsv
touch $R/R2MAIN.DONE
exit 0
