#!/bin/bash
# R5 D2a re-launch (69c guard fix, DEVIATIONS D-R5-1): simulations at -P 16 under memguard; done jobs skipped by 63a3.
B=~/deepsea_conv; O=$B/work/addon/R5; cd $B
sim1() { nice -n 19 "$@" > /dev/null 2>&1 < /dev/null; return 0; }
export -f sim1
rm -f $O/SIM.DONE; date +"%F %T sim2 start" >> $O/run.log
bash scripts/memguard.sh 48 bash -c "xargs -P 16 -L 1 bash -c 'sim1 \"\$@\"' _ < $O/sim_jobs.txt"; echo "sim2 rc=$?" >> $O/run.log
date +"%F %T sim2 done" >> $O/run.log; touch $O/SIM.DONE
exit 0
