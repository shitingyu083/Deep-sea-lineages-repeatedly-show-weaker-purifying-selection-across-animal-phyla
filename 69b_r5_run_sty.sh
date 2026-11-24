#!/bin/bash
# R5 runner (PREREG_r5.md): D1+D3 RELAX queue at -P 48, D2a simulations at -P 16, both under memguard (48 GB floor).
# Jobs whose output exists are skipped. Markers: R5/RELAX.DONE, R5/SIM.DONE, R5/R5.DONE.
B=~/deepsea_conv; O=$B/work/addon/R5; E=$B/envs/dc/bin; cd $B
relax1() { [ -s "$1" ] && return 0; nice -n 19 $E/hyphy CPU=1 relax --alignment "$2" --tree "$3" --test FG --reference BG \
  --models Minimal --output "$1" > "$4" 2>&1 < /dev/null; return 0; }
export -f relax1; export E
sim1() { nice -n 19 "$@" > /dev/null 2>&1 < /dev/null; return 0; }
export -f sim1
date +"%F %T start" > $O/run.log
( bash scripts/memguard.sh 48 bash -c "xargs -P 16 -L 1 bash -c 'sim1 \"\$@\"' _ < $O/sim_jobs.txt"; echo "sim rc=$?" >> $O/run.log
  date +"%F %T sim done" >> $O/run.log; touch $O/SIM.DONE ) &
bash scripts/memguard.sh 48 bash -c "xargs -P 48 -L 1 bash -c 'relax1 \"\$@\"' _ < $O/relax_jobs.txt"; echo "relax rc=$?" >> $O/run.log
date +"%F %T relax done" >> $O/run.log; touch $O/RELAX.DONE
wait; touch $O/R5.DONE; date +"%F %T all done" >> $O/run.log
exit 0
