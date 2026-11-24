#!/bin/bash
# R5 single retry (DEVIATIONS D-R5-2): re-run once every RELAX job without a parseable json and every D2a sim without res.tsv.
# First-attempt json/log kept as *.try1. Marker: R5/RETRY.DONE.
B=~/deepsea_conv; O=$B/work/addon/R5; E=$B/envs/dc/bin; cd $B
$E/python -c "
import json,os
fo=open('$O/retry_relax.txt','w')
for l in open('$O/relax_jobs.txt'):
    o=l.split()[0]
    try: json.load(open(o))['test results']['p-value']
    except Exception: fo.write(l)
fs=open('$O/retry_sim.txt','w')
for l in open('$O/sim_jobs.txt'):
    if not os.path.exists(l.split()[-1]+'/res.tsv'): fs.write(l)
"
while read o a t g; do [ -e "$o" ] && mv "$o" "$o.try1"; [ -e "$g" ] && mv "$g" "$g.try1"; done < $O/retry_relax.txt
while read -a w; do d=${w[-1]}; for f in o.json hy.log errors.log; do [ -e $d/$f ] && mv $d/$f $d/$f.try1; done; done < $O/retry_sim.txt
date +"%F %T retry start relax=$(wc -l < $O/retry_relax.txt) sim=$(wc -l < $O/retry_sim.txt)" >> $O/run.log
relax1() { [ -s "$1" ] && return 0; nice -n 19 $E/hyphy CPU=1 relax --alignment "$2" --tree "$3" --test FG --reference BG \
  --models Minimal --output "$1" > "$4" 2>&1 < /dev/null; return 0; }
sim1() { nice -n 19 "$@" > /dev/null 2>&1 < /dev/null; return 0; }
export -f relax1 sim1; export E
cat $O/retry_sim.txt | sed 's/^/S /' > $O/retry_all.txt; sed 's/^/R /' $O/retry_relax.txt >> $O/retry_all.txt
one() { k=$1; shift; if [ $k = R ]; then relax1 "$@"; else sim1 "$@"; fi; }
export -f one
bash scripts/memguard.sh 48 bash -c "xargs -P 48 -L 1 bash -c 'one \"\$@\"' _ < $O/retry_all.txt"; echo "retry rc=$?" >> $O/run.log
date +"%F %T retry done" >> $O/run.log; touch $O/RETRY.DONE
exit 0
