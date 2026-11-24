#!/bin/bash
# Restart the fish sub77 OrthoFinder -b pass with OrthoFinder's worker queue timeout raised.
# 2026-10-21 05:30: the 30_fish_sub_of.sh -b pass (started 04:05) ran "Initial processing of each
# species" with ONE live worker of -a 8 (~13 min/species, ~17 h for 77): OrthoFinder 2.5.5 workers
# return after cmd_queue.get(True, 1) times out, and at load ~660/256 cores the queue feeder
# (large pickled seqsInfo/Lengths per item) is slower than 1 s, so 7 workers quit at once.
# Fix: get(True, 1) -> get(True, 30) in the env's scripts_of (backup kept; only delays idle-worker
# exit by <=30 s per phase), move the partial Results_sub77 aside, rerun the identical -b command.
# Usage: bash scripts/30c_fish_of_b_restart.sh WRAPPER_PID OF_MAIN_PID
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH
WP=$1; OP=$2; TS=$(date +%Y%m%d_%H%M%S)
W=$R/work/orthofinder/Actinopteri/of/Results_Oct19/WorkingDirectory
D=$R/work/orthofinder/Actinopteri/sub_v2
own() { [ "$(stat -c %u /proc/$1 2>/dev/null)" = 1024 ] && tr '\0' ' ' 2>/dev/null < /proc/$1/cmdline | grep -q "$2"; }
own $WP 30_fish_sub_of.sh || { echo "wrapper $WP not verified"; exit 1; }
own $OP "orthofinder -b $W" || { echo "OF $OP not verified"; exit 1; }
KIDS=$(ps --ppid $OP -o pid= | tr -d ' ')
kill -TERM $WP; sleep 2
for p in $OP $KIDS; do own $p "orthofinder -b" && kill -TERM $p; done; sleep 10
for p in $OP $KIDS; do own $p "orthofinder -b" && kill -KILL $p; done; sleep 2
for p in $WP $OP $KIDS; do [ -e /proc/$p ] && echo "still alive $p"; done
mkdir -p $R/backups/of_src_$TS
cp -p envs/dc/bin/scripts_of/{__main__,orthologues,parallel_task_manager}.py $R/backups/of_src_$TS/
sed -i 's/\.get(True, 1)/.get(True, 30)/' envs/dc/bin/scripts_of/{__main__,orthologues,parallel_task_manager}.py
echo "patched: $(grep -c 'get(True, 30)' envs/dc/bin/scripts_of/{__main__,orthologues,parallel_task_manager}.py | tr '\n' ' ')"
mkdir -p $W/OrthoFinder/_superseded
[ -d $W/OrthoFinder/Results_sub77 ] && mv $W/OrthoFinder/Results_sub77 $W/OrthoFinder/_superseded/Results_sub77_serial_$TS
mv $D/of_b.log $D/of_b.serial_$TS.log
rm -f $D/SUB.FAIL
cat > /tmp/fish_ofb_run_sty.sh <<EOS
cd $R; export PATH=$R/envs/dc/bin:\$PATH; ulimit -n 65536
echo "\$(date +%F_%T) 30c restart: orthofinder -b (queue timeout 30 s)" >> $D/run.log
if orthofinder -b $W -t 48 -a 8 -X -n sub77 > $D/of_b.log 2>&1; then
  touch $D/SUB.DONE; echo "\$(date +%F_%T) orthofinder -b done" >> $D/run.log
else
  touch $D/SUB.FAIL; echo "\$(date +%F_%T) orthofinder -b FAILED" >> $D/run.log
fi
exit 0
EOS
(setsid nohup bash /tmp/fish_ofb_run_sty.sh > $R/logs/fish_ofb_restart.log 2>&1 < /dev/null 7>&- 8>&- &)
sleep 60; tail -n 3 $D/of_b.log
exit 0
