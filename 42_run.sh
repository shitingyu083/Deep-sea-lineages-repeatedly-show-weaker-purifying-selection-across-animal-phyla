#!/bin/bash
# 42 habitat temperature (2026-10-20): extract (duckdb, py310env) -> temp (envs/dc, h5py).
# Own session, autogroup nice 10; resumable (EXTRACT.DONE); log logs/habitat42.log.
# Usage: cd ~/deepsea_conv; (setsid nohup bash scripts/42_run.sh > /dev/null 2>&1 < /dev/null &)
R=~/deepsea_conv; cd $R
echo 10 > /proc/self/autogroup 2>/dev/null
H=work/habitat; L=logs/habitat42.log; mkdir -p $H
log() { echo "$(date +%F_%T) $*" >> $L; }
echo $$ > $H/run.pid
if [ -e $H/EXTRACT.DONE ]; then
  log "extract already done"
else
  log "extract start (pid $$)"
  ~/anaconda3/envs/py310env/bin/python scripts/42_habitat_temp.py extract >> $L 2>&1
  rc=$?; [ $rc = 0 ] || { log "extract failed (rc $rc)"; exit 1; }
  touch $H/EXTRACT.DONE; rm -rf $H/duck_tmp
  log "extract done: $(($(wc -l < $H/obis_panel_records.tsv) - 1)) rows"
fi
log "temp start"
envs/dc/bin/python scripts/42_habitat_temp.py temp >> $L 2>&1
rc=$?; [ $rc = 0 ] || { log "temp failed (rc $rc)"; exit 1; }
touch $H/TEMP.DONE; log "temp done"
