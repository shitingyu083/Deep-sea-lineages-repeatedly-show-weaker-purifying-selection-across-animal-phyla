#!/bin/bash
# A1 re-annotation control driver (2026-10-25; PREREG_addon.md A1). 59_route_reannot.py --select (unless
# OUT/SELECT.DONE exists), --one per selected species (NJ in parallel, T threads each, memguard 48 GB, nice 10; a code
# with REANNOT.DONE or REANNOT.FAIL is skipped), one serial retry of memguard kills and crashes, --collect, then
# 59b_reannot_refit.R -> OUT/A1R.DONE (OUT/A1R.FAIL names the failed step). --one exit codes: 0 done, 3-5
# REANNOT.FAIL written by 59, 99 memguard kill (FAIL "memguard" written here), other = crash (FAIL "crash" after
# the retry). Per-species logs OUT/logs/<code>.log, exit codes OUT/logs/rc.txt.
# Usage: bash scripts/59_run.sh [NJ 4] [T 16] [OUT work/addon/A1_reannot]
cd ~/deepsea_conv || exit 1
echo 10 > /proc/self/autogroup 2>/dev/null
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
NJ=${1:-4}; T=${2:-16}; OUT=$(readlink -f ${3:-work/addon/A1_reannot}); mkdir -p $OUT/logs
PY=envs/dc/bin/python; RS=envs/rphy/bin/Rscript
log() { echo "$(date +%F_%T) $*" >> $OUT/logs/run.log; }
die() { log "$1 failed ($2)"; echo "$(date +%F_%T) $1" > $OUT/A1R.FAIL; exit 1; }
nd() { echo "DONE $(ls $OUT/*/REANNOT.DONE 2>/dev/null | wc -l) FAIL $(ls $OUT/*/REANNOT.FAIL 2>/dev/null | wc -l)"; }
rm -f $OUT/A1R.FAIL
export OUT T PY
log "start NJ $NJ T $T; md5 59 $(md5sum scripts/59_route_reannot.py | cut -c1-32) 59b $(md5sum scripts/59b_reannot_refit.R | cut -c1-32) 30 $(md5sum scripts/30_physchem.py | cut -c1-32)"
if [ ! -e $OUT/SELECT.DONE ]; then
  $PY scripts/59_route_reannot.py --select --out $OUT >> $OUT/logs/select.log 2>&1
  [ -e $OUT/SELECT.DONE ] || die select logs/select.log
fi
CODES=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "code") c = i; next } c { print $c }' $OUT/selection.tsv)
[ -n "$CODES" ] || die select "selection.tsv has no code"
log "selected $(echo $CODES | wc -w) species"
one() {   # CODE: one re-annotation unless DONE / FAIL exists; rc -> logs/rc.txt; memguard kill -> REANNOT.FAIL
  local c=$1 d=$OUT/$1 rc
  if [ -e $d/REANNOT.DONE ] || [ -e $d/REANNOT.FAIL ]; then return 0; fi
  nice -n 10 bash scripts/memguard.sh 48 $PY scripts/59_route_reannot.py --one $c --threads $T --out $OUT >> $OUT/logs/$c.log 2>&1
  rc=$?
  echo -e "$(date +%F_%T)\t$c\t$rc" >> $OUT/logs/rc.txt
  if [ $rc = 99 ]; then
    mkdir -p $d; rm -f $d/genome.fa $d/genome.fa.tmp $d/genome.fa.fai $d/ref.faa
    echo "{\"code\": \"$c\", \"why\": \"memguard\", \"rc\": 99, \"ended\": \"$(date '+%F %T')\"}" > $d/REANNOT.FAIL
  fi
  return 0
}
export -f one
echo "$CODES" | xargs -P $NJ -I{} bash -c 'one "$@"' _ {}
log "pass 1: $(nd)"
for c in $CODES; do   # one serial retry of memguard kills and crashes (neither DONE nor FAIL)
  d=$OUT/$c
  [ -e $d/REANNOT.DONE ] && continue
  if [ -e $d/REANNOT.FAIL ]; then
    grep -q '"why": "memguard"' $d/REANNOT.FAIL || continue
    rm -f $d/REANNOT.FAIL; log "retry $c (memguard)"
  else log "retry $c (crash)"; fi
  one $c
  if [ ! -e $d/REANNOT.DONE ] && [ ! -e $d/REANNOT.FAIL ]; then
    rc=$(awk -F'\t' -v c=$c '$2 == c { r = $3 } END { print r }' $OUT/logs/rc.txt)
    rm -f $d/genome.fa $d/genome.fa.tmp $d/genome.fa.fai $d/ref.faa
    echo "{\"code\": \"$c\", \"why\": \"crash\", \"rc\": ${rc:-null}, \"ended\": \"$(date '+%F %T')\"}" > $d/REANNOT.FAIL
    log "crash $c rc $rc"
  fi
done
log "pass 2: $(nd)"
$PY scripts/59_route_reannot.py --collect --out $OUT >> $OUT/logs/collect.log 2>&1
[ -e $OUT/COLLECT.DONE ] || die collect logs/collect.log
log "collect: $(cat $OUT/COLLECT.DONE)"
NCORES=8 nice -n 10 bash scripts/memguard.sh 48 $RS scripts/59b_reannot_refit.R $OUT >> $OUT/logs/refit.log 2>&1
[ -e $OUT/REFIT.DONE ] || die refit logs/refit.log
date +%F_%T > $OUT/A1R.DONE; log "A1R.DONE"
exit 0
