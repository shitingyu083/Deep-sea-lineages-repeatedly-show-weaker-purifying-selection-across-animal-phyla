#!/bin/bash
# R6-1 step b (PREREG_r6.md): uniform miniprot re-annotation of the annotation-route species in the R6 sample
# (work/addon/R6/reannot_codes.txt), with the unchanged A1 route script 59_route_reannot.py --one (same reference
# rule, own species excluded, miniprot --outn=3, same filters). NJ jobs x T threads, memguard 48 GB, nice 10.
# A code with REANNOT.DONE / REANNOT.FAIL is skipped; crashes (rc not 0/3/4/5/99) get one serial retry.
# Usage: bash scripts/70b_r6_reannot_sty.sh [NJ 4] [T 16]
cd ~/deepsea_conv || exit 1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
NJ=${1:-4}; T=${2:-16}; OUT=$(readlink -f work/addon/R6/reannot); mkdir -p $OUT/logs
PY=envs/dc/bin/python; export OUT T PY
echo "$(date +%F_%T) start NJ $NJ T $T md5 59 $(md5sum scripts/59_route_reannot.py | cut -c1-32)" >> $OUT/logs/run.log
one() {
  local c=$1 d=$OUT/$1 rc
  if [ -e $d/REANNOT.DONE ] || [ -e $d/REANNOT.FAIL ]; then return 0; fi
  nice -n 10 bash scripts/memguard.sh 48 $PY scripts/59_route_reannot.py --one $c --threads $T --out $OUT >> $OUT/logs/$c.log 2>&1
  rc=$?; echo "$c $rc" >> $OUT/logs/rc.txt
  [ $rc = 99 ] && { mkdir -p $d; echo memguard > $d/REANNOT.FAIL; }
  return 0
}
export -f one
xargs -P $NJ -n 1 bash -c 'one "$0"' < work/addon/R6/reannot_codes.txt
for c in $(cat work/addon/R6/reannot_codes.txt); do
  [ -e $OUT/$c/REANNOT.DONE ] || [ -e $OUT/$c/REANNOT.FAIL ] || one $c
  [ -e $OUT/$c/REANNOT.DONE ] || [ -e $OUT/$c/REANNOT.FAIL ] || { mkdir -p $OUT/$c; echo crash > $OUT/$c/REANNOT.FAIL; }
done
echo "$(date +%F_%T) end DONE $(ls $OUT/*/REANNOT.DONE 2>/dev/null | wc -l) FAIL $(ls $OUT/*/REANNOT.FAIL 2>/dev/null | wc -l)" >> $OUT/logs/run.log
touch work/addon/R6/REANNOT.END
exit 0
