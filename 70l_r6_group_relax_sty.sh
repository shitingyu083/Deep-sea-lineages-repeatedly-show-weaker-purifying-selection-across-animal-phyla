#!/bin/bash
# R6-2 fits (ADDENDUM_R6_2.md md5 0e9d4fe6): RELAX Group mode, each group {Li} vs reference BG, on the R6-1 U-arm
# alignment; 2 fits per HOG (real, null1). Budget: after 300 completed fits, projected wall > 30 h -> CUT flag,
# later jobs outside cut100.txt are skipped (done ones outside the cut are not analysed); at 30 h wall -> STOP flag,
# remaining jobs skipped. Output R6_2/sel/<c>/<hog>.<set>.grp.relax.json; logs R6_2/sel/logs.
# Usage: bash 70l_r6_group_relax_sty.sh [NJ 64]   (run under scripts/memguard.sh)
R=~/deepsea_conv; O=$R/work/addon/R6; G=$O/R6_2; NJ=${1:-64}
[ -f $G/PREP.DONE ] || { echo "no PREP.DONE"; exit 2; }
mkdir -p $G/sel/logs
T0=$(date +%s); NT=$(wc -l < $G/jobs.txt)
echo "$(date +%F_%T) start NJ $NJ jobs $NT t0 $T0" >> $G/sel/logs/run.log
fit() {
  c=$1; h=$2; s=$3; d=$G/sel/$c; j=$d/$h.$s.grp.relax.json
  [ -s $j ] && return 0
  [ -f $G/STOP ] && return 0
  [ $(( $(date +%s) - T0 )) -gt 108000 ] && { touch $G/STOP; echo "$(date +%F_%T) STOP 30h" >> $G/sel/logs/run.log; return 0; }
  [ -f $G/CUT ] && ! grep -qx "$c $h" $G/cut100.txt && return 0
  mkdir -p $d
  t1=$(date +%s)
  nice -n 19 $R/envs/dc/bin/hyphy CPU=1 relax --alignment $O/aln/$c/$h.U.fas --tree $G/trees/$c/$h.$s.grp.nwk \
    --mode "Group mode" --reference-group BG --models Minimal --output $j.tmp > $G/sel/logs/$c.$h.$s.log 2>&1
  rc=$?
  if [ $rc -eq 0 ] && [ -s $j.tmp ]; then mv $j.tmp $j; else rm -f $j.tmp; fi
  echo "$c $h $s $rc $(( $(date +%s) - t1 ))" >> $G/sel/logs/rc.txt
  n=$(wc -l < $G/sel/logs/rc.txt)
  if [ $n -ge 300 ] && mkdir $G/sel/logs/PROJ.lock 2>/dev/null; then
    e=$(( $(date +%s) - T0 )); p=$(( e * NT / 300 ))
    echo "$(date +%F_%T) projection after 300: elapsed $e s, projected $p s" >> $G/sel/logs/run.log
    [ $p -gt 108000 ] && { touch $G/CUT; echo "$(date +%F_%T) CUT to cut100.txt" >> $G/sel/logs/run.log; }
  fi
}
export -f fit; export O R G T0 NT
xargs -P $NJ -L 1 bash -c 'fit $0 $1 $2' < $G/jobs.txt
echo "$(date +%F_%T) end ok $(ls $G/sel/*/*.grp.relax.json 2>/dev/null | wc -l)" >> $G/sel/logs/run.log
touch $G/RELAX.END
exit 0
