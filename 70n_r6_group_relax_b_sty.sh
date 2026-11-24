#!/bin/bash
# R6-2b fits (ADDENDUM_R6_4.md part 1): copy of 70l_r6_group_relax_sty.sh with OUT = R6_2b and the production
# original-route alignment work/codon/<c>/aln/<h>.fas. RELAX Group mode, {Li} vs reference BG; 2 fits per HOG.
# Budget: after 300 fits projected wall > 30 h -> CUT (cut100.txt); at 30 h wall -> STOP.
# Usage: bash 70n_r6_group_relax_b_sty.sh [NJ 64]   (run under scripts/memguard.sh)
R=~/deepsea_conv; O=$R/work/addon/R6; G=$O/R6_2b; NJ=${1:-64}
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
  nice -n 19 $R/envs/dc/bin/hyphy CPU=1 relax --alignment $R/work/codon/$c/aln/$h.fas --tree $G/trees/$c/$h.$s.grp.nwk \
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
