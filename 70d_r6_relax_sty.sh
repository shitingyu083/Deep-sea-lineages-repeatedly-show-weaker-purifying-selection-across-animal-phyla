#!/bin/bash
# R6-1 step d (PREREG_r6.md): RELAX on the paired O/U alignments built by 70c, 4 fits per HOG
# (real-O, real-U, null1-O, null1-U), interleaved by HOG; production command (Minimal, CPU=1).
# Output work/addon/R6/sel/<clade>/<hog>.<set>.<arm>.relax.json; existing JSONs skipped; logs in R6/sel/logs.
# Usage: bash 70d_r6_relax_sty.sh [NJ 64]   (run under scripts/memguard.sh)
R=~/deepsea_conv; O=$R/work/addon/R6; NJ=${1:-64}
[ -f $O/BUILD.DONE ] || { echo "no BUILD.DONE"; exit 2; }
mkdir -p $O/sel/logs
awk -F'\t' 'NR>1 && $7=="ok"{for(s=1;s<=2;s++)for(a=1;a<=2;a++)print $1, $2, (s==1?"real":"null1"), (a==1?"O":"U")}' \
  $O/build_report.tsv > $O/sel/jobs.txt
echo "$(date +%F_%T) start NJ $NJ jobs $(wc -l < $O/sel/jobs.txt)" >> $O/sel/logs/run.log
fit() {
  c=$1; h=$2; s=$3; a=$4; d=$O/sel/$c; j=$d/$h.$s.$a.relax.json
  [ -s $j ] && return 0
  mkdir -p $d
  t0=$(date +%s)
  nice -n 19 $R/envs/dc/bin/hyphy CPU=1 relax --alignment $O/aln/$c/$h.$a.fas --tree $O/trees/$c/$h.$s.nwk \
    --test FG --reference BG --models Minimal --output $j.tmp > $O/sel/logs/$c.$h.$s.$a.log 2>&1
  rc=$?
  if [ $rc -eq 0 ] && [ -s $j.tmp ]; then mv $j.tmp $j; else rm -f $j.tmp; fi
  echo "$c $h $s $a $rc $(( $(date +%s) - t0 ))" >> $O/sel/logs/rc.txt
}
export -f fit; export O R
xargs -P $NJ -L 1 bash -c 'fit $0 $1 $2 $3' < $O/sel/jobs.txt
echo "$(date +%F_%T) end ok $(ls $O/sel/*/*.relax.json 2>/dev/null | wc -l)" >> $O/sel/logs/run.log
touch $O/RELAX.END
exit 0
