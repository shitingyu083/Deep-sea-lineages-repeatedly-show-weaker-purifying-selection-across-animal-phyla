#!/bin/bash
# R6-1 re-fit of the RELAX O/U jobs that failed in 70d (DEVIATIONS_r6, 2026-11-04 entry).
# Cause of failure: HyPhy 2.5.96 RELAX crashes in ancestral.build (RecoverAncestralSequences,
# "Transition matrix not computed") AFTER the alternative model is fitted. That loop only stores
# the substitution map in the JSON; it does not enter K or p. The patched copy RELAX_noanc.bf
# comments out exactly that loop; everything else, including the command line, is as in 70d.
# Inputs are unchanged (same aln/ and trees/). Outputs go to R6/refit/, never to R6/sel/.
# Usage: bash 70r_r6_refit_noanc_sty.sh JOBLIST OUTSUB [NJ 16]
#   JOBLIST lines: clade hog set arm      OUTSUB: subdirectory of R6/refit (e.g. equiv, fail)
R=~/deepsea_conv; O=$R/work/addon/R6; F=$O/refit; L=${1:?joblist}; S=${2:?outsub}; NJ=${3:-16}
SRC=$R/envs/dc/share/hyphy/TemplateBatchFiles/SelectionAnalyses/RELAX.bf
BF=$R/envs/dc/share/hyphy/TemplateBatchFiles/SelectionAnalyses/RELAX_noanc_r6.bf
mkdir -p $F/$S/logs
if [ ! -f $BF ]; then
  # comment out the 5-line ancestral loop (for ... ancestral.build ... DeleteObject ... })
  awk '/for \(_partition_, _selection_; in; relax.selected_branches\)/ && !done {skip=5; done=1}
       skip>0 {print "// R6-noanc " $0; skip--; next} {print}' $SRC > $BF
fi
md5sum $SRC $BF > $F/$S/bf_md5.txt
diff $SRC $BF > $F/$S/bf_patch.diff
echo "$(date +%F_%T) start $S NJ $NJ jobs $(wc -l < $L)" >> $F/$S/logs/run.log
fit() {
  c=$1; h=$2; s=$3; a=$4; d=$F/$S/$c; j=$d/$h.$s.$a.relax.json
  [ -s $j ] && return 0
  mkdir -p $d
  t0=$(date +%s)
  nice -n 19 $R/envs/dc/bin/hyphy CPU=1 $BF --alignment $O/aln/$c/$h.$a.fas --tree $O/trees/$c/$h.$s.nwk \
    --test FG --reference BG --models Minimal --output $j.tmp > $F/$S/logs/$c.$h.$s.$a.log 2>&1
  rc=$?
  if [ $rc -eq 0 ] && [ -s $j.tmp ]; then mv $j.tmp $j; else rm -f $j.tmp; fi
  echo "$c $h $s $a $rc $(( $(date +%s) - t0 ))" >> $F/$S/logs/rc.txt
}
export -f fit; export O R F S BF
xargs -P $NJ -L 1 bash -c 'fit $0 $1 $2 $3' < $L
echo "$(date +%F_%T) end ok $(ls $F/$S/*/*.relax.json 2>/dev/null | wc -l)" >> $F/$S/logs/run.log
touch $F/$S/END
