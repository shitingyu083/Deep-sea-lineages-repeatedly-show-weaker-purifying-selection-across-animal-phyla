#!/bin/bash
# DIAGNOSTIC ONLY (DEVIATIONS_r6, 2026-11-05 entry): run-to-run reproducibility of the UNPATCHED
# HyPhy 2.5.96 RELAX.bf on the 24 previously-ok equivalence jobs. Same command line as 70d/70r,
# same inputs. If two unpatched runs differ as much as patched-vs-original did, the gate failure
# is optimizer noise, not the patch. Nothing here enters any R6 estimate.
# Outputs: R6/refit/<OUTSUB>/ (never R6/sel/).
# Usage: bash 70s_r6_repro_orig_sty.sh JOBLIST OUTSUB [NJ 12]
R=~/deepsea_conv; O=$R/work/addon/R6; F=$O/refit; L=${1:?joblist}; S=${2:?outsub}; NJ=${3:-12}
BF=$R/envs/dc/share/hyphy/TemplateBatchFiles/SelectionAnalyses/RELAX.bf
mkdir -p $F/$S/logs
md5sum $BF > $F/$S/bf_md5.txt
echo "$(date +%F_%T) start $S NJ $NJ jobs $(wc -l < $L) UNPATCHED" >> $F/$S/logs/run.log
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
