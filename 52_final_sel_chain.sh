#!/bin/bash
# final selection chain (2026-10-25): wait rescued 25b (7 clades + old fish) and fish_m SEL.DONE ->
# fish_m no-rescue 25b snapshot -> fish_m RELAX rescue (patched RELAX_noanc.bf) -> merge -> 25b ->
# 25c primary (PREREG work/combine25c) + SI variants -> Fig 4 enrichment (47, 8 clades) -> 48 combine.
cd ~/deepsea_conv
L=work/rescue/chain492.log; log() { echo "$(date +%F_%T) $*" >> $L; }
PY=envs/dc/bin/python; E=$PWD/envs/dc/bin; P=$PWD/scripts/hyphy_patch/SelectionAnalyses/RELAX_noanc.bf
export E P
log "start; waiting 25B_RESCUED.DONE and Actinopteri_m/SEL.DONE"
while [ ! -e work/rescue/25B_RESCUED.DONE ] || [ ! -e work/codon/Actinopteri_m/SEL.DONE ]; do sleep 600; done
C=Actinopteri_m; W=work/codon/$C
log "fish_m SEL.DONE seen; no-rescue 25b snapshot"
mkdir -p $W/norescue $W/sel_rescue
nice -n 10 $PY scripts/25b_collect_sel.py $C --nperm 10000 --out $W/norescue > $W/norescue/25b.log 2>&1
log "fish_m norescue 25b rc=$?"
one() { C=${1%%/*}; ht=${1#*/}; W=$PWD/work/codon/$C; o=$W/sel_rescue/$ht.relax.json
  [ -s $o ] && return 0
  nice -n 19 $E/hyphy CPU=1 $P --alignment $W/aln/${ht%%.null*}.fas --tree $W/trees/$ht.hyphy.nwk --test FG --reference BG \
     --models Minimal --output $o > $W/sel_rescue/$ht.relax.log 2>&1
  echo "$(date +%F_%T) $C $ht rc=$? size=$(stat -c %s $o 2>/dev/null)"; }
export -f one
find $W/sel -name '*.relax.json' -size 0 | while read j; do echo "$C/$(basename $j .relax.json)"; done > /tmp/c492_sty.list
log "fish_m rescue n=$(wc -l < /tmp/c492_sty.list) P md5 $(md5sum < $P | cut -c1-8)"
xargs -P 24 -I{} bash -c 'one {}' < /tmp/c492_sty.list >> work/rescue/rescue_fishm.log 2>&1
log "fish_m rescue done"
nice -n 10 $PY scripts/51_rescue_merge.py "$C" --apply --out work/rescue/fishm > work/rescue/merge_fishm.log 2>&1
log "fish_m merge rc=$?"
nice -n 10 $PY scripts/25b_collect_sel.py $C > $W/25b_rescued.log 2>&1
log "fish_m 25b rc=$?"
# 25c (frozen script; primary = 8 clades with fish_m; SI: no-rescue, fish old full-deep)
CL8="Chondrichthyes Porifera Arthropoda Echinodermata Cnidaria Annelida Mollusca Actinopteri_m"
a=(); b=(); c=()
for X in $CL8; do
  a+=("$X=work/codon/$X/25b_rescued.log"); b+=("$X=work/codon/$X/norescue/25b.log")
  [ $X = Actinopteri_m ] && c+=("Actinopteri=work/codon/Actinopteri/25b_rescued.log") || c+=("$X=work/codon/$X/25b_rescued.log")
done
O=work/combine25c
$PY scripts/25c_combine.py $O/final_main "${a[@]}" > $O/final_main.log 2>&1; log "25c main rc=$?"
$PY scripts/25c_combine.py $O/si_norescue "${b[@]}" > $O/si_norescue.log 2>&1; log "25c norescue rc=$?"
$PY scripts/25c_combine.py $O/si_fishold "${c[@]}" > $O/si_fishold.log 2>&1; log "25c fish-old rc=$?"
touch work/rescue/SEL_FINAL.DONE
# Fig 4 enrichment, confirmatory (PREREG_fig4 frozen before this point)
mkdir -p work/func/logs
for X in $CL8; do echo $X; done | xargs -P 4 -I{} bash -c "nice -n 10 $PY scripts/47_func_enrich.py {} > work/func/logs/47_{}.log 2>&1; echo \"\$(date +%F_%T) 47 {} rc=\$?\" >> $L"
nice -n 10 $PY scripts/48_func_combine.py > work/func/logs/48.log 2>&1; log "48 rc=$?"
touch work/func/FUNC.DONE; log "FUNC.DONE"
exit 0
