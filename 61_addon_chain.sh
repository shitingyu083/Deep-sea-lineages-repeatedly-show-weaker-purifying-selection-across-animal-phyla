#!/bin/bash
# Add-on analysis chain (2026-10-26; PREREG_addon.md, work/addon/FROZEN_addon.txt). Five parallel streams, every heavy
# step under memguard 48 GB and nice 10; flags and rc in work/addon/logs/chain61.log, per-step logs work/addon/logs.
#  A  53a genome GC (8 workers) -> 53 route / GC / quality confounds (A1)                         -> work/addon/A1
#  B  50b MCMCTree backbone (W = work/tree_v3, the loci of the full-tree run) -> BB.DONE -> 50c interpolation
#     -> 58 timing with the ROBUST / POINT trees                                                  -> work/addon/A6
#  C  58 on the LSD2 tree only (preliminary, ROBUST / POINT absent) -> 55 origins (A2) -> 56 LOCO lab (A3),
#     lab_strict (A3_strict)
#  D  54a/54b rescued (A4) and norescue (A4/norescue, SI) -> 57 on the 17:48 annotation (A5_prelim, SI) -> 60 (T1)
#  E  59_run.sh re-annotation control (A1_reannot)
#  F  after work/eggnog/EMAPPER.DONE: 47 per clade -> work/func_v2/<C>, 48 -> work/func_v2/combined, 57 -> work/addon/A5
# ADDON.DONE when every stream has ended (each stream writes STREAM_<X>.END with its failures).
# Usage: bash scripts/61_addon_chain.sh
cd ~/deepsea_conv || exit 1
echo 10 > /proc/self/autogroup 2>/dev/null
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=envs/dc/bin/python; RS=envs/rphy/bin/Rscript
A=work/addon; L=$A/logs; mkdir -p $L; rm -f $A/STREAM_*.END $A/ADDON.DONE
log() { echo "$(date +%F_%T) $*" >> $L/chain61.log; }
# run NAME CMD...: memguard + nice, rc logged; returns rc
run() { local n=$1; shift; log "$n start"; nice -n 10 bash scripts/memguard.sh 48 "$@" > $L/$n.log 2>&1; local rc=$?
        log "$n rc=$rc"; [ $rc = 0 ] || FAILS="$FAILS $n"; return $rc; }
fin() { echo "$(date +%F_%T) fails:${FAILS:- none}" > $A/STREAM_$1.END; log "stream $1 end; fails:${FAILS:- none}"; }
log "start; FROZEN md5 $(md5sum $A/FROZEN_addon.txt | cut -c1-32)"
md5chk() { local bad=0; while read m f; do case "$m" in \#*|md5sum:) continue;; esac
    if [ -e scripts/$f ]; then [ "$(md5sum scripts/$f | cut -c1-32)" = "$m" ] || { log "md5 changed: $f"; bad=1; }
    elif [ -e $f ]; then [ "$(md5sum $f | cut -c1-32)" = "$m" ] || { log "md5 changed: $f"; bad=1; }; fi
  done < $A/FROZEN_addon.txt; return $bad; }
md5chk || { log "FROZEN mismatch, stop"; exit 1; }

sA() { FAILS=
  run 53a $PY scripts/53a_genome_gc.py 8 && \
  NCORES=16 REF44=work/final_v3/pgls_controls/controls_species.tsv REF32=work/final_v3/pgls/pgls_gene.tsv \
    run 53 $RS scripts/53_route_confound.R work/final_v3/date/lsd.timetree.nwk $A/A1
  fin A; }

sB() { FAILS=
  run 50b bash scripts/50b_mcmctree_bb.sh work/tree_v3 work/mcmctree_bb
  if [ -e work/mcmctree_bb/BB.DONE ]; then
    run 50c $RS scripts/50c_bb_interp.R work/mcmctree_bb 8 && \
    NCORES=8 run 58 $RS scripts/58_origin_timing.R work/final_v3/date/lsd.timetree.nwk $A/A6
  else log "BB.DONE missing, 50c / 58 full not run"; FAILS="$FAILS BB"; fi
  fin B; }

sC() { FAILS=
  NCORES=8 ROBUST=/nonexistent_sty POINT=/nonexistent_sty run 58lsd $RS scripts/58_origin_timing.R \
    work/final_v3/date/lsd.timetree.nwk $A/A6_lsd
  NCORES=8 run 55 $RS scripts/55_origin_contrasts.R work/final_v3/date/lsd.timetree.nwk $A/A2
  run 56 $RS scripts/56_loco_predict.R work/final_v3/date/lsd.timetree.nwk $A/A3
  LABEL_COL=lab_strict run 56s $RS scripts/56_loco_predict.R work/final_v3/date/lsd.timetree.nwk $A/A3_strict
  fin C; }

sD() { FAILS=
  run 54a $PY scripts/54a_sel_features.py --out $A/A4 --np 4 && \
    run 54b $RS scripts/54b_sel_confound.R $A/A4/sel_features.tsv $A/A4 work/combine25c/final_main.combined.tsv
  run 54an $PY scripts/54a_sel_features.py --sub norescue --out $A/A4/norescue --np 4 && \
    run 54bn $RS scripts/54b_sel_confound.R $A/A4/norescue/sel_features.tsv $A/A4/norescue work/combine25c/si_norescue.combined.tsv
  run 57p $PY scripts/57_apriori_sets.py --func work/func --out $A/A5_prelim
  run 60 $PY scripts/60_dose_series.py --out $A/T1
  fin D; }

sE() { FAILS=
  run 59 bash scripts/59_run.sh 2 12 $A/A1_reannot
  [ -e $A/A1_reannot/A1R.DONE ] || FAILS="$FAILS A1R"
  fin E; }

sF() { FAILS=; local n=0
  until [ -e work/eggnog/EMAPPER.DONE ]; do sleep 1800; n=$((n+1)); [ $n -gt 672 ] && break; done   # <= 14 d
  if [ -e work/eggnog/EMAPPER.DONE ]; then
    log "EMAPPER.DONE seen"; mkdir -p work/func_v2/logs
    for C in Chondrichthyes Porifera Arthropoda Echinodermata Cnidaria Annelida Mollusca Actinopteri_m; do
      run 47_$C $PY scripts/47_func_enrich.py $C --out work/func_v2/$C; done
    run 48 $PY scripts/48_func_combine.py --func work/func_v2 && touch work/func_v2/FUNC.DONE
    run 57 $PY scripts/57_apriori_sets.py --func work/func_v2 --out $A/A5
  else log "EMAPPER.DONE not seen in 14 d"; FAILS="$FAILS EMAPPER"; fi
  fin F; }

sA & sB & sC & sD & sE & sF &
wait
touch $A/ADDON.DONE; log "ADDON.DONE"
exit 0
