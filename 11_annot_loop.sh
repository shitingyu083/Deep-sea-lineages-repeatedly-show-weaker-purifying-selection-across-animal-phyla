#!/bin/bash
# Annotation daemon:
#  - ncbi_annot species: run 08_prep_ncbi_annot.py as soon as the download is .done
#  - miniprot species: start only after ALL ncbi_annot species of the phylum are prepped (best refs),
#    then run 10_miniprot_annot.py NJ in parallel, T threads each.
# 2026-10-19: continuous queue. The old loop ran batches of NJ and waited for the slowest species of each
#   batch (one large genome idled the other slots). Now every eligible species goes into one xargs -P NJ
#   queue, a species runs under flock on work/miniprot/<code>/.lock (held species are skipped), and the
#   loop holds work/annot_loop.lock so a watchdog restart cannot start a second copy.
# Usage: nohup bash scripts/11_annot_loop.sh 8 24 > logs/annot_loop.log 2>&1 &
R=~/deepsea_conv; cd $R
P=~/anaconda3/envs/py310env/bin/python
NJ=${1:-6}; T=${2:-24}
mkdir -p logs/miniprot work/miniprot
exec 8> work/annot_loop.lock
flock -n 8 || { echo "$(date +%F_%T) another annotation loop holds work/annot_loop.lock; exiting"; exit 0; }
export R P T
mp_one() {
  code=$1
  [ -s $R/data/proteomes/$code.faa ] && return 0
  [ -f $R/work/miniprot/$code/.fail ] && return 0
  mkdir -p $R/work/miniprot/$code
  exec 9> $R/work/miniprot/$code/.lock
  flock -n 9 || return 0
  [ -s $R/data/proteomes/$code.faa ] && return 0
  [ -f $R/work/miniprot/$code/.fail ] && return 0
  $P $R/scripts/10_miniprot_annot.py $code $T > $R/logs/miniprot/$code.log 2>&1 \
    || touch $R/work/miniprot/$code/.fail
}
export -f mp_one
while true; do
  pending_annot=0
  pend_phyla=" "
  while IFS=$'\t' read -r code sp acc src role phylum rest; do
    [ "$src" = ncbi_annot ] || continue
    [ -s data/proteomes/$code.faa ] && continue
    [ -f work/annotfail/$code ] && continue
    grep -q "FAIL $acc" data/ncbi/failed.txt 2>/dev/null && continue
    if [ -f data/ncbi/$acc/.done ]; then
      $P scripts/08_prep_ncbi_annot.py $code $acc >> logs/annot_loop.log 2>&1
      # unusable annotation (missing files / too few genes) -> fall back to miniprot
      [ -s data/proteomes/$code.faa ] || { mkdir -p work/annotfail; touch work/annotfail/$code; }
    else
      pending_annot=$((pending_annot+1)); pend_phyla="$pend_phyla$phylum "
    fi
  done < <(tail -n +2 meta/species_codes.tsv)
  # miniprot refs are phylum-restricted: a species may start once its phylum's annotated references are all
  # prepped. Species held by a running job (lock) are left out of this round.
  todo=(); nheld=0
  while IFS=$'\t' read -r code sp acc src role phylum rest; do
    [ "$src" = miniprot ] || [ -f work/annotfail/$code ] || continue
    [[ "$pend_phyla" == *" $phylum "* ]] && continue
    [ -s data/proteomes/$code.faa ] && continue
    [ -f work/miniprot/$code/.fail ] && continue
    [ -f data/ncbi/$acc/.done ] || continue
    if [ -e work/miniprot/$code/.lock ] && ! flock -n work/miniprot/$code/.lock true; then
      nheld=$((nheld+1)); continue
    fi
    todo+=($code)
  done < <(tail -n +2 meta/species_codes.tsv)
  n_mp_left=$( (awk -F'\t' 'NR>1 && $4=="miniprot"' meta/species_codes.tsv | cut -f1; ls work/annotfail 2>/dev/null) | \
     while read c; do [ -s data/proteomes/$c.faa ] || [ -f work/miniprot/$c/.fail ] || echo $c; done | wc -l)
  echo "$(date +%F_%T) pending_annot=$pending_annot pend_phyla=[${pend_phyla}] mp_todo=${#todo[@]} held=$nheld mp_left=$n_mp_left"
  if [ ${#todo[@]} -gt 0 ]; then
    printf '%s\n' "${todo[@]}" | xargs -P $NJ -I{} bash -c 'mp_one {}'
  elif [ $pending_annot = 0 ] && [ $n_mp_left = 0 ]; then
    break
  elif [ $nheld -gt 0 ]; then
    sleep 60
  else
    sleep 300
  fi
done
echo ANNOTDONE > work/ANNOT.DONE
