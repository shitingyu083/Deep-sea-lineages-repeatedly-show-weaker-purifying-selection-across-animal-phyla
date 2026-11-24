#!/bin/bash
# Staged downloader (link is capped ~5 MB/s total, so order by what unblocks work):
#  phase 1: annotated accessions -> protein,cds,gff3 only (small); marks .done + .nogenome
#  phase 2: unannotated accessions -> genome (needed by miniprot)
#  phase 3: genomes for .nogenome accessions (only needed for genome-mode QC / SI)
# Skips .done, accessions another downloader is writing (pkg.zip touched <5 min),
# and uses an atomic mkdir lock. Writes DOWNLOAD.DONE when all phases finish.
# Usage: bash 05d_download_staged.sh data/ncbi/acc_list.tsv NPAR
set -u
B=~/deepsea_conv/envs/dc/bin
OUT=~/deepsea_conv/data/ncbi
list=$1; NP=${2:-8}

fetch() {  # acc include zipname
  acc=$1; inc=$2; z=$3; d="$OUT/$acc"
  for try in 1 2 3 4 5; do
    rm -f "$d/$z"
    if $B/datasets download genome accession "$acc" --include "$inc" \
         --filename "$d/$z" --no-progressbar > "$d/dl_$z.log" 2>&1; then
      (cd "$d" && unzip -oq "$z" && rm "$z")
      find "$d/ncbi_dataset/data/$acc" -type f \( -name "*.fna" -o -name "*.faa" -o -name "*.gff" \) -exec gzip -f {} \;
      return 0
    fi
    sleep $((try * 30))
  done
  return 1
}
busy() {  # another downloader active on this accession?
  local d="$OUT/$1"
  for z in pkg.zip pkg_g.zip; do
    [ -f "$d/$z" ] && [ -n "$(find "$d/$z" -mmin -5 2>/dev/null)" ] && return 0
  done
  return 1
}
dl() {
  phase=$1; acc=$2; ann=$3; d="$OUT/$acc"
  case $phase in
    1) [ "$ann" = none ] && return 0; [ -f "$d/.done" ] && return 0 ;;
    2) [ "$ann" != none ] && return 0; [ -f "$d/.done" ] && return 0 ;;
    3) [ -f "$d/.nogenome" ] || return 0 ;;
  esac
  busy $acc && { echo "SKIP_ACTIVE $phase $acc"; return 0; }
  mkdir -p "$d"; mkdir "$d/.lockd" 2>/dev/null || { echo "SKIP_LOCK $phase $acc"; return 0; }
  case $phase in
    1) if fetch $acc protein,cds,gff3 pkg.zip; then touch "$d/.nogenome" "$d/.done"; echo "OK1 $acc"
       else echo "FAIL $acc" >> "$OUT/failed.txt"; fi ;;
    2) if fetch $acc genome pkg.zip; then touch "$d/.done"; echo "OK2 $acc"
       else echo "FAIL $acc" >> "$OUT/failed.txt"; fi ;;
    3) if fetch $acc genome pkg_g.zip; then rm -f "$d/.nogenome"; echo "OK3 $acc"
       else echo "GFAIL $acc" >> "$OUT/failed.txt"; fi ;;
  esac
  rmdir "$d/.lockd"
}
export -f fetch busy dl; export OUT B
waitall() {  # wait for inherited downloads skipped as active
  while cut -f1 "$list" | while read a; do busy $a && echo x; done | grep -q x; do sleep 60; done
}
for ph in 1 2; do
  sed "s/^/$ph\t/" "$list" | xargs -P $NP -L 1 bash -c 'dl "$0" "$1" "$2"'
done
waitall
# retry phase 2 for anything an inherited downloader left unfinished
sed "s/^/2\t/" "$list" | xargs -P $NP -L 1 bash -c 'dl "$0" "$1" "$2"'
sed "s/^/3\t/" "$list" | xargs -P $NP -L 1 bash -c 'dl "$0" "$1" "$2"'
echo ALLDONE > "$OUT/DOWNLOAD.DONE"
