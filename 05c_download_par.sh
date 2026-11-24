#!/bin/bash
# Higher-parallelism downloader (bandwidth is capped per connection, not total).
# Same output layout as 05_download.sh. Skips accessions that are .done or that
# another downloader is actively writing (pkg.zip touched in the last 5 min).
# Uses an atomic mkdir lock so two instances never take the same accession.
# Usage: bash 05c_download_par.sh data/ncbi/acc_list.tsv NPAR
set -u
B=~/deepsea_conv/envs/dc/bin
OUT=~/deepsea_conv/data/ncbi
list=$1; NP=${2:-16}

dl() {
  acc=$1; ann=$2; OUT=$3; B=$4
  d="$OUT/$acc"
  [ -f "$d/.done" ] && return 0
  if [ -f "$d/pkg.zip" ] && [ -n "$(find "$d/pkg.zip" -mmin -5 2>/dev/null)" ]; then
    echo "SKIP_ACTIVE $acc"; return 0
  fi
  mkdir -p "$d"
  mkdir "$d/.lockc" 2>/dev/null || { echo "SKIP_LOCK $acc"; return 0; }
  inc="genome"
  [ "$ann" != "none" ] && inc="genome,protein,cds,gff3"
  for try in 1 2 3 4 5; do
    rm -f "$d/pkg.zip"
    if $B/datasets download genome accession "$acc" --include "$inc" \
         --filename "$d/pkg.zip" --no-progressbar > "$d/dl.log" 2>&1; then
      (cd "$d" && unzip -oq pkg.zip && rm pkg.zip)
      find "$d/ncbi_dataset/data/$acc" -type f \( -name "*.fna" -o -name "*.faa" -o -name "*.gff" \) -exec gzip -f {} \;
      touch "$d/.done"; rmdir "$d/.lockc"; echo "OK $acc"; return 0
    fi
    sleep $((try * 30))
  done
  rmdir "$d/.lockc"
  echo "FAIL $acc" >> "$OUT/failed.txt"
}
export -f dl
cat "$list" | xargs -P $NP -L 1 bash -c 'dl "$0" "$1" '"$OUT $B"
# wait for downloads inherited from the previous downloader (skipped as active)
while cut -f1 "$list" | while read a; do [ -f "$OUT/$a/pkg.zip" ] && echo x; done | grep -q x; do sleep 60; done
echo ALLDONE > "$OUT/DOWNLOAD.DONE"
