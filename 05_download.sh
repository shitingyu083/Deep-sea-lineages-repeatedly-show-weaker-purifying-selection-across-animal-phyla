#!/bin/bash
# Download panel genomes from NCBI (genome + annotation where present).
# Usage: bash 05_download.sh meta/panel_v1.tsv  (run from ~/deepsea_conv)
set -u
B=~/deepsea_conv/envs/dc/bin
OUT=~/deepsea_conv/data/ncbi
mkdir -p "$OUT"
panel=$1
tail -n +2 "$panel" | awk -F'\t' '{print ($9=="none"?1:0)"\t"$8"\t"$9}' | sort -u | sort -k1,1n -s | cut -f2- > "$OUT/acc_list.tsv"

dl() {
  acc=$1; ann=$2; OUT=$3; B=$4
  d="$OUT/$acc"
  [ -f "$d/.done" ] && return 0
  mkdir -p "$d"
  inc="genome"
  [ "$ann" != "none" ] && inc="genome,protein,cds,gff3"
  for try in 1 2 3; do
    if $B/datasets download genome accession "$acc" --include "$inc" \
         --filename "$d/pkg.zip" --no-progressbar > "$d/dl.log" 2>&1; then
      (cd "$d" && unzip -oq pkg.zip && rm pkg.zip)
      find "$d/ncbi_dataset/data/$acc" -type f \( -name "*.fna" -o -name "*.faa" -o -name "*.gff" \) -exec gzip -f {} \;
      touch "$d/.done"; return 0
    fi
    sleep $((try * 30))
  done
  echo "FAIL $acc" >> "$OUT/failed.txt"
}
export -f dl
cat "$OUT/acc_list.tsv" | xargs -P 6 -L 1 bash -c 'dl "$0" "$1" '"$OUT $B"
echo ALLDONE > "$OUT/DOWNLOAD.DONE"
