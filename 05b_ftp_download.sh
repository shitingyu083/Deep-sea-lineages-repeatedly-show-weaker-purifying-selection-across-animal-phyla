#!/bin/bash
# Faster NCBI download via FTP-over-HTTPS + aria2c (many files in parallel).
# Writes flat files into data/ncbi/<acc>/ : genomic.fna.gz [protein.faa.gz genomic.gff.gz cds.fna.gz]
# Skips accessions that already have .done (from 05_download.sh or this script).
# Usage (from ~/deepsea_conv): bash scripts/05b_ftp_download.sh data/ncbi/acc_list.tsv
set -u
A=~/deepsea_conv/envs/dc/bin/aria2c
OUT=~/deepsea_conv/data/ncbi
list=$1
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all
inp=$OUT/aria2_input.txt; : > "$inp"
: > "$OUT/ftp_pending.txt"
while IFS=$'\t' read -r acc ann; do
  d=$OUT/$acc
  [ -f "$d/.done" ] && continue
  p=${acc:0:3}/${acc:4:3}/${acc:7:3}/${acc:10:3}
  dir=$(curl -s --retry 3 "$BASE/$p/" | grep -o "\"${acc}_[^\"/]*/\"" | tr -d '"/' | head -1)
  if [ -z "$dir" ]; then echo "NOFTP $acc" >> "$OUT/failed.txt"; continue; fi
  rm -rf "$d"; mkdir -p "$d"
  echo -e "$acc\t$ann\t$dir" >> "$OUT/ftp_pending.txt"
  u=$BASE/$p/$dir/$dir
  printf '%s\n  dir=%s\n  out=genomic.fna.gz\n' "${u}_genomic.fna.gz" "$d" >> "$inp"
  if [ "$ann" != "none" ]; then
    printf '%s\n  dir=%s\n  out=protein.faa.gz\n' "${u}_protein.faa.gz" "$d" >> "$inp"
    printf '%s\n  dir=%s\n  out=genomic.gff.gz\n' "${u}_genomic.gff.gz" "$d" >> "$inp"
    printf '%s\n  dir=%s\n  out=cds.fna.gz\n' "${u}_cds_from_genomic.fna.gz" "$d" >> "$inp"
  fi
done < "$list"
echo "queued $(wc -l < "$OUT/ftp_pending.txt") assemblies"
$A -i "$inp" -j 24 -x 4 -s 4 -k 2M --file-allocation=none --max-tries=10 --retry-wait=20 \
   --auto-file-renaming=false --allow-overwrite=true --console-log-level=warn --summary-interval=0
# verify gzip integrity, then mark done
while IFS=$'\t' read -r acc ann dir; do
  d=$OUT/$acc; ok=1
  for f in "$d"/*.gz; do gzip -t "$f" 2>/dev/null || { ok=0; echo "BADGZ $acc $(basename "$f")" >> "$OUT/failed.txt"; }; done
  [ -s "$d/genomic.fna.gz" ] || ok=0
  [ $ok = 1 ] && touch "$d/.done"
done < "$OUT/ftp_pending.txt"
echo ALLDONE > "$OUT/DOWNLOAD.DONE"
