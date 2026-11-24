#!/bin/bash
# Retry accessions whose download failed because two downloaders wrote the same
# pkg.zip (orphaned 05c retries vs 05d, 2026-10-17). Sequential, one connection each.
# Annotated -> protein,cds,gff3 + .nogenome (05d phase 3 adds the genome later).
# Removes the accession's FAIL lines on success.
# Usage: bash 05e_retry_failed.sh ACC [ACC ...]
set -u
B=~/deepsea_conv/envs/dc/bin; OUT=~/deepsea_conv/data/ncbi
for acc in "$@"; do
  d=$OUT/$acc
  if [ -f $d/.done ] && ls $d/ncbi_dataset/data/$acc/protein.faa.gz >/dev/null 2>&1; then
    sed -i "/FAIL $acc$/d" $OUT/failed.txt; echo "ALREADY $acc"; continue
  fi
  rm -rf $d; mkdir -p $d
  ok=0
  for try in 1 2 3 4 5; do
    if $B/datasets download genome accession $acc --include protein,cds,gff3 \
         --filename $d/pkg_r.zip --no-progressbar > $d/dl_r.log 2>&1 && (cd $d && unzip -oq pkg_r.zip); then
      rm -f $d/pkg_r.zip
      find $d/ncbi_dataset/data/$acc -type f \( -name "*.fna" -o -name "*.faa" -o -name "*.gff" \) -exec gzip -f {} \;
      touch $d/.nogenome $d/.done; ok=1; break
    fi
    sleep $((try * 30))
  done
  if [ $ok = 1 ]; then sed -i "/FAIL $acc$/d" $OUT/failed.txt; echo "OK $acc"; else echo "STILLFAIL $acc"; fi
done
