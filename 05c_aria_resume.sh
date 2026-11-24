#!/bin/bash
# resume aria2 on existing aria2_input.txt with moderate concurrency, then verify
OUT=~/deepsea_conv/data/ncbi
A=~/deepsea_conv/envs/dc/bin/aria2c
J=${1:-8}; X=${2:-2}
$A -i $OUT/aria2_input.txt -j $J -x $X -s $X -k 2M --file-allocation=none --max-tries=20 --retry-wait=30 \
   --connect-timeout=30 --timeout=120 --auto-file-renaming=false -c --console-log-level=warn --summary-interval=60 --show-console-readout=false
while IFS=$'\t' read -r acc ann dir; do
  d=$OUT/$acc; ok=1
  for f in "$d"/*.gz; do gzip -t "$f" 2>/dev/null || { ok=0; echo "BADGZ $acc $(basename "$f")" >> "$OUT/failed.txt"; }; done
  [ -s "$d/genomic.fna.gz" ] || ok=0
  [ $ok = 1 ] && touch "$d/.done"
done < $OUT/ftp_pending.txt
echo ALLDONE > $OUT/DOWNLOAD.DONE
