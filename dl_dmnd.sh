#!/bin/bash
# Resumable download of eggNOG 5.0.2 DIAMOND db; restarts on stall, verifies size + gzip.
cd ~/deepsea_conv/db/eggnog
U=http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
F=eggnog_proteins.dmnd.gz; N=5208806170
for i in $(seq 1 200); do
  [ "$(stat -c %s $F 2>/dev/null)" = "$N" ] && break
  curl -s -C - --speed-limit 20000 --speed-time 120 --connect-timeout 60 -o $F $U
  sleep 10
done
[ "$(stat -c %s $F)" = "$N" ] || { echo "size mismatch"; exit 1; }
gzip -t $F && gunzip $F && touch dmnd.done && echo DMND_OK
