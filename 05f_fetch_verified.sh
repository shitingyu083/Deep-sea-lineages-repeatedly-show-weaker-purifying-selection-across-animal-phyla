#!/bin/bash
# Verified NCBI datasets download for accessions that never entered the
# acc_list*.tsv queues (2026-10-19: Ampamp GCF_019059575.1, Saccar GCA_916048095.2,
# both in panel_v3 / species_codes). Each package is downloaded into
# data/ncbi/.tmp_<acc>, the zip is tested, md5sum.txt is checked on the unzipped
# files, .fna/.faa/.gff are gzipped, and only then is the directory moved to
# data/ncbi/<acc> with .done (same layout as 05d). 3 tries; FAIL lines as 05d.
# Usage: bash scripts/05f_fetch_verified.sh LIST.tsv     (acc<TAB>include)
set -u
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH
LIST=$1
while IFS=$'\t' read -r acc inc; do
  [ -z "$acc" ] && continue
  d=data/ncbi/$acc
  [ -f $d/.done ] && { echo "$(date +%F_%T) $acc already done"; continue; }
  [ -e $d ] && { echo "$(date +%F_%T) $acc: $d exists without .done, left alone"; continue; }
  t=data/ncbi/.tmp_$acc; ok=0
  for try in 1 2 3; do
    rm -rf $t; mkdir -p $t
    if datasets download genome accession $acc --include $inc --filename $t/pkg.zip \
         --no-progressbar > $t/dl.log 2>&1 \
       && (cd $t && unzip -tq pkg.zip > /dev/null && unzip -oq pkg.zip && rm pkg.zip) \
       && (cd $t && md5sum --quiet -c md5sum.txt); then
      ok=1; break
    fi
    echo "$(date +%F_%T) $acc try $try failed"; sleep $((try * 30))
  done
  if [ $ok = 1 ]; then
    find $t -type f \( -name '*.fna' -o -name '*.faa' -o -name '*.gff' \) -exec gzip {} +
    mv $t $d && touch $d/.done && echo "$(date +%F_%T) $acc OK ($inc)"
  else
    echo "FAIL $acc" >> data/ncbi/failed.txt; echo "$(date +%F_%T) $acc FAILED"
  fi
done < $LIST
echo "$(date +%F_%T) fetch list done"
