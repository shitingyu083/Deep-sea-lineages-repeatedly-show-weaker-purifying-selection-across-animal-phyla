#!/bin/bash
# Genome identity audit for every panel assembly (species_codes.tsv accessions):
#   Mash sketches (k=21, s=10000) -> all-vs-all distances -> 35_identity_report.py
#   flags (i) different species with near-identical genomes (D < 0.01; duplicate or
#   mislabelled submissions) and (ii) species whose nearest panel genome is outside its
#   family while a confamilial genome exists (possible misidentification/contamination).
# Usage: setsid nohup bash scripts/34_genome_identity.sh NJ > logs/identity.log 2>&1 < /dev/null &
set -u
R=~/deepsea_conv; export PATH=$R/envs/dc/bin:$PATH; cd $R
NJ=${1:-16}; W=$R/work/identity; mkdir -p $W/sk
export W
sk() {
  a=$1; [ -s $W/sk/$a.msh ] && return 0
  f=$(find data/ncbi/$a -name "*genomic.fna.gz" 2>/dev/null | head -1)
  [ -n "$f" ] || { echo "NOFILE $a"; return 0; }
  nice -n 19 mash sketch -k 21 -s 10000 -o $W/sk/$a $f > /dev/null 2>&1 && echo "OK $a"
}
export -f sk
tail -n +2 meta/species_codes.tsv | cut -f3 | grep ^GC | sort -u | xargs -P $NJ -I{} bash -c 'sk {}'
ls $W/sk/*.msh > $W/sketches.txt
mash paste -l $W/all $W/sketches.txt > /dev/null
mash dist -p $NJ $W/all.msh $W/all.msh > $W/dist.tsv
python3 scripts/35_identity_report.py
echo IDENTITY_DONE
