#!/bin/bash
# Pass 2 of the genome identity audit: waits for pass 1 (IDENTITY_DONE), then
# sketches accessions that had no genome file at pass-1 time (re-downloaded later,
# or 'local' panel species whose proteomes came from figshare/papers -> genome fetched
# here only for the identity check), and recomputes dist/report on all sketches.
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH
W=$R/work/identity; NJ=${1:-16}
until grep -q IDENTITY_DONE logs/identity.log; do sleep 120; done
tail -n +2 meta/species_codes.tsv | cut -f3 | grep ^GC | sort -u | while read a; do
  [ -s $W/sk/$a.msh ] && continue
  f=$(find data/ncbi/$a -name "*genomic.fna.gz" 2>/dev/null | head -1)
  if [ -z "$f" ]; then
    mkdir -p data/ncbi_idonly/$a
    f=$(find data/ncbi_idonly/$a -name "*genomic.fna.gz" 2>/dev/null | head -1)
    if [ -z "$f" ]; then
      ( cd data/ncbi_idonly/$a && datasets download genome accession $a --include genome \
          --filename pkg.zip > dl.log 2>&1 && unzip -oq pkg.zip && rm -f pkg.zip \
          && find . -name "*_genomic.fna" -exec gzip {} \; )
      f=$(find data/ncbi_idonly/$a -name "*genomic.fna.gz" 2>/dev/null | head -1)
    fi
  fi
  [ -n "$f" ] || { echo "NOFILE2 $a"; continue; }
  nice -n 19 mash sketch -k 21 -s 10000 -o $W/sk/$a $f > /dev/null 2>&1 && echo "OK2 $a"
done
ls $W/sk/*.msh > $W/sketches.txt
rm -f $W/all.msh
mash paste -l $W/all $W/sketches.txt > /dev/null
mash dist -p $NJ $W/all.msh $W/all.msh > $W/dist.tsv
envs/dc/bin/python scripts/35_identity_report.py
echo IDENTITY2_DONE $(ls $W/sk/*.msh | wc -l)
