#!/bin/bash
# Actinopteri OrthoFinder on a 77-species subsample (design change 2026-10-20).
# The full 157-species all-vs-all (24,649 DIAMOND searches) was 6.6% done after 17 h, so the
# clade is restricted to the species in work/orthofinder/Actinopteri/sub_v2/selection.tsv
# (1-3 representatives per maximal deep unit of the draft_v2 tree, the closest shallow species
# of each unit's sister clade, one shallow species per remaining family, polar shallow
# notothenioids, two quality backups and the outgroups; all 38 deep units are preserved).
# Searches already finished by the full run are reused after gzip -t; missing ones are run
# with OrthoFinder 2.5.5's own DIAMOND command line; then OrthoFinder -b (removed species
# commented out in SpeciesIDs.txt, the documented way to drop species) builds orthogroups,
# gene trees and HOGs.
# Usage: bash scripts/30_fish_sub_of.sh NJ   (xargs: kill -USR1/-USR2 its PID to add/remove a slot)
R=~/deepsea_conv; cd $R; export PATH=$R/envs/dc/bin:$PATH
echo 10 > /proc/self/autogroup 2>/dev/null
ulimit -n 65536
NJ=${1:-48}
W=$R/work/orthofinder/Actinopteri/of/Results_Oct19/WorkingDirectory
D=$R/work/orthofinder/Actinopteri/sub_v2
IDS=$(grep -v '^#' $D/SpeciesIDs.sub.txt | cut -d: -f1)
echo "$(date +%F_%T) start; species $(echo $IDS | wc -w); NJ $NJ"
nbad=0
for i in $IDS; do for j in $IDS; do f=$W/Blast${i}_${j}.txt.gz
  if [ -f $f ] && ! gzip -t $f 2>/dev/null; then echo "corrupt $f"; rm -f $f; nbad=$((nbad+1)); fi
done; done
rm -f $W/Blast*_*.sub_tmp.txt $W/Blast*_*.sub_tmp.txt.gz
awk -F'\t' 'NR>1{n[$2]=$8} END{for(i in n) for(j in n) print n[i]*n[j], i, j}' $D/selection.tsv \
  | sort -k1,1nr | while read c i j; do [ -f $W/Blast${i}_${j}.txt.gz ] || echo "$i $j"; done > $D/jobs.txt
echo "$(date +%F_%T) corrupt removed $nbad; searches to run $(wc -l < $D/jobs.txt)"
one() {
  i=$1; j=$2; t=$W/Blast${i}_${j}.sub_tmp.txt
  diamond blastp --ignore-warnings -d $W/diamondDBSpecies$j -q $W/Species$i.fa -o $t \
    --more-sensitive -p 1 --quiet -e 0.001 --compress 1
  if [ $? = 0 ] && [ -f $t.gz ]; then mv $t.gz $W/Blast${i}_${j}.txt.gz; echo "$i $j" >> $D/done.txt
  else echo "FAIL $i $j"; rm -f $t $t.gz; fi
}
export -f one; export W D
xargs -P $NJ -L 1 bash -c 'one "$0" "$1"' < $D/jobs.txt
miss=0
for i in $IDS; do for j in $IDS; do [ -f $W/Blast${i}_${j}.txt.gz ] || miss=$((miss+1)); done; done
echo "$(date +%F_%T) searches finished; missing $miss"
if [ $miss -gt 0 ]; then touch $D/SUB.FAIL; exit 1; fi
[ -f $W/SpeciesIDs.full157.txt ] || cp $W/SpeciesIDs.txt $W/SpeciesIDs.full157.txt
cp $D/SpeciesIDs.sub.txt $W/SpeciesIDs.txt
if orthofinder -b $W -t 48 -a 8 -X -n sub77 > $D/of_b.log 2>&1; then
  touch $D/SUB.DONE; echo "$(date +%F_%T) orthofinder -b done"
else
  touch $D/SUB.FAIL; echo "$(date +%F_%T) orthofinder -b FAILED"
fi
