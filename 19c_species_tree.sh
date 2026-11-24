#!/bin/bash
# Global species tree v2 (draft = QC-pass species at launch; final = rerun on the completed panel).
#   loci   W/loci from 18c_busco_loci_v2.py (metazoa_odb12 single-copy BUSCOs, >= 50 % ingroup occupancy,
#          choanoflagellate outgroups Monbre + Salros)
#   1 per locus (19c_locus.sh): MAFFT L-INS-i, trimAl -automated1, drop sequences < 50 % of trimmed
#     columns, IQ-TREE ModelFinder (LG/WAG/JTT/Q.pfam x FU/F x G4/R2-R6) + aBayes gene tree  GENETREES.DONE
#   2-6 19c_post.sh: TreeShrink, wASTRAL, concatenated branch lengths on the wASTRAL topology
#     (DRAFT_BL.DONE), concatenated ML + UFBoot, gCF/sCF, rooting, comparisons (TREE.DONE)
# Usage: setsid nohup bash scripts/19c_species_tree.sh W [NJ=64] [TPJ=3] > logs/tree_v2.log 2>&1 &
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
W=$(readlink -f $1); NJ=${2:-64}; TPJ=${3:-3}
mkdir -p $W/aln $W/trim $W/filt $W/gt $W/skip $W/log
log() { echo "$(date '+%F %T') $*"; }

# 1 gene trees; largest loci first so the longest jobs do not form the tail
if [ ! -e $W/GENETREES.DONE ]; then
  for pass in 1 2; do
    log "gene trees pass $pass ($NJ jobs x $TPJ threads)"
    ls -S $W/loci | sed 's/\.faa$//' | xargs -P $NJ -I{} nice -n 10 bash scripts/19c_locus.sh $W {} $TPJ
    : > $W/failed_loci.txt
    for f in $W/loci/*.faa; do
      b=$(basename $f .faa); [ -s $W/gt/$b.treefile ] || [ -e $W/skip/$b ] || echo $b >> $W/failed_loci.txt
    done
    nf=$(wc -l < $W/failed_loci.txt)
    log "pass $pass: $(ls $W/gt | grep -c '\.treefile$') trees, $(ls $W/skip | wc -l) skipped, $nf failed"
    [ $nf = 0 ] && break
  done
  [ $nf = 0 ] || log "continuing without the $nf loci in failed_loci.txt"
  touch $W/GENETREES.DONE
fi
exec bash scripts/19c_post.sh $W
