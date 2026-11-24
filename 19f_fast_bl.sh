#!/bin/bash
# Stage 4 of the fast DRAFT species tree (exec'd by 19f_fast_post.sh); resumable via DRAFT_BL.DONE.
# Branch lengths on the fixed wASTRAL topology: FastTree LG+Gamma ML lengths (-nome -mllen) on the concatenated
# final alignments. Draft scaffold only: the final 19c tree gets partitioned IQ-TREE lengths (iqtree3 rejects
# -mem with partition models, and a 396-taxon partitioned run would take hours). The topology is checked
# unchanged (RF = 0), rooted on the choanoflagellates; gCF from the gene trees; RF against the per-clade trees.
# Outputs: concat.faa, concat_partition.txt, concat1.nex (one LG+G4 charset, for dating with fixed lengths),
#   bl_astral.treefile (+ bl_method.txt), species_tree_astral{_og,}.nwk, gcf_astral.cf.tree,
#   species_tree_astral_gcf{_og,}.nwk                                                              DRAFT_BL.DONE
# Usage: bash scripts/19f_fast_bl.sh W
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
W=$(readlink -f $1)
TOOLS="$PWD/envs/dc/bin/python scripts/19c_tools.py"
FTOOLS="$PWD/envs/dc/bin/python scripts/19f_tools.py"
log() { echo "$(date '+%F %T') [bl] $*"; }
[ -e $W/DRAFT_BL.DONE ] && exit 0

if [ ! -s $W/bl_astral.treefile ]; then
  log "concatenating the final alignments"
  $FTOOLS concat $W || { log "concat failed"; exit 1; }
  log "FastTree ML branch lengths on the fixed wASTRAL topology (memguard 64G)"
  OMP_NUM_THREADS=1 bash scripts/memguard.sh 64 FastTree -lg -gamma -nome -mllen -nosupport \
      -intree $W/wastral_topo_unrooted.tre $W/concat.faa > $W/bl_astral.tmp 2> $W/log/bl_fasttree.log \
    && [ -s $W/bl_astral.tmp ] || { log "FastTree branch lengths failed"; exit 1; }
  $FTOOLS samesplits $W/bl_astral.tmp $W/wastral_topo.tre || { log "FastTree changed the topology"; exit 1; }
  mv $W/bl_astral.tmp $W/bl_astral.treefile
  echo "FastTree 2.2.0 -lg -gamma -nome -mllen -nosupport on concat.faa, fixed wASTRAL topology" > $W/bl_method.txt
  log "branch lengths done"
fi

$TOOLS root $W/bl_astral.treefile $W/species_tree_astral_og.nwk $W/species_tree_astral.nwk \
  || { log "rooting failed"; exit 1; }
iqtree3 -t $W/bl_astral.treefile --gcf $W/genetrees_shrunk.tre -T 4 --prefix $W/gcf_astral -redo \
  > $W/log/gcf_astral.log 2>&1 || log "gcf failed (non-fatal)"
[ -s $W/gcf_astral.cf.tree ] \
  && $TOOLS root $W/gcf_astral.cf.tree $W/species_tree_astral_gcf_og.nwk $W/species_tree_astral_gcf.nwk
for c in work/tree_clade/*/; do
  cn=$(basename $c); [ -s $c/astral.tre ] || continue
  $TOOLS compare $W $W/species_tree_astral.nwk $c/astral.tre astral_vs_clade_$cn
done
touch $W/DRAFT_BL.DONE
log "DRAFT_BL.DONE"
