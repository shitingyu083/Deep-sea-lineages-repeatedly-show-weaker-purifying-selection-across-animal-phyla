#!/bin/bash
# Stages 2-6 of the global species tree (after the gene trees of 19c_species_tree.sh); resumable via *.DONE.
#   1b wait for extra gene-tree workers (19c_worker.sh) and running claimed loci, then one retry of
#     unfinished loci (none after W/CUT, the time-boxed draft cut-off)                       W/failed_loci.txt
#   2 TreeShrink per-species q = 0.05 on the aBayes gene trees; removals applied to the original trees
#     and alignments -> W/final_aln, W/genetrees_shrunk.tre, W/partition.nex              SHRINK.DONE
#   3 wASTRAL hybrid weighting (aBayes support + branch length; -B aBayes support range, -R more
#     search rounds); the tree re-scored with full quartet annotation (-u 2) and frequencies (-u 3)
#                                                                                          WASTRAL.DONE
#   4 partitioned IQ-TREE branch lengths on the wASTRAL topology (each locus with the model ModelFinder
#     chose for its gene tree), rooted on the choanoflagellates                            DRAFT_BL.DONE
#   5 partitioned concatenated IQ-TREE ML search, UFBoot 1000 resampling genes then sites   ML.DONE
#   6 gCF (gene trees) and likelihood sCF (100 quartets per branch) on both trees, RF comparisons with
#     each other and with the per-clade ASTRAL trees used by the selection analyses        TREE.DONE
# W/DRAFT (2026-10-19): draft run, superseded by the final-panel tree -> no stage 5, gCF only in stage 6.
# Heavy stages run under setsid (own session = own CPU share under autogroup scheduling).
# Usage: bash scripts/19c_post.sh W        (exec'd by 19c_species_tree.sh)
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
P=~/anaconda3/envs/py310env/bin/python                 # TreeShrink lives here
W=$(readlink -f $1)
TOOLS="$PWD/envs/dc/bin/python scripts/19c_tools.py"   # ete3 3.1.3 lives here
log() { echo "$(date '+%F %T') $*"; }
DRAFT=0; [ -e $W/DRAFT ] && DRAFT=1

# 1b gene trees may still be running in worker sessions: wait for the workers and for every unfinished locus
#    whose claim lock is held, then retry what is left once (not after a draft cut-off)
unfinished() { for f in $W/loci/*.faa; do b=$(basename $f .faa)
  [ -s $W/gt/$b.treefile ] || [ -e $W/skip/$b ] || echo $b; done; }
held() { unfinished | while read b; do flock -n $W/claim/$b.lock true 2>/dev/null || echo $b; done; }
busy() { pgrep -f "19c_worker.sh $W " > /dev/null || [ -n "$(held)" ]; }
if [ ! -e $W/SHRINK.DONE ]; then
  mkdir -p $W/claim
  while busy; do sleep 120; done
  if [ ! -e $W/CUT ] && [ -n "$(unfinished)" ]; then
    log "retrying $(unfinished | wc -l) unfinished loci"
    unfinished | xargs -r -P 48 -I{} nice -n 10 bash scripts/19c_locus.sh $W {} 3
    while busy; do sleep 120; done
  fi
  unfinished > $W/failed_loci.txt
  log "gene trees: $(ls $W/gt | grep -c '\.treefile$') trees, $(ls $W/skip | wc -l) skipped, $(wc -l < $W/failed_loci.txt) unfinished$([ -e $W/CUT ] && echo ' (draft cut-off)')"
fi

# 2 TreeShrink
if [ ! -e $W/SHRINK.DONE ]; then
  log "TreeShrink"
  ls $W/gt | grep '\.treefile$' | sed 's/\.treefile$//' > $W/genetrees.list
  while read b; do cat $W/gt/$b.treefile; done < $W/genetrees.list > $W/genetrees.tre
  rm -rf $W/treeshrink
  nice -n 10 $P ~/anaconda3/envs/py310env/bin/run_treeshrink.py -t $W/genetrees.tre -o $W/treeshrink \
    -m per-species -q 0.05 > $W/log/treeshrink.log 2>&1 || { log "treeshrink failed"; exit 1; }
  $TOOLS shrink $W || { log "shrink failed"; exit 1; }
  touch $W/SHRINK.DONE
fi

# 3 wASTRAL
if [ ! -e $W/WASTRAL.DONE ]; then
  log "wASTRAL"
  setsid nice -n 10 wastral -B --mode 1 -R -u 1 -t 32 --seed 20261019 -i $W/genetrees_shrunk.tre -o $W/wastral.tre \
    > $W/log/wastral.log 2>&1 || { log "wastral failed"; exit 1; }
  $TOOLS topo $W/wastral.tre $W/wastral_topo.tre
  mkdir -p $W/wastral_score
  for u in 2 3; do
    ( cd $W/wastral_score && setsid nice -n 10 wastral -B --mode 1 -t 32 -u $u -C -c $W/wastral_topo.tre \
        -i $W/genetrees_shrunk.tre -o $W/wastral_score/wastral_u$u.tre > $W/log/wastral_u$u.log 2>&1 ) \
      || log "wastral -u $u scoring failed (non-fatal)"
  done
  touch $W/WASTRAL.DONE
fi

# memory caps from what is available now: other users' jobs share the machine, never invite the OOM killer
avail=$(free -g | awk '/^Mem:/{print $7}')
MBL=$(( avail * 25 / 100 )); [ $MBL -gt 200 ] && MBL=200
MML=$(( avail * 35 / 100 )); [ $MML -gt 300 ] && MML=300
log "available ${avail}G -> bl_astral/cf -mem ${MBL}G, ml -mem ${MML}G"

# 4 + 5 in parallel: branch lengths on the wASTRAL topology (unblocks the downstream chain),
#   independent concatenated ML search (final runs only)
if [ ! -e $W/DRAFT_BL.DONE ]; then
  ( [ -s $W/bl_astral.treefile ] || setsid nice -n 10 iqtree3 -p $W/partition.nex -te $W/wastral_topo.tre -T 64 \
        -mem ${MBL}G -seed 20261019 --prefix $W/bl_astral > $W/log/bl_astral.log 2>&1
    [ -s $W/bl_astral.treefile ] \
      && $TOOLS root $W/bl_astral.treefile $W/species_tree_astral_og.nwk $W/species_tree_astral.nwk \
      && touch $W/DRAFT_BL.DONE
    log "bl_astral finished; DRAFT_BL.DONE $([ -e $W/DRAFT_BL.DONE ] && echo written || echo NOT written)" ) &
fi
if [ ! -e $W/ML.DONE ] && [ $DRAFT = 0 ]; then
  ( [ -s $W/ml.treefile ] || setsid nice -n 10 iqtree3 -p $W/partition.nex -B 1000 --sampling GENESITE -T 96 \
        -mem ${MML}G -seed 20261019 --prefix $W/ml > $W/log/ml.log 2>&1
    [ -s $W/ml.treefile ] \
      && $TOOLS root $W/ml.treefile $W/species_tree_ml_og.nwk $W/species_tree_ml.nwk \
      && touch $W/ML.DONE
    log "ml finished; ML.DONE $([ -e $W/ML.DONE ] && echo written || echo NOT written)" ) &
fi
wait
[ -e $W/DRAFT_BL.DONE ] && { [ -e $W/ML.DONE ] || [ $DRAFT = 1 ]; } || { log "concatenated runs incomplete"; exit 1; }

# 6 concordance factors (non-fatal): gCF needs only the gene trees; sCFL re-fits the partition model on
#   the fixed tree and samples 100 quartets per branch (final runs only)
cf() {
  local t=$1 tree=$2 bm=$3
  [ -s $bm ] || bm=$W/partition.nex
  [ -s $W/gcf_$t.cf.tree ] || iqtree3 -t $tree --gcf $W/genetrees_shrunk.tre -T 4 --prefix $W/gcf_$t -redo \
      > $W/log/gcf_$t.log 2>&1 || log "gcf_$t failed (non-fatal)"
  [ $DRAFT = 1 ] && return 0
  [ -s $W/scfl_$t.cf.tree ] || setsid nice -n 10 iqtree3 -te $tree -p $bm --scfl 100 -T 48 -mem ${MBL}G -seed 20261019 \
      --prefix $W/scfl_$t -redo > $W/log/scfl_$t.log 2>&1 || log "scfl_$t failed (non-fatal)"
}
cf astral $W/bl_astral.treefile $W/bl_astral.best_model.nex &
[ $DRAFT = 1 ] || cf ml $W/ml.treefile $W/ml.best_model.nex &
wait
for t in astral ml; do
  for k in gcf scfl; do
    [ -s $W/${k}_$t.cf.tree ] && $TOOLS root $W/${k}_$t.cf.tree $W/species_tree_${t}_${k}_og.nwk $W/species_tree_${t}_$k.nwk
  done
done

# comparisons: ML vs wASTRAL, and both against the per-clade ASTRAL trees used by HyPhy/csubst
[ -s $W/species_tree_ml.nwk ] && $TOOLS compare $W $W/species_tree_ml.nwk $W/species_tree_astral.nwk ml_vs_astral
for c in work/tree_clade/*/; do
  cn=$(basename $c)
  [ -s $c/astral.tre ] || continue
  $TOOLS compare $W $W/species_tree_astral.nwk $c/astral.tre astral_vs_clade_$cn
  [ -s $W/species_tree_ml.nwk ] && $TOOLS compare $W $W/species_tree_ml.nwk $c/astral.tre ml_vs_clade_$cn
done
touch $W/TREE.DONE
log "TREE.DONE"
