#!/bin/bash
# Corrected DRAFT species tree: the 19f wASTRAL topology with the documented fixes of meta/topo_fix.tsv
# (19g_topo_fix.py), FastTree LG+Gamma ML branch lengths on the concatenated final alignments of the 19f run
# (fixed topology), rooted on the choanoflagellates, gCF, and the quartet-score cost of the fixes
# (wastral -C on the same gene trees for the fixed and the unconstrained topology). Resumable via BL.DONE.
# Outputs in W: topo.tre / topo_unrooted.tre, fix.log.tsv, concat.faa, concat_partition.txt, genetrees.tre,
#   bl.treefile, species_tree_og.nwk, species_tree.nwk, gcf.cf.tree, species_tree_gcf{_og,}.nwk,
#   qc.monophyly.tsv, score_{fixed,unconstrained}.tre + score.tsv                               BL.DONE
# Usage: bash scripts/19g_fixed_bl.sh SRC W [NT]      (SRC = 19f work dir with WASTRAL.DONE)
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
SRC=$(readlink -f $1); mkdir -p $2/log; W=$(readlink -f $2); NT=${3:-16}
PY=$PWD/envs/dc/bin/python
TOOLS="$PY scripts/19c_tools.py"; GT="$PY scripts/19g_tools.py"
log() { echo "$(date '+%F %T') [19g] $*"; }
[ -e $W/BL.DONE ] && exit 0
[ -e $SRC/WASTRAL.DONE ] || { log "no WASTRAL.DONE in $SRC"; exit 1; }
DROP=$($GT dropped meta/topo_fix.tsv)

if [ ! -s $W/topo_unrooted.tre ]; then
  $PY scripts/19g_topo_fix.py $SRC/wastral_topo_unrooted.tre $W/fix Monbre,Salros > $W/log/fix.txt 2>&1 \
    || { log "topology fix failed"; cat $W/log/fix.txt; exit 1; }
  mv $W/fix.tre $W/topo.tre; mv $W/fix_unrooted.tre $W/topo_unrooted.tre
  log "fixed topology: $(tail -1 $W/log/fix.txt)"
fi
$PY scripts/19g_topo_qc.py $W/topo.tre $W/qc Monbre,Salros > $W/log/qc.txt 2>&1; log "QC: $(tail -1 $W/log/qc.txt)"
[ -s $W/concat.faa ] || { $GT dropfa $SRC/concat.faa $W/concat.faa.tmp $DROP && mv $W/concat.faa.tmp $W/concat.faa; }
cp -p $SRC/concat_partition.txt $W/concat_partition.txt
[ -s $W/genetrees.tre ] || { $GT prunegt $SRC/genetrees_shrunk.tre $W/genetrees.tmp $DROP && mv $W/genetrees.tmp $W/genetrees.tre; }
[ -s $W/unconstrained.tre ] || $GT prunetopo $SRC/wastral_topo.tre $W/unconstrained.tre $DROP

# quartet-score cost of the fixes (runs alongside the branch lengths)
if [ ! -s $W/score.tsv ]; then
  ( for k in fixed unconstrained; do
      c=$W/topo.tre; [ $k = unconstrained ] && c=$W/unconstrained.tre
      wastral -C -L --mode 1 -u 2 -t 4 --seed 20261019 -c $c -i $W/genetrees.tre -o $W/score_$k.tre \
        > $W/log/score_$k.log 2>&1
    done
    sf=$(grep -E "^Score" $W/log/score_fixed.log | tail -1 | awk '{print $2}')
    su=$(grep -E "^Score" $W/log/score_unconstrained.log | tail -1 | awk '{print $2}')
    printf "topology\tscore\nfixed\t%s\nunconstrained\t%s\n" "$sf" "$su" > $W/score.tsv ) &
  SCPID=$!
fi

if [ ! -s $W/bl.treefile ]; then
  log "FastTreeMP ML branch lengths on the fixed topology ($NT threads, memguard 64G)"
  OMP_NUM_THREADS=$NT bash scripts/memguard.sh 64 FastTreeMP -lg -gamma -nome -mllen -nosupport \
      -intree $W/topo_unrooted.tre $W/concat.faa > $W/bl.tmp 2> $W/log/bl_fasttree.log \
    && [ -s $W/bl.tmp ] || { log "FastTree branch lengths failed"; exit 1; }
  $PY scripts/19f_tools.py samesplits $W/bl.tmp $W/topo.tre || { log "FastTree changed the topology"; exit 1; }
  mv $W/bl.tmp $W/bl.treefile
  echo "FastTree 2.2.0 (OpenMP) -lg -gamma -nome -mllen -nosupport on concat.faa (19f loci minus $DROP), fixed topology (meta/topo_fix.tsv)" > $W/bl_method.txt
  log "branch lengths done"
fi
$TOOLS root $W/bl.treefile $W/species_tree_og.nwk $W/species_tree.nwk || { log "rooting failed"; exit 1; }
iqtree3 -t $W/bl.treefile --gcf $W/genetrees.tre -T 4 --prefix $W/gcf -redo > $W/log/gcf.log 2>&1 || log "gcf failed (non-fatal)"
[ -s $W/gcf.cf.tree ] && $TOOLS root $W/gcf.cf.tree $W/species_tree_gcf_og.nwk $W/species_tree_gcf.nwk
[ -n "${SCPID:-}" ] && wait $SCPID
log "quartet scores: $(tr '\n' ' ' < $W/score.tsv)"
touch $W/BL.DONE
log "BL.DONE"
