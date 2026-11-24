#!/bin/bash
# FINAL species tree, stages 2-7 (2026-10-20): exec'd by 19h_final_tree.sh after GENETREES.DONE (keeps its driver
# lock, fd 8), or run by hand to resume (takes the lock). Every stage is resumable (markers, IQ-TREE checkpoints).
#   2 TreeShrink (per-species, q 0.05; all-genes if per-species fails) on the gene trees of the finished IQ-TREE
#     runs (unfinished ones are run again first) -> W/genetrees_shrunk.tre, final_aln/, partition.nex (per-locus
#     ModelFinder models), locus_stats.tsv, locus_methods.tsv, treeshrink_method.txt                  W/SHRINK.DONE
#   3 wASTRAL (-B aBayes weights, hybrid weighting, -R) -> W/wastral_topo.tre; quartet details (-u 2, -u 3);
#     concatenation (19f_tools.py concat)                                                              W/WASTRAL.DONE
#   4 fixed topology: 19g_topo_fix.py SPEC on the unconstrained tree, then 19g_topo_qc.py; a non-monophyletic
#     group above order/family level that is not in ACC -> F/REVIEW_NEEDED and a wait until SPEC changes or
#     F/TOPO.OK exists (a failed fix waits for a SPEC change only)                                      F/TOPO.DONE
#   5 concatenated ML (branch ml, from stage 2 on): IQ-TREE -p partition.nex, UFBoot 1000 (sites resampled within
#     partitions) + SH-aLRT 1000 -> W/ml.treefile, species_tree_ml{,_og,_sup,_sup_og}.nwk, ml_topo.tre   W/ML.DONE
#   6 branch lengths on the fixed topology (IQ-TREE -p partition.nex -te), rooting, gCF; quartet scores of the
#     fixed and the unconstrained topology (wastral -C, the cost of the fixes) -> F/score.tsv            F/BL.DONE
#   7 downstream chains: 36_rer_trees.sh (branch rer36, from stage 4 on), then 38c_chain.sh on the fixed tree
#     (dating, ASR, PGLS, RERconverge; all label sets, NPERM_MAIN permulations)                       F/CHAIN.DONE
#     branch unc: the unconstrained wASTRAL tree in F_unc (branch lengths, 36, 38c with the main labels and
#       NPERM_SENS; SKIPPED if it equals the fixed topology)                                      F_unc/CHAIN.DONE
#     branch cf: sCFL, gCF of the unconstrained and ML trees, quartet score of the ML topology, AU test (fixed,
#       unconstrained, ML), RF comparisons incl. the clade trees of work/tree_clade                     F/CF.DONE
#     branch mlsens: the ML tree in F_ml (36, 38c as unc; SKIPPED if it equals the fixed or the unconstrained
#       topology)                                                                                  F_ml/CHAIN.DONE
#   then, once every branch has finished                                                              W/TREE.DONE
# Branches run in their own sessions, one per W/locks/MODE.lock (held by all their processes); restart one with
# bash scripts/19h_post.sh W F RV MODE. Big IQ-TREE runs share SLOTS slots, start only with >= GATE GB available
# and run under memguard.sh GUARD (killed below GUARD GB MemAvailable, restarted later from the checkpoint).
# Usage: bash scripts/19h_post.sh W [F=work/final_v3] [RV=work/rer_v3] [MODE=main|ml|rer36|unc|cf|mlsens]
#   env SPEC (meta/topo_fix_v3.tsv), NPERM_MAIN (20000), NPERM_SENS (10000), CHAIN (scripts/38c_chain.sh; smoke tests)
set -u
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
P=~/anaconda3/envs/py310env/bin/python; PY=$PWD/envs/dc/bin/python
TOOLS="$PY scripts/19c_tools.py"; FT="$PY scripts/19f_tools.py"
mkdir -p ${2:-work/final_v3}/log ${3:-work/rer_v3}
W=$(readlink -f $1); F=$(readlink -f ${2:-work/final_v3}); RV=$(readlink -f ${3:-work/rer_v3}); MODE=${4:-main}
U=${F}_unc; RU=${RV}_unc; M=${F}_ml; RM=${RV}_ml
SPEC=${SPEC:-meta/topo_fix_v3.tsv}; ACC=meta/topo_qc_accepted_v3.tsv; OG=Monbre,Salros
SLOTS=4; GATE=200; GUARD=100; NPERM_MAIN=${NPERM_MAIN:-20000}; NPERM_SENS=${NPERM_SENS:-10000}
CHAIN=${CHAIN:-scripts/38c_chain.sh}; GTS=$W/genetrees_shrunk.tre
NP=$(nproc); thr() { [ $1 -le $NP ] && echo $1 || echo $NP; }   # IQ-TREE refuses more threads than (allowed) cores
mkdir -p $W/locks $W/log
log() { echo "$(date '+%F %T') [19h post $MODE] $*"; }
memavail() { awk '/^MemAvailable/ {print int($2 / 1048576)}' /proc/meminfo; }
iqdone() { grep -q "^Date and Time:" $1.log 2>/dev/null && [ -s $1.${2:-treefile} ]; }   # PREFIX [OUTPUT_SUFFIX]
running() { ! flock -n $W/locks/$1.lock true 2>/dev/null; }
specmd5() { md5sum < $SPEC | cut -c1-32; }
iqver() { iqtree3 --version 2>&1 | head -1 | awk '{for (i = 1; i < NF; i++) if ($i == "version") {print $(i + 1); exit}}'; }
qscore() { grep -E "^Score" $1 2>/dev/null | tail -1 | awk '{print $2}'; }
# mgrun LOG CMD...: run CMD under memguard in a free slot (its lock is inherited, so held while CMD lives); a
# memguard stop (rc 99) is retried later (IQ-TREE resumes from its checkpoint); returns the rc of CMD
mgrun() { local lg=$1 s rc; shift
  while :; do
    for s in $(seq 1 $SLOTS); do
      [ $(memavail) -ge $GATE ] || break
      flock -n -E 75 $W/locks/big$s.lock bash scripts/memguard.sh $GUARD "$@" >> $lg 2>&1; rc=$?
      [ $rc = 75 ] && continue
      [ $rc = 99 ] && { log "memguard stopped $(basename $lg .run.log) (slot $s); resuming later"
                        sleep $((300 + RANDOM % 900)); continue 2; }
      return $rc
    done
    sleep 60
  done; }
# iqrun PREFIX OUTPUT_SUFFIX ARGS...: IQ-TREE run, resumed from its checkpoint, finally retried with -redo
iqrun() { local p=$1 e=$2 rc a; shift 2
  for a in "" "" -redo; do
    iqdone $p $e && return 0
    [ -n "$a" ] && log "$(basename $p): restarting with -redo"
    mgrun $p.run.log iqtree3 "$@" --prefix $p $a; rc=$?
    iqdone $p $e && return 0
    log "$(basename $p): IQ-TREE rc $rc, see $p.log"
  done
  log "$(basename $p) failed"; return 1; }
# branch MODE [DONE_MARKER...]: start MODE in its own session unless it runs or one of its markers exists
branch() { local m=$1 f; shift
  for f in "$@"; do [ -e $f ] && return 0; done
  running $m && return 0
  (setsid nohup bash scripts/19h_post.sh $W $F $RV $m >> $W/log/post_$m.log 2>&1 < /dev/null 8>&- 9>&- &)
  log "branch $m started ($W/log/post_$m.log)"; sleep 5; }
# need WHAT BRANCH TEST...: wait until TEST succeeds; fail if BRANCH (its producer) is no longer running
need() { local what=$1 prod=$2; shift 2
  until "$@"; do
    running $prod || { sleep 10; "$@" && return 0; log "$what missing and branch $prod is not running"; return 1; }
    sleep 300
  done; }
# the mode lock; -w 5 rides out the instant a running() test holds it
lock_mode() { exec 9> $W/locks/$MODE.lock; flock -w 5 9 || { log "already running"; exit 0; }; }
# rootcf CFTREE NAME: rooted copies F/species_tree_NAME{_og,}.nwk (rc 3, a non-monophyletic ingroup, is reported)
rootcf() { local rc; $TOOLS root $1 $F/species_tree_$2_og.nwk $F/species_tree_$2.nwk >> $F/log/root_cf.txt 2>&1; rc=$?
  [ $rc = 3 ] && log "$2: ingroup not monophyletic once rooted"; [ $rc = 0 ] || [ $rc = 3 ]; }
# gcf TREE NAME: gene concordance factors of TREE (the shrunk gene trees) -> F/NAME.cf.{tree,branch,stat}
gcf() { iqtree3 -t $1 --gcf $GTS -T 4 --prefix $F/$2 -redo > $F/log/$2.log 2>&1 && [ -s $F/$2.cf.tree ] && rootcf $F/$2.cf.tree $2; }
TS=~/anaconda3/envs/py310env/bin/run_treeshrink.py

case $MODE in
ml)        # ---- concatenated ML tree (from stage 2 on)
  lock_mode; [ -e $W/ML.DONE ] && exit 0; [ -e $W/SHRINK.DONE ] || { log "no SHRINK.DONE"; exit 1; }
  # UFBoot resamples sites within partitions: with taxa missing from partitions, IQ-TREE 3.1.4 --sampling GENESITE
  # (and GENE) crashes while it assigns the supports, after the whole search (4-locus toy: work/si/iqtree_genesite)
  log "concatenated ML: IQ-TREE -p partition.nex, UFBoot 1000 (sites resampled within partitions) + SH-aLRT 1000"
  iqrun $W/ml treefile -p $W/partition.nex -B 1000 -alrt 1000 -T $(thr 64) -seed 20261019 || exit 1
  echo "IQ-TREE $(iqver) -p partition.nex (per-locus ModelFinder models, edge-proportional), UFBoot 1000 (sites resampled within partitions), SH-aLRT 1000, seed 20261019" > $W/ml_method.txt
  sed -E 's/\)[^:,();]+/)/g' $W/ml.treefile > $W/ml_nolab.treefile
  $TOOLS root $W/ml_nolab.treefile $W/species_tree_ml_og.nwk $W/species_tree_ml.nwk > $W/log/ml_root.txt 2>&1; rc=$?
  if [ $rc = 3 ]; then touch $W/ML_INGROUP.WARN; log "WARNING: ingroup not monophyletic in the ML tree ($(tail -n 1 $W/log/ml_root.txt))"
  elif [ $rc != 0 ]; then log "ML tree: rooting failed"; cat $W/log/ml_root.txt; exit 1; fi
  $TOOLS root $W/ml.treefile $W/species_tree_ml_sup_og.nwk $W/species_tree_ml_sup.nwk >> $W/log/ml_root.txt 2>&1
  $TOOLS topo $W/ml_nolab.treefile $W/ml_topo.tre
  touch $W/ML.DONE; log "ML.DONE: $(grep -m1 "BEST SCORE FOUND" $W/ml.log)"; exit 0 ;;
rer36)     # ---- per-locus branch lengths on the fixed topology for RERconverge (from stage 4 on)
  lock_mode; [ -e $RV/TREES.DONE ] && exit 0; [ -e $F/TOPO.DONE ] || { log "no TOPO.DONE"; exit 1; }
  bash scripts/36_rer_trees.sh $W $F/topo.tre $RV 96 >> $F/log/36.log 2>&1 || { log "36_rer_trees.sh failed ($F/log/36.log)"; exit 1; }
  log "$(tail -n 1 $F/log/36.log)"; exit 0 ;;
unc)       # ---- sensitivity: the unconstrained wASTRAL topology (from stage 4 on)
  lock_mode; mkdir -p $U/log
  { [ -e $U/CHAIN.DONE ] || [ -e $U/SKIPPED ]; } && exit 0
  [ -e $F/TOPO.DONE ] || { log "no TOPO.DONE"; exit 1; }
  if $FT samesplits $W/wastral_topo.tre $F/topo.tre > $U/log/samesplits.txt 2>&1; then
    echo "unconstrained wASTRAL topology = fixed topology ($(cat $U/log/samesplits.txt))" > $U/SKIPPED
    log "SKIPPED: $(cat $U/SKIPPED)"; exit 0
  fi
  log "unconstrained wASTRAL topology ($(cat $U/log/samesplits.txt) against the fixed one)"
  iqrun $U/bl treefile -p $W/partition.nex -te $W/wastral_topo_unrooted.tre -T $(thr 32) -seed 20261019 || exit 1
  $FT samesplits $U/bl.treefile $W/wastral_topo.tre > $U/log/bl_check.txt 2>&1 || { log "IQ-TREE changed the topology"; exit 1; }
  $TOOLS root $U/bl.treefile $U/tree_og.nwk $U/tree.nwk > $U/log/root.txt 2>&1; rc=$?
  if [ $rc = 3 ]; then echo "ingroup not monophyletic in the unconstrained tree ($(tail -n 1 $U/log/root.txt))" > $U/SKIPPED
    log "SKIPPED: $(cat $U/SKIPPED)"; exit 0
  elif [ $rc != 0 ]; then log "rooting failed: $(tail -n 1 $U/log/root.txt)"; exit 1; fi
  bash scripts/36_rer_trees.sh $W $U/tree_og.nwk $RU 48 >> $U/log/36.log 2>&1 || { log "RER gene trees failed ($U/log/36.log)"; exit 1; }
  NPERM=$NPERM_SENS NCORES=32 bash $CHAIN $U $U/tree_og.nwk $W/concat.faa $RU $RU/rer main >> $U/log/38c.log 2>&1 \
    || { log "chain failed ($U/log/38c.log)"; exit 1; }
  log "CHAIN.DONE"; exit 0 ;;
cf)        # ---- concordance factors, quartet score of the ML topology, AU test, RF comparisons (from stage 6 on)
  lock_mode; [ -e $F/CF.DONE ] && exit 0; [ -e $F/BL.DONE ] || { log "no BL.DONE"; exit 1; }; fail=0
  [ -s $F/gcf.cf.tree ] || gcf $F/bl.treefile gcf || { log "gCF (fixed) failed"; fail=$((fail + 1)); }
  if iqrun $F/scfl cf.tree -te $F/bl.treefile -p $W/partition.nex --scfl 100 -T $(thr 48) -seed 20261019; then
    rootcf $F/scfl.cf.tree scfl || fail=$((fail + 1))
  else fail=$((fail + 1)); fi
  ubl() { iqdone $U/bl || [ -e $U/SKIPPED ]; }
  if need "unconstrained branch lengths" unc ubl; then
    if iqdone $U/bl && [ ! -s $F/gcf_unc.cf.tree ]; then gcf $U/bl.treefile gcf_unc || { log "gCF (unconstrained) failed"; fail=$((fail + 1)); }; fi
  else fail=$((fail + 1)); fi
  if need "ML tree" ml test -e $W/ML.DONE; then
    [ -s $F/gcf_ml.cf.tree ] || gcf $W/ml_nolab.treefile gcf_ml || { log "gCF (ML) failed"; fail=$((fail + 1)); }
    if [ ! -s $F/score_ml.tsv ]; then
      wastral -C -B --mode 1 -u 2 -t 8 --seed 20261019 -c $W/ml_topo.tre -i $GTS -o $F/score_ml.tre > $F/log/score_ml.log 2>&1 \
        && [ -n "$(qscore $F/log/score_ml.log)" ] && printf "topology\tscore\nml\t%s\n" "$(qscore $F/log/score_ml.log)" > $F/score_ml.tsv \
        || { log "quartet score (ML) failed"; fail=$((fail + 1)); }
    fi
    if ! iqdone $F/au iqtree; then
      for f in $F/topo_unrooted.tre $W/wastral_topo_unrooted.tre $W/ml_topo.tre; do cat $f; echo; done | sed '/^$/d' > $F/au_trees.tre
      printf "tree\ttopology\tfile\n1\tfixed\t$F/topo_unrooted.tre\n2\tunconstrained wASTRAL\t$W/wastral_topo_unrooted.tre\n3\tconcatenated ML\t$W/ml_topo.tre\n" > $F/au_trees.tsv
    fi
    iqrun $F/au iqtree -p $W/partition.nex -t $W/ml_nolab.treefile -n 0 -z $F/au_trees.tre -zb 10000 -au -T $(thr 32) -seed 20261019 \
      || fail=$((fail + 1))
  else fail=$((fail + 1)); fi
  { $TOOLS compare $F $F/topo.tre $W/wastral_topo.tre fixed_vs_unc
    if [ -e $W/ML.DONE ]; then
      $TOOLS compare $F $F/topo.tre $W/ml_nolab.treefile fixed_vs_ml
      $TOOLS compare $F $W/wastral_topo.tre $W/ml_nolab.treefile unc_vs_ml
    fi
    for c in work/tree_clade/*/astral.tre; do
      [ -s "$c" ] || continue; k=$(basename $(dirname $c))
      $TOOLS compare $F $F/species_tree.nwk $c fixed_vs_clade_$k
      [ -s $W/species_tree_ml.nwk ] && $TOOLS compare $F $W/species_tree_ml.nwk $c ml_vs_clade_$k
    done; } > $F/log/compare.txt 2>&1
  [ $fail = 0 ] || { log "$fail concordance/test step(s) failed; rerun: bash scripts/19h_post.sh $W $F $RV cf"; exit 1; }
  touch $F/CF.DONE; log "CF.DONE"; exit 0 ;;
mlsens)    # ---- sensitivity: the concatenated ML tree (from stage 6 on, once ML.DONE)
  lock_mode; mkdir -p $M/log
  { [ -e $M/CHAIN.DONE ] || [ -e $M/SKIPPED ]; } && exit 0
  need "ML tree" ml test -e $W/ML.DONE || exit 1
  if [ -e $W/ML_INGROUP.WARN ]; then
    echo "ingroup not monophyletic in the ML tree ($W/log/ml_root.txt)" > $M/SKIPPED; log "SKIPPED: $(cat $M/SKIPPED)"; exit 0
  fi
  for x in fixed:$F/topo.tre unconstrained:$W/wastral_topo.tre; do
    if $FT samesplits $W/ml_nolab.treefile ${x#*:} > $M/log/samesplits_${x%%:*}.txt 2>&1; then
      echo "ML topology = ${x%%:*} topology ($(cat $M/log/samesplits_${x%%:*}.txt))" > $M/SKIPPED; log "SKIPPED: $(cat $M/SKIPPED)"; exit 0
    fi
  done
  log "ML topology: $(cat $M/log/samesplits_fixed.txt) against the fixed one"
  bash scripts/36_rer_trees.sh $W $W/species_tree_ml_og.nwk $RM 48 >> $M/log/36.log 2>&1 || { log "RER gene trees failed ($M/log/36.log)"; exit 1; }
  NPERM=$NPERM_SENS NCORES=32 bash $CHAIN $M $W/species_tree_ml_og.nwk $W/concat.faa $RM $RM/rer main >> $M/log/38c.log 2>&1 \
    || { log "chain failed ($M/log/38c.log)"; exit 1; }
  log "CHAIN.DONE"; exit 0 ;;
main) ;;
*) log "unknown mode $MODE"; exit 1 ;;
esac

# ---- main: stages 2-7 in order (holds W/driver.lock: inherited from 19h_final_tree.sh, or taken here)
[ "$(readlink /proc/$$/fd/8 2>/dev/null)" = "$W/driver.lock" ] \
  || { exec 8> $W/driver.lock; flock -n 8 || { log "another driver holds $W/driver.lock"; exit 1; }; }
[ -e $W/TREE.DONE ] && { log "TREE.DONE exists"; exit 0; }
[ -e $W/GENETREES.DONE ] || { log "no GENETREES.DONE in $W"; exit 1; }
log "start: W $W, F $F, RV $RV; SPEC $SPEC ($(specmd5)); NPERM main $NPERM_MAIN, sensitivity $NPERM_SENS; IQ-TREE $(iqver)"
# ---- stage 2: TreeShrink
tm() { local x; x=$(tail -n 1 $W/log/$1 2>/dev/null | awk -v k=$2 'NF >= 2 {print $k}'); echo ${x:-NA}; }   # FILE FIELD of /usr/bin/time "%e %M"
if [ ! -e $W/SHRINK.DONE ]; then
  # IQ-TREE writes intermediate treefiles: a locus whose run did not finish (no "Date and Time:" in its log) is run
  # again (19h_locus.sh resumes it from its checkpoint) and left out if it still does not finish
  ls $W/gt | sed -n 's/\.treefile$//p' | sort | while read b; do iqdone $W/gt/$b || echo $b; done > $W/gt_rerun.txt
  if [ -s $W/gt_rerun.txt ]; then
    log "$(wc -l < $W/gt_rerun.txt) loci with an unfinished IQ-TREE run: running them again"
    awk -F'\t' 'NR == FNR {r[$1]; next} FNR > 1 && ($1 in r) {print $1 "\t" $5}' $W/gt_rerun.txt $W/cost.tsv \
      | xargs -P 48 -L1 bash scripts/19h_locus.sh $W > $W/log/gt_rerun.log 2>&1 8>&-
  fi
  ls $W/gt | sed -n 's/\.treefile$//p' | sort | while read b; do iqdone $W/gt/$b && echo $b; done > $W/genetrees.list
  grep -vxFf $W/genetrees.list $W/gt_rerun.txt > $W/gt_excluded.txt
  [ -s $W/gt_excluded.txt ] && log "left out (IQ-TREE run still unfinished): $(tr '\n' ' ' < $W/gt_excluded.txt)"
  while read b; do head -n 1 $W/gt/$b.treefile; done < $W/genetrees.list > $W/genetrees.tre
  rm -rf $W/treeshrink; tsm=per-species
  log "TreeShrink ($tsm, q 0.05) on $(wc -l < $W/genetrees.list) gene trees"
  if ! $P $TS -t $W/genetrees.tre -o $W/treeshrink -m $tsm -q 0.05 > $W/log/treeshrink.log 2>&1; then
    mv $W/log/treeshrink.log $W/log/treeshrink_per-species.log; rm -rf $W/treeshrink; tsm=all-genes
    log "TreeShrink per-species failed ($(grep -E '^[A-Za-z.]*(Error|Exception)' $W/log/treeshrink_per-species.log | tail -n 1)); retrying with -m $tsm"
    $P $TS -t $W/genetrees.tre -o $W/treeshrink -m $tsm -q 0.05 > $W/log/treeshrink.log 2>&1 || { log "TreeShrink failed"; exit 1; }
  fi
  echo "TreeShrink $(grep -m1 -oiE 'version [0-9.]+' $W/log/treeshrink.log | awk '{print $2}') -m $tsm -q 0.05 on $(wc -l < $W/genetrees.list) gene trees ($(wc -l < $W/gt_rerun.txt) loci rerun after an unfinished IQ-TREE run, $(wc -l < $W/gt_excluded.txt) left out)" > $W/treeshrink_method.txt
  $TOOLS shrink $W > $W/log/shrink.log 2>&1 || { log "shrink failed"; cat $W/log/shrink.log; exit 1; }
  log "$(tail -n 1 $W/log/shrink.log)"
  { printf "locus\tn_taxa\tn_removed\tlength\tmodel\taligner\tiqtree_threads\tmafft_s\tiqtree_s\tiqtree_maxrss_kb\n"
    tail -n +2 $W/locus_stats.tsv | while IFS=$'\t' read b n r L m; do
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $b $n $r $L $m "$(head -n 1 $W/aln/$b.method 2>/dev/null || echo NA)" \
        "$(cat $W/log/$b.iqtree.threads 2>/dev/null || echo NA)" $(tm $b.mafft.time 1) $(tm $b.iqtree.time 1) $(tm $b.iqtree.time 2)
    done; } > $W/locus_methods.tsv
  touch $W/SHRINK.DONE
fi
branch ml $W/ML.DONE
# ---- stage 3: wASTRAL species tree, quartet details, concatenation
if [ ! -e $W/WASTRAL.DONE ]; then
  log "wASTRAL (-B aBayes, hybrid weighting, -R) on $(wc -l < $GTS) gene trees"
  setsid wastral -B --mode 1 -R -u 1 -t 64 --seed 20261019 -i $GTS -o $W/wastral.tre > $W/log/wastral.log 2>&1 \
    || { log "wastral failed"; tail -n 5 $W/log/wastral.log; exit 1; }
  log "wASTRAL score $(qscore $W/log/wastral.log)"
  $TOOLS topo $W/wastral.tre $W/wastral_topo.tre
  mkdir -p $W/wastral_score
  for u in 2 3; do
    ( cd $W/wastral_score && setsid wastral -C -B --mode 1 -u $u -t 32 --seed 20261019 -c $W/wastral_topo.tre -i $GTS \
        -o $W/wastral_score/u$u.tre > $W/log/wastral_u$u.log 2>&1 ) || log "wastral -u $u failed (non-fatal)"
  done
  $FT concat $W > $W/log/concat.txt 2>&1 || { log "concat failed"; cat $W/log/concat.txt; exit 1; }
  log "$(tail -n 1 $W/log/concat.txt)"
  touch $W/WASTRAL.DONE
fi
# ---- stage 4: fixed topology (SPEC ops on the unconstrained tree) + taxonomic QC gate
if [ ! -e $F/TOPO.DONE ]; then
  $PY scripts/19g_topo_qc.py $W/wastral_topo.tre $W/qc_unconstrained $OG > $W/log/qc_unconstrained.txt 2>&1
  log "QC of the unconstrained wASTRAL tree: $(tail -n 1 $W/log/qc_unconstrained.txt)"
  while :; do
    m0=$(specmd5); ok=0
    if $PY scripts/19g_topo_fix.py $W/wastral_topo.tre $F/fix $OG $SPEC > $F/log/fix.txt 2>&1; then
      mv $F/fix.tre $F/topo.tre; mv $F/fix_unrooted.tre $F/topo_unrooted.tre; mv $F/fix.log.tsv $F/topo_fix.log.tsv
      log "fixed topology: $(tail -n 1 $F/log/fix.txt)"
      $PY scripts/19g_topo_qc.py $F/topo.tre $F/qc $OG > $F/log/qc.txt 2>&1
      # non-monophyletic groups that are not in ACC; those above order/family level stop the pipeline
      awk -F'\t' 'NR == FNR {if (FNR > 1) a[$1 "\t" $2]; next} FNR == 1 || ($4 == 0 && !(($1 "\t" $2) in a))' $ACC $F/qc.monophyly.tsv > $F/qc_new.tsv
      nh=$(awk -F'\t' 'NR > 1 && $1 != "order" && $1 != "family"' $F/qc_new.tsv | wc -l)
      log "QC: $(tail -n 1 $F/log/qc.txt); not accepted: $(($(wc -l < $F/qc_new.tsv) - 1)) ($nh above order/family)"
      if [ $nh = 0 ] || [ -e $F/TOPO.OK ]; then rm -f $F/REVIEW_NEEDED; break; fi
      { echo "non-monophyletic groups above order/family level that are not in $ACC:"
        awk -F'\t' 'NR > 1 && $1 != "order" && $1 != "family"' $F/qc_new.tsv; } > $F/REVIEW_NEEDED
      ok=1
    else
      { echo "topology fix failed ($SPEC):"; tail -n 5 $F/log/fix.txt; } > $F/REVIEW_NEEDED
    fi
    log "REVIEW_NEEDED ($F/REVIEW_NEEDED): edit $SPEC$([ $ok = 1 ] && echo " or touch $F/TOPO.OK to accept the QC")"
    while [ "$(specmd5)" = "$m0" ] && { [ $ok = 0 ] || [ ! -e $F/TOPO.OK ]; }; do sleep 600; done
  done
  printf "spec\t%s\t%s\naccepted\t%s\t%s\nwastral_topo\t%s\t%s\n" $(readlink -f $SPEC) $(specmd5) $(readlink -f $ACC) \
    $(md5sum < $ACC | cut -c1-32) $W/wastral_topo.tre $(md5sum < $W/wastral_topo.tre | cut -c1-32) > $F/topo_inputs.tsv
  touch $F/TOPO.DONE
fi
branch rer36 $RV/TREES.DONE; branch unc $U/CHAIN.DONE $U/SKIPPED
# ---- stage 6: quartet scores (background), branch lengths on the fixed topology, rooting, gCF
scores() { local k c
  for k in fixed unconstrained; do
    c=$F/topo.tre; [ $k = unconstrained ] && c=$W/wastral_topo.tre
    wastral -C -B --mode 1 -u 2 -t 8 --seed 20261019 -c $c -i $GTS -o $F/score_$k.tre > $F/log/score_$k.log 2>&1
  done
  [ -n "$(qscore $F/log/score_fixed.log)" ] && [ -n "$(qscore $F/log/score_unconstrained.log)" ] \
    && printf "topology\tscore\nfixed\t%s\nunconstrained\t%s\n" "$(qscore $F/log/score_fixed.log)" \
         "$(qscore $F/log/score_unconstrained.log)" > $F/score.tsv; }
SCPID=""; [ -s $F/score.tsv ] || { scores & SCPID=$!; }
if [ ! -e $F/BL.DONE ]; then
  log "branch lengths on the fixed topology (IQ-TREE -p partition.nex -te)"
  iqrun $F/bl treefile -p $W/partition.nex -te $F/topo_unrooted.tre -T $(thr 48) -seed 20261019 || exit 1
  $FT samesplits $F/bl.treefile $F/topo.tre > $F/log/bl_check.txt 2>&1 || { log "IQ-TREE changed the fixed topology"; exit 1; }
  echo "IQ-TREE $(iqver) -p partition.nex (per-locus ModelFinder models, edge-proportional) -te fixed topology ($SPEC md5 $(specmd5)), seed 20261019" > $F/bl_method.txt
  $TOOLS root $F/bl.treefile $F/species_tree_og.nwk $F/species_tree.nwk > $F/log/root.txt 2>&1 || { log "rooting failed: $(tail -n 1 $F/log/root.txt)"; exit 1; }
  gcf $F/bl.treefile gcf || log "gCF failed (non-fatal; branch cf retries)"
  [ -n "$SCPID" ] && { wait $SCPID; SCPID=""; }
  [ -s $F/score.tsv ] || log "quartet scores failed (non-fatal; $F/log/score_*.log)"
  touch $F/BL.DONE
  log "BL.DONE: $(grep "Log-likelihood of the tree" $F/bl.iqtree); quartet scores: $(tail -n +2 $F/score.tsv 2>/dev/null | tr '\n\t' '; ')"
fi
[ -n "$SCPID" ] && wait $SCPID
branch cf $F/CF.DONE; branch mlsens $M/CHAIN.DONE $M/SKIPPED
# ---- stage 7: downstream chain on the fixed tree (dating, ASR, PGLS, RERconverge; all label sets)
if [ ! -e $F/CHAIN.DONE ]; then
  while running rer36; do sleep 300; done
  bash scripts/36_rer_trees.sh $W $F/topo.tre $RV 96 >> $F/log/36.log 2>&1 || { log "RER gene trees failed ($F/log/36.log)"; exit 1; }
  log "main chain on the fixed tree (all label sets, NPERM $NPERM_MAIN)"
  NPERM=$NPERM_MAIN NCORES=48 bash $CHAIN $F $F/species_tree_og.nwk $W/concat.faa $RV $RV/rer all >> $F/log/38c.log 2>&1 \
    || { log "main chain failed ($F/log/38c.log)"; exit 1; }
  log "CHAIN.DONE"
fi
[ -s $F/score.tsv ] || scores
# ---- wait for the branches; restart a failed one once
sleep 30
for round in 1 2; do
  while running ml || running rer36 || running unc || running cf || running mlsens; do sleep 600; done
  miss=""
  [ -e $W/ML.DONE ] || miss="$miss ml"
  [ -e $RV/TREES.DONE ] || miss="$miss rer36"
  [ -e $U/CHAIN.DONE ] || [ -e $U/SKIPPED ] || miss="$miss unc"
  [ -e $F/CF.DONE ] || miss="$miss cf"
  [ -e $M/CHAIN.DONE ] || [ -e $M/SKIPPED ] || miss="$miss mlsens"
  [ -z "$miss" ] && break
  [ $round = 2 ] && { log "not finished:$miss (see $W/log/post_*.log); fix, then rerun this script"; exit 1; }
  log "restarting once:$miss"
  for m in $miss; do branch $m; done
done
[ -s $F/score.tsv ] || log "WARNING: no quartet scores ($F/log/score_*.log)"
touch $W/TREE.DONE
log "TREE.DONE"
