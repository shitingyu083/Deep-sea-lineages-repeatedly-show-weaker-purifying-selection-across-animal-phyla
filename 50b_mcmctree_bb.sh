#!/bin/bash
# MCMCTree on a calibration-preserving backbone (2026-10-25; PREREG_addon A6; DEVIATIONS_fig4 2026-10-25b). Copy of
# 50_mcmctree_run.sh (unchanged ctl except sampfreq 50 and the explicit 'usedata = 2 in.BV'); the full-tree run needed
# 2.1 s per iteration (~49 d per chain) and was stopped 2026-10-25 13:10.
# 1 per N (NSERIES, default 120 100 80) in OUT/n<N>: 50b_backbone_select.py -> bb.timetree.nwk; 50_mcmctree_prep.py on
#   it; the applied calibration set must equal the full tree's (work/mcmctree/calibrations.tsv; else CAL_MISMATCH)
#   -> PREP.DONE; Hessian as 50_mcmctree_run.sh (usedata 3, codeml LG+G4 per partition) -> HESS.DONE
# 2 speed probe: 200 iterations of the c2_r1 ctl (wall time incl. start-up = conservative); projected chain time =
#   s/iter x (burnin + sampfreq x nsample) x 1.5 (five chains on a shared server). The first N with a projection
#   <= 7 d is used; if none, N = 80 with sampfreq 25 -> OUT/N_USED ("N sampfreq"), OUT/use -> n<N>
# 3 usedata 2: clock 2 chains r1 r2 (seeds 20261025, 20261026), clock 3 chains r1 r2 (20261027, 20261029), prior
#   (usedata 0, clock 2, 20261028); burnin 20000, nsample 20000 -> <run>/RUN.DONE
# 4 50b_mcmctree_collect.py; if the clock-2 chains fail PREREG_fig4 C, one extension: the four posterior chains again
#   with doubled sampfreq and new seeds (<run>x), collected with --suffix x (reported either way) -> OUT/BB.DONE
# Usage: bash 50b_mcmctree_bb.sh W [OUT]
cd ~/deepsea_conv || exit 1
echo 10 > /proc/self/autogroup 2>/dev/null
W=$(readlink -f $1); OUT=$(readlink -f ${2:-work/mcmctree_bb}); mkdir -p $OUT/logs
PY=envs/dc/bin/python; BIN=$(readlink -f envs/dc/bin); FULLCAL=$(readlink -f work/mcmctree/calibrations.tsv)
log() { echo "$(date +%F_%T) $*" >> $OUT/logs/run.log; }
LG=$(find envs/dc -name lg.dat 2>/dev/null | head -1); LG=${LG:+$(readlink -f $LG)}
[ -n "$LG" ] || { log "lg.dat not found"; exit 1; }
[ -s "$FULLCAL" ] || { log "full-tree calibrations.tsv missing"; exit 1; }
NSERIES=${NSERIES:-"120 100 80"}; LIMIT=${LIMIT_S:-604800}; LOADF=1.5
BURN=20000; SF=50; NSAMP=20000
ctl() {  # usedata clock seed
  cat <<EOF
seed = $3
seqfile = ../aln.phy
treefile = ../tree.tre
mcmcfile = mcmc.txt
outfile = out.txt
ndata = $NDATA
seqtype = 2
usedata = $1$( [ "$1" = 2 ] && echo "$BVSFX" )
clock = $2
$( [ -n "$ROOTAGE" ] && echo "RootAge = $ROOTAGE" )
model = 0
alpha = 0
ncatG = 5
cleandata = 0
BDparas = 1 1 0.1 m
kappa_gamma = 6 2
alpha_gamma = 1 1
rgene_gamma = $RGENE
sigma2_gamma = 1 10 1
finetune = 1: .1 .1 .1 .1 .1 .1
print = 1
burnin = $BURN
sampfreq = $SF
nsample = $NSAMP
EOF
}
# the usedata 2 line exactly as in the full-tree chains (PAML reads in.BV from the run directory by default)
REFCTL=$(ls work/mcmctree/stopped_fulltree_*/c2_r1/mcmctree.ctl work/mcmctree/c2_r1/mcmctree.ctl 2>/dev/null | head -1)
BVSFX=$( [ -n "$REFCTL" ] && grep -E '^ *usedata *= *2' $REFCTL | head -1 | sed -E 's/^ *usedata *= *2//; s/ *$//' )
BVSFX=${BVSFX:-" in.BV"}   # without the file name mcmctree stops with "file name empty" (ck96)
log "start N series '$NSERIES'; usedata 2 line from ${REFCTL:-default}: 'usedata = 2$BVSFX'"
prep_n() {  # N -> OUT/n<N>/{PREP.DONE, HESS.DONE}; 1 on failure
  local N=$1 D=$OUT/n$1 H=$OUT/n$1/hess n c b
  if [ ! -e $D/PREP.DONE ]; then
    rm -rf $D; mkdir -p $D
    $PY scripts/50b_backbone_select.py $W --out $D --n $N > $D/select.log 2>&1 || { log "n$N select FAIL: $(tail -2 $D/select.log | tr '\n' ' ')"; return 1; }
    $PY scripts/50_mcmctree_prep.py $W --out $D --timetree $D/bb.timetree.nwk > $D/prep.log 2>&1 || { log "n$N prep FAIL"; return 1; }
    $PY scripts/50b_backbone_select.py --check $FULLCAL $D/calibrations.tsv > $D/calcheck.txt 2>&1 || {
      touch $D/CAL_MISMATCH; log "n$N CAL_MISMATCH: $(tail -1 $D/calcheck.txt)"; return 1; }
    touch $D/PREP.DONE; log "n$N prep: $(head -1 $D/prep_report.txt); $(tail -1 $D/calcheck.txt)"
  fi
  . $D/params.sh
  [ -e $D/HESS.DONE ] && return 0
  rm -rf $H; mkdir -p $H/fake; cd $H
  printf '#!/bin/sh\necho 0 > rst2\nexit 0\n' > fake/codeml; chmod +x fake/codeml
  ctl 3 2 20261025 > mcmctree.ctl
  PATH=$H/fake:$PATH $BIN/mcmctree mcmctree.ctl > mcmctree_hess.log 2>&1
  n=$(ls tmp0*.ctl 2>/dev/null | wc -l)
  [ "$n" = "$NDATA" ] || { log "n$N usedata=3 wrote $n tmp ctl (want $NDATA)"; cd ~/deepsea_conv; return 1; }
  for c in tmp0*.ctl; do
    b=${c%.ctl}; mkdir -p $b; mv $b.* $b/ 2>/dev/null; cd $b
    sed -i -E -e 's/^( *model *=).*/\1 2/' -e 's/^( *fix_alpha *=).*/\1 0/' -e 's/^( *alpha *=).*/\1 0.5/' \
      -e 's/^( *ncatG *=).*/\1 4/' -e 's/^( *method *=).*/\1 1/' $b.ctl
    grep -q aaRatefile $b.ctl && sed -i -E "s#^( *aaRatefile *=).*#\1 $LG#" $b.ctl || echo "aaRatefile = $LG" >> $b.ctl
    grep -q fix_alpha $b.ctl || echo "fix_alpha = 0" >> $b.ctl
    grep -q ncatG $b.ctl || echo "ncatG = 4" >> $b.ctl
    grep -qE '^ *alpha *=' $b.ctl || echo "alpha = 0.5" >> $b.ctl
    (nice -n 10 $BIN/codeml $b.ctl > codeml.log 2>&1; echo "$(date +%F_%T) n$N $b codeml rc=$?" >> $OUT/logs/run.log) &
    cd $H
  done
  wait
  rm -f in.BV
  for c in tmp0*/; do b=${c%/}; [ -s $b/rst2 ] || { log "n$N $b no rst2"; cd ~/deepsea_conv; return 1; }; cat $b/rst2 >> in.BV; done
  cd ~/deepsea_conv; touch $D/HESS.DONE; log "n$N in.BV from $NDATA partitions"
}
probe() {  # D -> wall seconds per iteration over 200 iterations incl. start-up (99999 on failure)
  local P=$1/probe t0 t1 rc; rm -rf $P; mkdir -p $P; cd $P; cp $1/hess/in.BV .
  ( BURN=10; SF=1; NSAMP=190; ctl 2 2 20261025 ) > mcmctree.ctl
  t0=$(date +%s.%N); timeout 7200 nice -n 10 $BIN/mcmctree mcmctree.ctl > mcmctree.log 2>&1; rc=$?; t1=$(date +%s.%N)
  cd ~/deepsea_conv
  if [ $rc = 0 ] && [ -s $P/mcmc.txt ]; then awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.4f\n", (b - a) / 200}'; else echo 99999; fi
}
if [ ! -e $OUT/N_USED ]; then
  [ -s $OUT/speed.tsv ] || printf 'N\ts_per_iter\tprojected_h\tlimit_h\n' > $OUT/speed.tsv
  for N in $NSERIES; do
    prep_n $N || continue
    s=$(probe $OUT/n$N)
    p=$(awk -v s=$s -v b=$BURN -v f=$SF -v n=$NSAMP -v l=$LOADF 'BEGIN{printf "%.0f", s * (b + f * n) * l}')
    printf '%s\t%s\t%.1f\t%.1f\n' $N $s $(awk -v p=$p 'BEGIN{print p / 3600}') $(awk -v p=$LIMIT 'BEGIN{print p / 3600}') >> $OUT/speed.tsv
    log "n$N probe $s s/iter -> projected ${p} s per chain (x$LOADF; limit $LIMIT s)"
    if [ "$p" -le "$LIMIT" ]; then echo "$N $SF" > $OUT/N_USED; break; fi
  done
  if [ ! -e $OUT/N_USED ]; then
    last=${NSERIES##* }
    [ -e $OUT/n$last/HESS.DONE ] || { log "no backbone prepared (see logs); stop"; exit 1; }
    echo "$last 25" > $OUT/N_USED; log "no N within the limit: N=$last with sampfreq 25"
  fi
fi
read N SF < $OUT/N_USED
D=$OUT/n$N; ln -sfn n$N $OUT/use; . $D/params.sh
log "using n$N, sampfreq $SF (burnin $BURN, nsample $NSAMP)"
run() {  # name usedata clock seed
  local X=$D/$1 rc; [ -e $X/RUN.DONE ] && return 0
  rm -rf $X; mkdir -p $X; cd $X; cp $D/hess/in.BV .
  ctl $2 $3 $4 > mcmctree.ctl
  nice -n 10 $BIN/mcmctree mcmctree.ctl > mcmctree.log 2>&1; rc=$?
  cd ~/deepsea_conv; log "n$N $1 rc=$rc"; [ $rc = 0 ] && [ -s $X/mcmc.txt ] && touch $X/RUN.DONE
}
run c2_r1 2 2 20261025 & run c2_r2 2 2 20261026 & run c3_r1 2 3 20261027 & run c3_r2 2 3 20261029 & run prior 0 2 20261028 &
wait
$PY scripts/50b_mcmctree_collect.py $D > $D/collect.log 2>&1; log "n$N collect rc=$?: $(head -1 $D/convergence.txt 2>/dev/null)"
USE=""
if ! head -1 $D/convergence.txt 2>/dev/null | grep -q '^PASS'; then
  # one extension (PREREG_fig4 C 'else extend'): the four posterior chains with doubled sampfreq and new seeds
  SF=$((SF * 2)); log "n$N clock-2 chains fail PREREG_fig4 C -> extension, sampfreq $SF"
  run c2_r1x 2 2 20261125 & run c2_r2x 2 2 20261126 & run c3_r1x 2 3 20261127 & run c3_r2x 2 3 20261129 &
  wait; touch $D/EXT.DONE
  $PY scripts/50b_mcmctree_collect.py $D --suffix x > $D/collect_x.log 2>&1
  log "n$N collect x rc=$?: $(head -1 $D/convergencex.txt 2>/dev/null)"; USE=x
fi
# 50c reads draws<USE>.tsv; the verdict is reported either way
echo "${USE:-none} $(head -1 $D/convergence$USE.txt 2>/dev/null | cut -d' ' -f1)" > $D/USE
touch $OUT/BB.DONE; log "BB.DONE n$N suffix '${USE}' $(cat $D/USE)"
exit 0
