#!/bin/bash
# MCMCTree dating, approximate likelihood (2026-10-25).
# 1 prep (50_mcmctree_prep.py W)                                         -> OUT/PREP.DONE
# 2 gradient + Hessian: mcmctree usedata=3 writes tmp000k.ctl (codeml is replaced by a stub that writes an empty rst2 while it does), each
#   tmp000k.ctl is set to LG+G4 (model=2 lg.dat, alpha 0.5 estimated, ncatG 4), codeml runs per partition in
#   parallel, rst2 files are concatenated to in.BV                     -> OUT/HESS.DONE
# 3 usedata=2 runs: clock 2 (independent rates) chains r1 r2 (seeds 20261025, 20261026), clock 3 (autocorrelated)
#   chain r1, and the joint prior (usedata 0, clock 2); burnin 20000, sampfreq 100, nsample 20000 -> OUT/<run>/RUN.DONE
# 4 50_mcmctree_collect.py                                              -> OUT/MCMC.DONE
# Usage: bash 50_mcmctree_run.sh W [OUT]
cd ~/deepsea_conv
W=$(readlink -f $1); OUT=$(readlink -f ${2:-work/mcmctree}); mkdir -p $OUT/logs
PY=envs/dc/bin/python; BIN=$(readlink -f envs/dc/bin)
log() { echo "$(date +%F_%T) $*" >> $OUT/logs/run.log; }
LG=$(find envs/dc -name lg.dat 2>/dev/null | head -1); LG=${LG:+$(readlink -f $LG)}
[ -n "$LG" ] || { log "lg.dat not found"; exit 1; }
if [ ! -e $OUT/PREP.DONE ]; then
  $PY scripts/50_mcmctree_prep.py $W --out $OUT > $OUT/logs/prep.log 2>&1 || { log "prep FAIL"; exit 1; }
  touch $OUT/PREP.DONE; log "prep: $(head -1 $OUT/prep_report.txt)"
fi
. $OUT/params.sh
ctl() {  # usedata clock seed
  cat <<EOF
seed = $3
seqfile = ../aln.phy
treefile = ../tree.tre
mcmcfile = mcmc.txt
outfile = out.txt
ndata = $NDATA
seqtype = 2
usedata = $1$( [ "$1" = 2 ] && echo " in.BV" )
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
burnin = 20000
sampfreq = 100
nsample = 20000
EOF
}
H=$OUT/hess
if [ ! -e $OUT/HESS.DONE ]; then
  rm -rf $H; mkdir -p $H/fake; cd $H
  printf '#!/bin/sh\necho 0 > rst2\nexit 0\n' > fake/codeml; chmod +x fake/codeml
  ctl 3 2 20261025 > mcmctree.ctl
  PATH=$H/fake:$PATH $BIN/mcmctree mcmctree.ctl > mcmctree_hess.log 2>&1
  n=$(ls tmp0*.ctl 2>/dev/null | wc -l)
  [ "$n" = "$NDATA" ] || { log "usedata=3 wrote $n tmp ctl (want $NDATA), see $H/mcmctree_hess.log"; exit 1; }
  for c in tmp0*.ctl; do
    b=${c%.ctl}; mkdir -p $b; mv $b.* $b/ 2>/dev/null; cd $b
    sed -i -E -e 's/^( *model *=).*/\1 2/' -e 's/^( *fix_alpha *=).*/\1 0/' -e 's/^( *alpha *=).*/\1 0.5/' \
      -e 's/^( *ncatG *=).*/\1 4/' -e 's/^( *method *=).*/\1 1/' $b.ctl
    grep -q aaRatefile $b.ctl && sed -i -E "s#^( *aaRatefile *=).*#\1 $LG#" $b.ctl || echo "aaRatefile = $LG" >> $b.ctl
    grep -q fix_alpha $b.ctl || echo "fix_alpha = 0" >> $b.ctl
    grep -q ncatG $b.ctl || echo "ncatG = 4" >> $b.ctl
    grep -qE '^ *alpha *=' $b.ctl || echo "alpha = 0.5" >> $b.ctl
    (nice -n 10 $BIN/codeml $b.ctl > codeml.log 2>&1; echo "$(date +%F_%T) $b codeml rc=$?" >> $OUT/logs/run.log) &
    cd $H
  done
  wait
  cd $H; rm -f in.BV
  for c in tmp0*/; do b=${c%/}; [ -s $b/rst2 ] || { log "$b no rst2"; exit 1; }; cat $b/rst2 >> in.BV; done
  cd ~/deepsea_conv; touch $OUT/HESS.DONE; log "in.BV from $NDATA partitions"
fi
run() {  # name usedata clock seed
  local D=$OUT/$1; [ -e $D/RUN.DONE ] && return 0
  rm -rf $D; mkdir -p $D; cd $D; cp $H/in.BV .
  ctl $2 $3 $4 > mcmctree.ctl
  nice -n 10 $BIN/mcmctree mcmctree.ctl > mcmctree.log 2>&1; local rc=$?
  cd ~/deepsea_conv; log "$1 rc=$rc"; [ $rc = 0 ] && [ -s $D/mcmc.txt ] && touch $D/RUN.DONE
}
run c2_r1 2 2 20261025 & run c2_r2 2 2 20261026 & run c3_r1 2 3 20261027 & run prior 0 2 20261028 &
wait
$PY scripts/50_mcmctree_collect.py $OUT > $OUT/logs/collect.log 2>&1
touch $OUT/MCMC.DONE; log "MCMC.DONE"
exit 0
