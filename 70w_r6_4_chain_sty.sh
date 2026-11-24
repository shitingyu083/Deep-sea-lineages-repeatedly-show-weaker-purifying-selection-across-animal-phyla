#!/bin/bash
# R6-2b + R6-4 chain (ADDENDUM_R6_4.md, md5 in work/addon/FROZEN_r6_4.txt).
# 70p covariates -> 70m prep -> 70n fits (memguard 48, NJ 64) -> 70m collect -> 70o stats B / AB / A -> 70q R6-4.
# Log: work/addon/R6/R6_4_chain.log; ends with R6_4.DONE (or R6_4.FAIL).
R=~/deepsea_conv; cd $R; O=$R/work/addon/R6; L=$O/R6_4_chain.log
PY=$R/envs/dc/bin/python; RS=$R/envs/rphy/bin/Rscript
say() { echo "$(date +%F_%T) $*" >> $L; }
die() { say "FAIL $*"; touch $O/R6_4.FAIL; exit 1; }
say "start"
$PY scripts/70p_r6_lineage_cov_sty.py >> $L 2>&1 || die "70p"
$PY scripts/70m_r6_group_b_sty.py prep >> $L 2>&1 || die "70m prep"
say "prep ok jobs $(wc -l < $O/R6_2b/jobs.txt)"
bash scripts/memguard.sh 48 bash scripts/70n_r6_group_relax_b_sty.sh 64 >> $L 2>&1
[ -f $O/R6_2b/RELAX.END ] || die "70n (no RELAX.END; memguard?)"
$PY scripts/70m_r6_group_b_sty.py collect >> $L 2>&1 || die "70m collect"
$RS scripts/70o_r6_2_stats_b_sty.R $O/R6_2b/stats B $O/R6_2b >> $L 2>&1 || die "70o B"
$RS scripts/70o_r6_2_stats_b_sty.R $O/R6_2ab AB $O/R6_2 $O/R6_2b >> $L 2>&1 || die "70o AB"
$RS scripts/70o_r6_2_stats_b_sty.R $O/R6_4/setA A $O/R6_2 >> $L 2>&1 || die "70o A"
# set A via 70o must reproduce the frozen R6-2 summary (70k)
cmp -s $O/R6_4/setA/r6_2_summary.tsv $O/R6_2/r6_2_summary.tsv && say "setA reproduces R6_2 summary: IDENTICAL" \
  || say "setA vs R6_2 summary: DIFFERENT (see diff below)"
diff $O/R6_4/setA/r6_2_summary.tsv $O/R6_2/r6_2_summary.tsv >> $L 2>&1
$RS scripts/70q_r6_4_stats_sty.R $O/R6_4/lineage_covariates.tsv $O/R6_4 \
  A=$O/R6_4/setA B=$O/R6_2b/stats AB=$O/R6_2ab >> $L 2>&1 || die "70q"
touch $O/R6_4.DONE; say "R6_4.DONE"
exit 0
