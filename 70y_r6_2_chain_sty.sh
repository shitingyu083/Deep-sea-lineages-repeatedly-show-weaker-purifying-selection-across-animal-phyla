#!/bin/bash
# R6-2 chain (ADDENDUM_R6_2.md md5 0e9d4fe6): wait R6_1.DONE -> 70j prep -> 70l fits (memguard 48, NJ 64)
# -> 70j collect -> 70k stats -> R6_2.DONE. Log R6/R6_2/chain.log.
R=~/deepsea_conv; O=$R/work/addon/R6; G=$O/R6_2; cd $R
mkdir -p $G; L=$G/chain.log
echo "$(date +%F_%T) chain start pid $$" >> $L
while [ ! -f $O/R6_1.DONE ]; do sleep 300; done
echo "$(date +%F_%T) R6_1.DONE seen" >> $L
md5sum $O/ADDENDUM_R6_2.md >> $L
$R/envs/dc/bin/python scripts/70j_r6_group_sty.py prep >> $L 2>&1 || { echo "prep FAIL" >> $L; exit 1; }
echo "$(date +%F_%T) prep ok jobs $(wc -l < $G/jobs.txt)" >> $L
bash scripts/memguard.sh 48 bash scripts/70l_r6_group_relax_sty.sh 64 >> $L 2>&1
echo "$(date +%F_%T) fits end rc $?" >> $L
$R/envs/dc/bin/python scripts/70j_r6_group_sty.py collect >> $L 2>&1 || { echo "collect FAIL" >> $L; exit 1; }
$R/envs/rphy/bin/Rscript scripts/70k_r6_2_stats_sty.R $O >> $L 2>&1 || { echo "stats FAIL" >> $L; exit 1; }
touch $G/R6_2.DONE; echo "$(date +%F_%T) R6_2.DONE" >> $L
exit 0
