#!/bin/bash
# R6-1 chain: wait for REANNOT.END -> 70c build -> 70d RELAX (4 fits/HOG) -> 70e collect -> 70h stats.
# Stops (no stats) if a step fails or memguard kills it. Log: work/addon/R6/chain.log
cd ~/deepsea_conv; O=work/addon/R6; L=$O/chain.log
say() { echo "$(date +%F_%T) $*" >> $L; }
say "chain start pid $$"
while [ ! -f $O/REANNOT.END ]; do sleep 120; done
say "REANNOT.END seen; done $(ls $O/reannot/*/REANNOT.DONE 2>/dev/null | wc -l) fail $(ls $O/reannot/*/REANNOT.FAIL 2>/dev/null | wc -l)"
if [ ! -f $O/BUILD.DONE ]; then
  bash scripts/memguard.sh 48 envs/dc/bin/python scripts/70c_r6_build_sty.py 32 > $O/build_nohup.txt 2>&1
  rc=$?; say "70c rc $rc"; [ $rc -eq 0 ] && [ -f $O/BUILD.DONE ] || { say "STOP build"; exit 0; }
fi
say "built ok $(awk -F'\t' 'NR>1 && $7=="ok"' $O/build_report.tsv | wc -l)"
bash scripts/memguard.sh 48 bash scripts/70d_r6_relax_sty.sh 64 > $O/relax_nohup.txt 2>&1
rc=$?; say "70d rc $rc"; [ $rc -eq 0 ] || { say "STOP relax"; exit 0; }
envs/dc/bin/python scripts/70e_r6_collect_sty.py > $O/collect_nohup.txt 2>&1; say "70e rc $?"
envs/rphy/bin/Rscript scripts/70h_r6_route_stats_sty.R $O > $O/stats_nohup.txt 2>&1; say "70h rc $?"
touch $O/R6_1.DONE; say "chain end"
exit 0
