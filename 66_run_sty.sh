#!/bin/bash
# R1 depth-threshold, R2 GO spike-in, R3 fish RER GO + RTF1 aa scan (PREREG_nature items not yet run).
cd ~/deepsea_conv || exit 0
FUNC=work/func_v2; [ -d $FUNC/Mollusca ] || FUNC=work/func
echo "FUNC=$FUNC"
L=work/addon/logs/66_run.log
{
date; echo "== 66a"; bash scripts/memguard.sh 48 python scripts/66a_fgdepth_sty.py
echo "== 66b"; bash scripts/memguard.sh 48 Rscript scripts/66b_depth_thresh_sty.R
echo "== 66c"; bash scripts/memguard.sh 48 python scripts/66c_fishrer_go_sty.py $FUNC
echo "== 66e"; bash scripts/memguard.sh 48 python scripts/66e_rtf1_aascan_sty.py
echo "== 66d"; bash scripts/memguard.sh 48 python scripts/66d_go_spikein_sty.py $FUNC
date; echo END
} > $L 2>&1
touch work/addon/R66.DONE
exit 0
