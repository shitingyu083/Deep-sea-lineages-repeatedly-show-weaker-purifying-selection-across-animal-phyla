#!/bin/bash
# Keep the panel-v2 daemons alive until annotation of every species_codes row is done.
#  - 07_qc_loop_v2 exits once DOWNLOAD*.DONE exist and nothing is queued, but the v2
#    downloader (05d ... acc_list_v2.tsv) is still adding genomes -> restart it.
#  - 14_protqc / 17_emapper stop only after work/ANNOT.DONE -> restart if they died early.
#  - 11_annot_loop rereads species_codes.tsv -> restart only if it died without ANNOT.DONE.
# Usage: setsid nohup bash scripts/27_watchdog_v2.sh > logs/watchdog_v2.log 2>&1 < /dev/null &
R=~/deepsea_conv; cd $R
up() { pgrep -f "$1" > /dev/null; }
while true; do
  dl=0; up "05d_download_staged.sh data/ncbi/acc_list_v" && dl=1
  # downloads still landing, or finished downloads without a compleasm summary
  need_qc=$dl
  [ $need_qc = 0 ] && for a in $(cut -f1 data/ncbi/acc_list_v2.tsv data/ncbi/acc_list_v3.tsv); do
    [ -f data/ncbi/$a/.done ] && [ ! -s work/qc/$a/summary.txt ] && [ ! -f work/qc/$a/.fail ] && { need_qc=1; break; }
  done
  if [ $need_qc = 1 ] && ! up "07_qc_loop_v2.sh"; then
    rm -f work/qc/QC.DONE
    echo "$(date +%F_%T) restart 07_qc_loop_v2"
    setsid nohup bash scripts/07_qc_loop_v2.sh 6 16 >> logs/qc.log 2>&1 < /dev/null &
  fi
  if [ ! -f work/ANNOT.DONE ]; then
    up "11_annot_loop.sh" || { echo "$(date +%F_%T) restart 11_annot_loop";
      setsid nohup bash scripts/11_annot_loop.sh 6 24 >> logs/annot_loop.log 2>&1 < /dev/null & }
    up "14_protqc_loop.sh" || { echo "$(date +%F_%T) restart 14_protqc_loop";
      setsid nohup bash scripts/14_protqc_loop.sh 4 8 >> logs/protqc.log 2>&1 < /dev/null & }
    up "17_emapper_loop.sh" || { echo "$(date +%F_%T) restart 17_emapper_loop";
      setsid nohup bash scripts/17_emapper_loop.sh 4 32 >> logs/emapper.log 2>&1 < /dev/null & }
  elif [ $dl = 0 ] && [ $need_qc = 0 ]; then
    echo "$(date +%F_%T) all v2 download/QC/annotation finished"; break
  fi
  sleep 600
done
