#!/bin/bash
# R2 pilot gate rerun with v2 (per-HOG tree-length calibration; DEVIATIONS_addon.md D-R2-1).
# Same 3 Echinodermata HOGs x {H0r,K50}; gate: K-hat(K50) in [0.3,0.75], refit H0r tree length within +-20% of fitL.
B=$HOME/deepsea_conv; PY=$B/envs/dc/bin/python; IQ=$B/envs/dc/bin/iqtree
R=$B/work/addon/R2/pilot2; mkdir -p $R; cd $R || exit 0
cp $B/work/addon/R2/pilot/hogs.txt hogs.txt
HOGS=$(cat hogs.txt)
for h in $HOGS; do echo "$PY $B/scripts/63a2_powersim_sty.py calib Echinodermata $h $R/$h/_fit"; done | xargs -P 3 -I{} bash -c "{} < /dev/null" > calib.log 2>&1
for h in $HOGS; do for s in H0r K50; do
  echo "$PY $B/scripts/63a2_powersim_sty.py Echinodermata $h $s 0 $R/$h/$s"
done; done > jobs.txt
cat jobs.txt | xargs -P 6 -I{} bash -c "{} < /dev/null" > run.log 2>&1
for h in $HOGS; do
  gzip -dc $R/$h/H0r/sim.fas.gz > $R/$h/H0r/sim_refit.fas
  $IQ -s $R/$h/H0r/sim_refit.fas -te $R/$h/_fit/topo.nwk -m GY+F3X4 -st CODON -T 1 --prefix $R/$h/H0r/refit -redo -quiet > /dev/null 2>&1
  rm -f $R/$h/H0r/sim_refit.fas
  a=$(grep "Total tree length" $R/$h/_fit/fit.iqtree | awk '{print $NF}'); b=$(grep "Total tree length" $R/$h/H0r/refit.iqtree | awk '{print $NF}')
  sc=$(cut -f1 $R/$h/_fit/scale.txt)
  echo -e "$h\tfitL\t$a\tsimL\t$b\tscale\t$sc"
done > treelen.tsv
cat */*/res.tsv | awk 'NR==1||$1!="clade"' > pilot_res.tsv
touch PILOT2.DONE
exit 0
