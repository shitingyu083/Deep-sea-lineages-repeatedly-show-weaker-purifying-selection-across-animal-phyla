#!/bin/bash
# R2 pilot gate (PREREG_nature.md): 3 Echinodermata HOGs x {H0r,K50}; check K-hat and tree length.
B=$HOME/deepsea_conv; PY=$B/envs/dc/bin/python; IQ=$B/envs/dc/bin/iqtree
R=$B/work/addon/R2/pilot; mkdir -p $R; cd $R || exit 0
W=$B/work/codon/Echinodermata
HOGS=$(ls $W/sel/*.relax.json | grep -v null | xargs -n1 basename | sed 's/.relax.json//' | head -3)
echo "$HOGS" > hogs.txt
for h in $HOGS; do for s in H0r K50; do
  echo "$PY $B/scripts/63a_powersim_sty.py Echinodermata $h $s 0 $R/$h/$s"
done; done > jobs.txt
cat jobs.txt | xargs -P 6 -I{} bash -c "{} < /dev/null" > run.log 2>&1
for h in $HOGS; do
  $IQ -s $R/$h/H0r/sim.fas -te $R/$h/_fit/topo.nwk -m GY+F3X4 -st CODON -T 1 --prefix $R/$h/H0r/refit -redo -quiet > /dev/null 2>&1
  a=$(grep "Total tree length" $R/$h/_fit/fit.iqtree | awk '{print $NF}'); b=$(grep "Total tree length" $R/$h/H0r/refit.iqtree | awk '{print $NF}')
  echo -e "$h\tfitL\t$a\tsimL\t$b"
done > treelen.tsv
cat */*/res.tsv | awk 'NR==1||$1!="clade"' > pilot_res.tsv
touch PILOT.DONE
exit 0
