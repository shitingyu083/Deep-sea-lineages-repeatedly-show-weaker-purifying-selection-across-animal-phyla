#!/bin/bash
# eggNOG-mapper 2.1.13 (eggNOG 5.0.2, DIAMOND) on every final proteome.
# Waits for db/eggnog/dmnd.done. Re-runs if proteome newer than annotations.
# Usage: nohup bash scripts/17_emapper_loop.sh 4 32 > logs/emapper.log 2>&1 &
set -u
R=~/deepsea_conv; export PATH=$R/envs/dc/bin:$PATH
NJ=${1:-4}; T=${2:-32}
D=$R/db/eggnog; O=$R/work/eggnog; mkdir -p $O
export R D O T
em_one() {
  c=$1; f=$R/data/proteomes/$c.faa; t=$O/$c.tmp; mkdir -p $t
  emapper.py -i $f --itype proteins -m diamond --data_dir $D --cpu $T \
    --output $c --output_dir $t --temp_dir $t --override > $t/run.log 2>&1 \
    && mv $t/$c.emapper.annotations $t/$c.emapper.seed_orthologs $O/ \
    && rm -rf $t || echo "EMFAIL $c" >> $O/failed.txt
}
export -f em_one
until [ -f $D/dmnd.done ]; do sleep 300; done
while true; do
  done_flag=0; [ -f $R/work/ANNOT.DONE ] && done_flag=1
  for f in $R/data/proteomes/*.faa; do
    c=$(basename $f .faa); a=$O/$c.emapper.annotations
    [ -n "$(find $f -mmin -3)" ] && continue
    grep -qx "EMFAIL $c" $O/failed.txt 2>/dev/null && continue
    { [ ! -s $a ] || [ $f -nt $a ]; } && echo $c
  done > $O/todo.txt
  if [ -s $O/todo.txt ]; then
    xargs -P $NJ -I{} bash -c 'em_one {}' < $O/todo.txt
  elif [ $done_flag = 1 ]; then
    break
  else
    sleep 300
  fi
done
echo DONE > $O/EMAPPER.DONE
