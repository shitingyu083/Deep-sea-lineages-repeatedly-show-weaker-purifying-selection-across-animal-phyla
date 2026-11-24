#!/bin/bash
# R2 CAFE permutation nulls nullp11-40, pass 2 (PREREG_nature.md): Mollusca and Actinopteri_m, fitted on the
# rescued inputs (work/cafe/<C>_rescue: filtered counts, tree, error model of the successful base fit).
# Replaces the pass-2 branch of 65a (wrong dir names). Lambda trees: 49_cafe_prep.py --nperm 40 into scratch;
# deep and nullp1-10 must match the rescue copies byte for byte, else the clade is skipped.
B=$HOME/deepsea_conv; cd $B || exit 0
PY=envs/dc/bin/python; CAFE=envs/dc/bin/cafe5; O=work/addon/R2/cafe_perm; mkdir -p $O
while [ ! -e work/cafe/CAFE_RESCUE.DONE ] || [ ! -e $O/PASS1.DONE ]; do sleep 300; done
: > $O/jobs2.txt
for C in Mollusca Actinopteri_m; do
  D=work/cafe/${C}_rescue; P=$O/${C}_rescue; mkdir -p $P/prep $P/logs
  [ -e $D/CAFE.FAIL ] && { echo "$C rescue failed, skipped" >> $O/check2.txt; continue; }
  [ -s $P/prep/lambda_nullp40.nwk ] || $PY scripts/49_cafe_prep.py $C --out $P/prep --nperm 40 > $P/logs/prep.log 2>&1
  ok=1
  for f in lambda_deep.nwk $(cd $D && ls lambda_nullp*.nwk); do
    cmp -s $D/$f $P/prep/$f || { echo "$C $f differs" >> $O/check2.txt; ok=0; }
  done
  [ $ok = 1 ] || continue
  E=""; [ -s $D/base_e/Base_error_model.txt ] && E="-e$D/base_e/Base_error_model.txt"
  for k in $(seq 11 40); do
    echo "[ -s $P/nullp$k/Base_results.txt ] || nice -n 10 $CAFE -i $D/counts.tsv -t $D/tree.nwk -c 2 -o $P/nullp$k -y $P/prep/lambda_nullp$k.nwk $E > $P/logs/nullp$k.log 2>&1" >> $O/jobs2.txt
  done
done
echo "$(date +%F_%T) pass2 jobs $(wc -l < $O/jobs2.txt)" >> $O/run.log
cat $O/jobs2.txt | xargs -P 24 -I{} bash -c "{} < /dev/null"
echo "$(date +%F_%T) pass2 done" >> $O/run.log
touch $O/PASS2.DONE
exit 0
