#!/bin/bash
# R3 (PREREG_nature.md): RELAX and BUSTED-PH on the RTF1 HOG of each clade without HyPhy results (real + null1).
B=$HOME/deepsea_conv; HY=$B/envs/dc/bin/hyphy; O=$B/work/addon/R3/rtf1; mkdir -p $O; cd $O || exit 0
LIST="Actinopteri_m:N0.HOG0011415 Annelida:N0.HOG0010811 Arthropoda:N0.HOG0011086 Echinodermata:N0.HOG0010916 Cnidaria:N0.HOG0010137 Chondrichthyes:N0.HOG0007318 Mollusca:N0.HOG0014485"
: > jobs.txt
for x in $LIST; do C=${x%%:*}; h=${x##*:}; A=$B/work/codon/${C%_m}/aln/$h.fas; T=$B/work/codon/$C/trees
  [ -f $A ] || { echo "$C $h no aln" >> missing.txt; continue; }
  for s in "" .null1; do tr=$T/$h$s.hyphy.nwk; [ -f $tr ] || { echo "$C $h$s no tree" >> missing.txt; continue; }
    p=$O/$C.$h${s:-.real}
    echo "[ -s $p.relax.json ] || $HY CPU=1 relax --alignment $A --tree $tr --test FG --reference BG --models Minimal --output $p.relax.json > $p.relax.log 2>&1" >> jobs.txt
    echo "[ -s $p.bph.json ] || $HY CPU=1 busted-ph --alignment $A --tree $tr --branches FG --srv No --output $p.bph.json > $p.bph.log 2>&1" >> jobs.txt
  done; done
cat jobs.txt | xargs -P 14 -I{} bash -c "{} < /dev/null"
touch RTF1.DONE
exit 0
