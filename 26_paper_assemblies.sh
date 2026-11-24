#!/bin/bash
# Resolve BioProject / assembly accessions cited in deep-sea genome papers (lit/epmc) to NCBI assemblies.
# Usage: bash scripts/26_paper_assemblies.sh meta/paper_accs.txt meta/paper_assemblies.jsonl
R=~/deepsea_conv; D=$R/envs/dc/bin/datasets; IN=$1; OUT=$2; mkdir -p $R/meta/cache_pa
one() { a=$1; f=$R/meta/cache_pa/$a.jsonl; [ -e $f ] && return 0
  timeout 120 $D summary genome accession $a --as-json-lines > $f.tmp 2>/dev/null; mv $f.tmp $f; }
export -f one; export R D
xargs -P 8 -I{} bash -c 'one {}' < $IN
cat $R/meta/cache_pa/*.jsonl > $OUT; wc -l $OUT
