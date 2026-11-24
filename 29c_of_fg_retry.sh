#!/bin/bash
# Resume a clade's failed v2 OrthoFinder -b pass from its orthogroups (-fg) instead of redoing the search.
# 2026-10-20: Cnidaria's add_v2 (29_of_add_v2_loop.sh) finished the all-vs-all DIAMOND search and MCL
# (Orthogroups/ written 01:27) and then died in "Calculating gene distances" with OSError [Errno 24] Too many
# open files (soft nofile limit 1024; multiprocessing heap arenas). This reruns only the orthologue workflow
# (gene trees, STAG species tree, STRIDE root, HOGs) with the soft limit raised. Results go beside the -b
# results (<v1 WorkingDirectory>/OrthoFinder/Results_v2fg), where 21_hog_select.py finds them (newest N0.tsv).
# Usage: bash scripts/29c_of_fg_retry.sh CLADE [THREADS=48] [ALIGN_THREADS=8]
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
c=$1; T=${2:-48}; A=${3:-8}
d=work/orthofinder/$c
ulimit -n 65536 || { echo "$c: cannot raise the nofile limit"; exit 1; }
# the -b results that have orthogroups but no HOGs
rd=$(ls -d $d/of/Results_*/WorkingDirectory/OrthoFinder/Results_* 2>/dev/null | while read r; do
       [ -s $r/Orthogroups/Orthogroups.tsv ] && [ ! -e $r/Phylogenetic_Hierarchical_Orthogroups/N0.tsv ] && echo $r
     done | tail -1)
[ -n "$rd" ] || { echo "$c: no -b results with orthogroups and without HOGs"; exit 1; }
rm -f $d/.v2_fail $d/.v2_done
touch $d/.v2_fg_started $d/.v2_try_stamp
echo "$(date '+%F %T') $c: orthofinder -fg $rd -t $T -a $A -X -n v2fg (nofile $(ulimit -n))"
nice -n 10 orthofinder -fg $rd -t $T -a $A -X -n v2fg > $d/of_fg_v2.log 2>&1; rc=$?
n0=$(find $(dirname $rd) -path "*Phylogenetic_Hierarchical_Orthogroups/N0.tsv" -newer $d/.v2_try_stamp 2>/dev/null | head -1)
if [ $rc = 0 ] && [ -n "$n0" ]; then
  touch $d/.v2_done; echo "$(date '+%F %T') $c v2 done: $n0"; exit 0
fi
touch $d/.v2_fail; echo "$(date '+%F %T') $c -fg failed rc=$rc: $(grep -m1 -i error $d/of_fg_v2.log)"; exit 1
