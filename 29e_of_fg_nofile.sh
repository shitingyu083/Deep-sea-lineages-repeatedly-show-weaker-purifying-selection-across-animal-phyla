#!/bin/bash
# Resume a clade's failed plain OrthoFinder -f run from its orthogroups (-fg) with the soft nofile limit raised.
# 2026-10-21: Mollusca (78 spp) finished DIAMOND + MCL (Orthogroups/ 01:05) and died in "Calculating gene
# distances" with OSError [Errno 24] Too many open files (soft nofile 1024, multiprocessing heap arenas);
# the 20b rescue (-fg without ulimit) failed the same way (.rescue_fail kept so 20b does not race this).
# Same cause/fix as Cnidaria (29c). On success: .fail -> .fail_resumed, touch .done (43_autochain -> prep).
# Usage: bash scripts/29e_of_fg_nofile.sh CLADE [THREADS=48] [ALIGN_THREADS=8]
cd ~/deepsea_conv
export PATH=$PWD/envs/dc/bin:$PATH
c=$1; T=${2:-48}; A=${3:-8}
d=work/orthofinder/$c
ulimit -n 65536 || { echo "$c: cannot raise the nofile limit"; exit 1; }
res=$(ls -d $d/of/Results_* 2>/dev/null | while read r; do
        [ -s $r/Orthogroups/Orthogroups.tsv ] && [ ! -e $r/Phylogenetic_Hierarchical_Orthogroups/N0.tsv ] && echo $r
      done | tail -1)
[ -n "$res" ] || { echo "$c: no results with orthogroups and without HOGs"; exit 1; }
touch $d/.fg_nofile_started $d/.fg_try_stamp
echo "$(date '+%F %T') $c: orthofinder -fg $res -t $T -a $A -X -n fgnf (nofile $(ulimit -n))"
nice -n 10 orthofinder -fg $res -t $T -a $A -X -n fgnf > $d/of_fg_nofile.log 2>&1; rc=$?
n0=$(find $d/of -path "*Phylogenetic_Hierarchical_Orthogroups/N0.tsv" -newer $d/.fg_try_stamp 2>/dev/null | head -1)
if [ $rc = 0 ] && [ -n "$n0" ]; then
  [ -e $d/.fail ] && mv $d/.fail $d/.fail_resumed
  touch $d/.done; echo "$(date '+%F %T') $c done: $n0"; exit 0
fi
echo "$(date '+%F %T') $c -fg failed rc=$rc: $(grep -m1 -i error $d/of_fg_nofile.log)"; exit 1
