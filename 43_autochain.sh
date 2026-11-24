#!/bin/bash
# Auto-chain for the steps that wait on other runs (2026-10-20). Polls every 5 min; each step runs once
# (marker work/chain43/<step>.launched); launches go to their own sessions. Log logs/chain43.log.
# Touches only my own sessions, identified by PID + uid + cmdline (never by name patterns alone).
#   sel_Echinodermata, sel_Annelida  PREP.DONE -> 43b_sel_launch.sh CLADE work/draft_v2/species_tree.nwk 48 10
#        (both clades are complete in the draft tree; HOGs whose labelled topology changes in the final
#        tree are rerun after BL.DONE, decided by hand)
#   up_gt          work/tree_v3/GENETREES.DONE -> fish sub77 session (30_fish_sub_of.sh) and the 27 prep
#                  session (247553) to autogroup nice 0; fish xargs +32 slots (SIGUSR1, only if caught)
#   prep_Actinopteri  sub_v2/SUB.DONE -> 27_sel_prep_v2.sh 16 Actinopteri (autogroup then set to 0)
#   prep_Mollusca     work/orthofinder/Mollusca/.done -> 27_sel_prep_v2.sh 16 Mollusca (same)
#   up_bl          work/final_v3/BL.DONE -> running draft-tree selection sessions to autogroup nice 0
#   sel_<clade>    BL.DONE + PREP.DONE -> 43b_sel_launch.sh CLADE work/final_v3/species_tree.nwk 48 0
#                  (Arthropoda Chondrichthyes Porifera Cnidaria Mollusca Actinopteri)
# Exits when every step has run. Usage: (setsid nohup bash scripts/43_autochain.sh ... &)
R=~/deepsea_conv; cd $R
D=work/chain43; mkdir -p $D; L=logs/chain43.log
echo $$ > $D/chain.pid
F=work/final_v3; DRAFT=work/draft_v2/species_tree.nwk
FINAL_CLADES="Arthropoda Chondrichthyes Porifera Cnidaria Mollusca Actinopteri"
log() { echo "$(date +%F_%T) $*" >> $L; }
did() { [ -e $D/$1.launched ]; }
mark() { touch $D/$1.launched; }
# mine PID STRING: PID alive, owned by me, cmdline contains STRING
mine() { [ -d /proc/$1 ] && [ "$(stat -c %u /proc/$1 2>/dev/null)" = "$(id -u)" ] \
  && tr '\0' ' ' 2>/dev/null < /proc/$1/cmdline | grep -qF -- "$2"; }
# agset PID NICE STRING: set the autogroup nice of PID's session (verified to be mine)
agset() { local i
  mine $1 "$3" || { log "agset: $1 is not '$3' (or not mine): unchanged"; return 1; }
  for i in 1 2 3; do echo $2 > /proc/$1/autogroup 2>/dev/null && break; sleep 2; done
  log "autogroup $1 ($3): $(cat /proc/$1/autogroup 2>/dev/null | tr -s ' ')"; }
# pid of my process whose cmdline contains STRING and that leads its own session
leader() { local p c
  for p in /proc/[0-9]*; do p=${p#/proc/}
    [ "$(stat -c %u /proc/$p 2>/dev/null)" = "$(id -u)" ] || continue
    c=$(tr '\0' ' ' 2>/dev/null < /proc/$p/cmdline); case "$c" in *"$1"*) ;; *) continue;; esac
    [ "$(ps -o sid= -p $p 2>/dev/null | tr -d ' ')" = "$p" ] && { echo $p; return 0; }
  done; return 1; }
launch() { local n=$1; shift
  ( setsid nohup "$@" > logs/chain43_$n.out 2>&1 < /dev/null 7>&- 8>&- & )
  mark $n; log "launched $n: $*"; }
prep() { local c=$1 p
  launch prep_$c bash scripts/27_sel_prep_v2.sh 16 $c
  sleep 30; p=$(leader "27_sel_prep_v2.sh 16 $c") && agset $p 0 "27_sel_prep_v2.sh 16 $c"; }

log "start (pid $$)"
while :; do
  for c in Echinodermata Annelida; do
    did sel_$c || { [ -e work/codon/$c/PREP.DONE ] && launch sel_$c bash scripts/43b_sel_launch.sh $c $DRAFT 48 10; }
  done
  if ! did up_gt && [ -e work/tree_v3/GENETREES.DONE ]; then
    mark up_gt; s=3654866
    if agset $s 0 30_fish_sub_of.sh; then
      for x in $(ps -o pid=,comm= -s $s | awk '$2 == "xargs" {print $1}'); do
        mine $x xargs || continue
        cg=$(awk '/^SigCgt/ {print $2}' /proc/$x/status)
        if [ $(( 0x$cg & 0x200 )) != 0 ]; then
          for i in $(seq 32); do kill -USR1 $x; sleep 0.3; done; log "fish xargs $x: +32 slots (SigCgt $cg)"
        else log "fish xargs $x does not catch SIGUSR1 (SigCgt $cg): slots unchanged"; fi
      done
    fi
    agset 247553 0 27_sel_prep_v2.sh
  fi
  did prep_Actinopteri || { [ -e work/orthofinder/Actinopteri/sub_v2/SUB.DONE ] && prep Actinopteri; }
  did prep_Mollusca || { [ -e work/orthofinder/Mollusca/.done ] && prep Mollusca; }
  if [ -e work/orthofinder/Mollusca/.fail ] && ! did fail_Mollusca; then mark fail_Mollusca; log "Mollusca OrthoFinder .fail"; fi
  if [ -e work/orthofinder/Actinopteri/sub_v2/SUB.FAIL ] && ! did fail_fish; then mark fail_fish; log "fish SUB.FAIL"; fi
  if [ -e $F/BL.DONE ] && [ -s $F/species_tree.nwk ]; then
    if ! did up_bl; then
      mark up_bl
      for c in Echinodermata Annelida; do
        [ -s $D/sel_$c.pid ] && agset $(cat $D/sel_$c.pid) 0 43b_sel_launch.sh
      done
    fi
    for c in $FINAL_CLADES; do
      did sel_$c || { [ -e work/codon/$c/PREP.DONE ] && launch sel_$c bash scripts/43b_sel_launch.sh $c $F/species_tree.nwk 48 0; }
    done
  fi
  left=0
  for s in sel_Echinodermata sel_Annelida up_gt prep_Actinopteri prep_Mollusca up_bl; do did $s || left=$((left+1)); done
  for c in $FINAL_CLADES; do did sel_$c || left=$((left+1)); done
  [ $left = 0 ] && { log "all steps launched: exit"; exit 0; }
  sleep 300
done
