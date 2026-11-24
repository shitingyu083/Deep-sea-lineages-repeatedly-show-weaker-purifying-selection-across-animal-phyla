#!/bin/bash
# Extra gene-tree worker for a 19c_species_tree.sh run on W (2026-10-19). Started in its own session, so it
# gets its own CPU share under autogroup scheduling. Smallest loci first (the main run goes largest first);
# the per-locus flock in 19c_locus.sh keeps every locus to one run at a time. 19c_post.sh waits for all
# workers of W before TreeShrink.
# Usage: cd ~/deepsea_conv; setsid nohup bash scripts/19c_worker.sh W NJ TPJ > logs/x.log 2>&1 < /dev/null &
cd ~/deepsea_conv
W=$(readlink -f $1); NJ=${2:-32}; TPJ=${3:-3}
echo "$(date '+%F %T') worker start W=$W NJ=$NJ TPJ=$TPJ"
ls -Sr $W/loci | sed 's/\.faa$//' | xargs -P $NJ -I{} nice -n 10 bash scripts/19c_locus.sh $W {} $TPJ
echo "$(date '+%F %T') worker finished"
