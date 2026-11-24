"""Root a per-clade species tree on the clade's backbone outgroups (21_hog_select.OUTG).
Input:  work/tree_clade/<clade>/concat_bl.treefile (ASTRAL topology, concatenated-ML branch lengths)
        falls back to astral.tre if the branch-length step failed.
Output: work/tree_clade/<clade>/rooted.nwk (+ prints ASTRAL local posterior of every ingroup node < 0.95)
Usage: python 19b_root_clade_tree.py CLADE"""
import os, sys
from importlib import import_module
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv")
sys.path.insert(0, f"{R}/scripts")
OUTG = import_module("21_hog_select").OUTG

def main():
    c = sys.argv[1]; d = f"{R}/work/tree_clade/{c}"
    src = f"{d}/concat_bl.treefile" if os.path.exists(f"{d}/concat_bl.treefile") else f"{d}/astral.tre"
    t = Tree(src, format=1)
    og = [x for x in OUTG[c] if t.search_nodes(name=x)]
    if not og:
        sys.exit(f"{c}: no outgroup in tree")
    ing = [x for x in t.get_leaf_names() if x not in og][0]
    t.set_outgroup(t & ing)  # temporary root inside the ingroup so the outgroup MRCA is a clade
    if len(og) > 1 and t.check_monophyly(values=og, target_attr="name", unrooted=True)[0]:
        t.set_outgroup(t.get_common_ancestor(og))
    else:
        t.set_outgroup(t & og[0])  # first listed = most distant backbone taxon
        print(f"{c}: outgroups {og} not a clade; rooted on {og[0]}")
    t.write(outfile=f"{d}/rooted.nwk", format=1)
    # support values from ASTRAL (local posterior) for the record
    a = Tree(f"{d}/astral.tre", format=0) if os.path.exists(f"{d}/astral.tre") else None
    low = []
    if a:
        for n in a.traverse():
            if not n.is_leaf() and not n.is_root() and n.support < 0.95:
                low.append((round(n.support, 3), ",".join(sorted(n.get_leaf_names()))))
    print(f"{c}: rooted on {og} from {os.path.basename(src)}; {len(t)} tips; low-support nodes: {low}")

if __name__ == "__main__":
    main()
