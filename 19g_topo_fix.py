#!/usr/bin/env python3
"""Minimal, documented corrections of a species-tree topology for artefacts found by 19g_topo_qc.py
(long-branch attraction of single fast-evolving tips, a mis-rooted clade, an identity-conflicted tip).
Ops are read from SPEC (tab-separated: op, subject, target, reason, refs), applied in order:
  drop    CODE   -         remove the tip
  sister  CODE   GROUP     regraft the tip on the stem of GROUP (GROUP minus the tip must be a clade)
  reroot  GROUP  SUBGROUP  re-root the clade GROUP on its sub-clade SUBGROUP (both must be bipartitions)
GROUP = a HIGHER name of 19g_topo_qc.py or rank:value (phylum/class/order/family:value of meta/species_codes.tsv).
Every reference clade is checked to be monophyletic before the move, so nothing else is rearranged.
Usage: python scripts/19g_topo_fix.py TREE_IN OUTPREFIX OG1,OG2 [SPEC]
  TREE_IN is rooted on the outgroups first. SPEC defaults to meta/topo_fix.tsv.
Writes OUTPREFIX.tre (rooted topology), OUTPREFIX_unrooted.tre, OUTPREFIX.log.tsv (one line per op)."""
import csv, importlib.util, os, sys

R = os.path.expanduser("~/deepsea_conv")
spec = importlib.util.spec_from_file_location("topo_qc", f"{R}/scripts/19g_topo_qc.py")
qc = importlib.util.module_from_spec(spec); spec.loader.exec_module(qc)


def members(g, tax, tips):
    if g in qc.HIGHER:
        return {c for c in tips if c in tax and qc.HIGHER[g](tax[c])}
    rank, v = g.split(":", 1)
    return {c for c in tips if c in tax and tax[c][rank] == v}


def summ(leaves, tax):
    c = {}
    for l in leaves:
        k = tax[l]["phylum"] if l in tax else l
        c[k] = c.get(k, 0) + 1
    return ",".join(f"{k}:{v}" for k, v in sorted(c.items(), key=lambda x: -x[1]))


def detach_tip(t, code):
    l = t & code
    p = l.up
    l.detach()
    if p is not t and len(p.children) == 1:
        p.delete(prevent_nondicotomic=False)
    return l


def clade_of(t, ms, what):
    ms = set(ms)
    if not ms:
        raise SystemExit(f"{what}: no members in the tree")
    anc = t.get_common_ancestor(list(ms)) if len(ms) > 1 else t & next(iter(ms))
    got = set(anc.get_leaf_names())
    if got != ms:
        raise SystemExit(f"{what}: not a clade (intruders {sorted(got - ms)})")
    return anc


def main():
    tin, outp, og = sys.argv[1:4]
    fspec = sys.argv[4] if len(sys.argv) > 4 else f"{R}/meta/topo_fix.tsv"
    og = og.split(",")
    tax = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}
    t = qc.reroot(qc.load(tin), og)
    t0 = t.copy("newick")
    ops = [r for r in csv.DictReader(open(fspec), delimiter="\t") if r["op"] and not r["op"].startswith("#")]
    log = open(f"{outp}.log.tsv", "w")
    log.write("op\tsubject\ttarget\told_sister\tnew_sister\treason\trefs\n")
    for r in ops:
        op, s, g = r["op"], r["subject"], r["target"]
        tips = set(t.get_leaf_names())
        if op == "drop":
            if s not in tips:
                raise SystemExit(f"drop {s}: not in the tree")
            old = summ(set((t & s).up.get_leaf_names()) - {s}, tax)
            detach_tip(t, s)
            new = "-"
        elif op == "sister":
            if s not in tips:
                raise SystemExit(f"sister {s}: not in the tree")
            old = summ(set((t & s).up.get_leaf_names()) - {s}, tax)
            leaf = detach_tip(t, s)
            tgt = clade_of(t, members(g, tax, tips) - {s}, f"sister {s} {g}")
            p = tgt.up
            if p is None:
                raise SystemExit(f"sister {s} {g}: target is the root")
            tgt.detach()
            n = p.add_child()
            n.add_child(tgt); n.add_child(leaf)
            new = summ(tgt.get_leaf_names(), tax)
        elif op == "reroot":
            ms = members(s, tax, tips)
            cl = clade_of(t, ms, f"reroot {s}")
            sub = members(g, tax, tips) & ms
            p = cl.up
            old = summ(set(cl.children[0].get_leaf_names()), tax) + " | " + summ(set(cl.children[1].get_leaf_names()), tax)
            cl.detach()
            cl.set_outgroup([l for l in cl.get_leaves() if l.name not in sub][0])
            a = clade_of(cl, sub, f"reroot {s} on {g}")
            cl.set_outgroup(a)
            p.add_child(cl)
            new = summ(set(cl.children[0].get_leaf_names()), tax) + " | " + summ(set(cl.children[1].get_leaf_names()), tax)
        else:
            raise SystemExit(f"unknown op {op}")
        log.write(f"{op}\t{s}\t{g}\t{old}\t{new}\t{r['reason']}\t{r['refs']}\n")
        print(f"{op} {s} {g}: {old}  ->  {new}")
    for n in t.traverse():
        if not n.is_leaf() and len(n.children) != 2:
            raise SystemExit(f"non-binary node with {len(n.children)} children")
    t.write(outfile=f"{outp}.tre", format=9)
    u = t.copy("newick"); u.unroot()
    u.write(outfile=f"{outp}_unrooted.tre", format=9)
    keep = set(t.get_leaf_names())
    t0.prune([l for l in t0.get_leaf_names() if l in keep])
    rf = t.robinson_foulds(t0, unrooted_trees=True)
    print(f"{len(t)} tips; RF(fixed, input pruned to the same tips) = {rf[0]} / {rf[1]}")
    log.close()


if __name__ == "__main__":
    main()
