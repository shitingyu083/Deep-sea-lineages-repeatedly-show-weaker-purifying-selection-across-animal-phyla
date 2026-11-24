#!/usr/bin/env python3
"""Topology QC of a species tree against the taxonomy in meta/species_codes.tsv.
Tests monophyly of every phylum/class/order/family with >=2 sampled species and of the higher clades in
HIGHER (consensus metazoan relationships, refs in the Methods). For each group it finds the best-fitting
clade (the clade maximising members - non-members inside it) and lists the members outside it
(misplaced members) and the non-members inside it (intruders).
Usage: python scripts/19g_topo_qc.py TREE OUTPREFIX [OG1,OG2]
  With outgroups the tree is rooted on them first (otherwise it is used as rooted in the file).
Writes OUTPREFIX.monophyly.tsv; prints the non-monophyletic groups."""
import csv, os, sys
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv")
NONBIL = {"Porifera", "Cnidaria", "Ctenophora", "Placozoa"}
ECDY = {"Arthropoda", "Nematoda", "Priapulida", "Tardigrada", "Onychophora", "Nematomorpha", "Kinorhyncha"}
SPIR = {"Mollusca", "Annelida", "Nemertea", "Brachiopoda", "Bryozoa", "Phoronida", "Entoprocta",
        "Rotifera", "Platyhelminthes", "Gastrotricha", "Chaetognatha"}
DEUT = {"Chordata", "Echinodermata", "Hemichordata"}
TUNIC = {"Ascidiacea", "Thaliacea", "Appendicularia"}
CYCLO = {"Myxini", "Hyperoartia"}
HOLOSTEI = {"Lepisosteiformes", "Semionotiformes", "Amiiformes"}   # gars are filed as Semionotiformes in species_codes
NONTEL = {"Acipenseriformes", "Polyodontiformes"} | HOLOSTEI
ELOPO = {"Anguilliformes", "Notacanthiformes", "Elopiformes", "Albuliformes", "Saccopharyngiformes"}
OSTEOGL = {"Osteoglossiformes", "Hiodontiformes"}
OTOCEPH = {"Alepocephaliformes", "Clupeiformes", "Cypriniformes", "Siluriformes", "Characiformes", "Gonorynchiformes",
           "Gymnotiformes"}
NONACANTHO = NONTEL | ELOPO | OSTEOGL | OTOCEPH | {"Argentiniformes", "Osmeriformes", "Salmoniformes", "Esociformes",
              "Galaxiiformes", "Stomiiformes", "Aulopiformes", "Myctophiformes", "Ateleopodiformes"}
CHELIC = {"Merostomata", "Arachnida", "Pycnogonida"}
OCTO = {"Scleralcyonacea", "Malacalcyonacea", "Alcyonacea", "Pennatulacea", "Helioporacea"}
HEXA = {"Actiniaria", "Scleractinia", "Antipatharia", "Zoantharia", "Corallimorpharia"}
BATOID = {"Rajiformes", "Myliobatiformes", "Torpediniformes", "Rhinopristiformes"}
ACULIFERA = {"Polyplacophora", "Solenogastres", "Caudofoveata"}
MEDUSO = {"Hydrozoa", "Scyphozoa", "Cubozoa", "Staurozoa"}
HIGHER = {
    "Metazoa": lambda t: t["phylum"] != "",
    "Bilateria": lambda t: t["phylum"] not in NONBIL | {""},
    "Protostomia": lambda t: t["phylum"] in ECDY | SPIR,
    "Ecdysozoa": lambda t: t["phylum"] in ECDY,
    "Spiralia": lambda t: t["phylum"] in SPIR,
    "Deuterostomia": lambda t: t["phylum"] in DEUT,
    "Ambulacraria": lambda t: t["phylum"] in ("Echinodermata", "Hemichordata"),
    "Olfactores": lambda t: t["phylum"] == "Chordata" and t["class"] != "Leptocardii",
    "Tunicata": lambda t: t["class"] in TUNIC,
    "Vertebrata": lambda t: t["phylum"] == "Chordata" and t["class"] not in TUNIC | {"Leptocardii"},
    "Cyclostomata": lambda t: t["class"] in CYCLO,
    "Gnathostomata": lambda t: t["phylum"] == "Chordata" and t["class"] not in TUNIC | CYCLO | {"Leptocardii"},
    "Osteichthyes": lambda t: t["phylum"] == "Chordata" and t["class"] not in TUNIC | CYCLO | {"Leptocardii", "Chondrichthyes"},
    "Actinopterygii": lambda t: t["class"] in ("Actinopteri", "Cladistia"),
    "Neopterygii": lambda t: t["class"] == "Actinopteri" and t["order"] not in {"Acipenseriformes", "Polyodontiformes", ""},
    "Teleostei": lambda t: t["class"] == "Actinopteri" and t["order"] not in NONTEL | {""},
    "Medusozoa": lambda t: t["class"] in MEDUSO,
    "Aculifera": lambda t: t["phylum"] == "Mollusca" and t["class"] in ACULIFERA,
    "Conchifera": lambda t: t["phylum"] == "Mollusca" and t["class"] not in ACULIFERA | {""},
    "Lophotrochozoa": lambda t: t["phylum"] in SPIR - {"Rotifera", "Platyhelminthes", "Gastrotricha", "Chaetognatha"},
    "Pancrustacea": lambda t: t["phylum"] == "Arthropoda" and t["class"] not in CHELIC | {""},
    "Eleutherozoa": lambda t: t["phylum"] == "Echinodermata" and t["class"] not in ("Crinoidea", ""),
    "Echinozoa": lambda t: t["class"] in ("Echinoidea", "Holothuroidea"),
    "Silicea": lambda t: t["class"] in ("Demospongiae", "Hexactinellida"),
    "Octocorallia": lambda t: t["class"] == "Anthozoa" and t["order"] in OCTO,
    "Hexacorallia": lambda t: t["class"] == "Anthozoa" and t["order"] in HEXA,
    "Elasmobranchii": lambda t: t["class"] == "Chondrichthyes" and t["order"] not in ("Chimaeriformes", ""),
    "Batoidea": lambda t: t["order"] in BATOID,
    "Holostei": lambda t: t["order"] in HOLOSTEI,
    "Elopomorpha": lambda t: t["class"] == "Actinopteri" and t["order"] in ELOPO,
    "Clupeocephala": lambda t: t["class"] == "Actinopteri" and t["order"] not in NONTEL | ELOPO | OSTEOGL | {""},
    "Euteleostei": lambda t: t["class"] == "Actinopteri" and t["order"] not in NONTEL | ELOPO | OSTEOGL | OTOCEPH | {""},
    "Acanthomorpha": lambda t: t["class"] == "Actinopteri" and t["order"] not in NONACANTHO | {""},
}


def load(f):
    t = Tree(open(f).read().strip(), format=1, quoted_node_names=True)
    for n in t.traverse():
        if not n.is_leaf():
            n.name = ""
    return t


def reroot(t, og):
    og = [o for o in og if o in set(t.get_leaf_names())]
    if not og:
        return t
    t.set_outgroup([l for l in t.get_leaves() if l.name not in og][0])
    t.set_outgroup(t.get_common_ancestor(og) if len(og) > 1 else t & og[0])
    return t


def best_clade(t, m):
    """Clade maximising (#members - #non-members) inside it."""
    best, bs = None, None
    for n in t.traverse("postorder"):
        if n.is_leaf():
            n.k = int(n.name in m); n.s = 1
        else:
            n.k = sum(c.k for c in n.children); n.s = sum(c.s for c in n.children)
        sc = 2 * n.k - n.s
        if bs is None or sc > bs or (sc == bs and n.s < best.s):
            best, bs = n, sc
    return best


def main():
    tree, outp = sys.argv[1:3]
    og = sys.argv[3].split(",") if len(sys.argv) > 3 and sys.argv[3] != "-" else []
    tax = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}
    t = reroot(load(tree), og)
    tp = [l for l in t.get_leaf_names() if l in tax]
    groups = []
    for rank in ("phylum", "class", "order", "family"):
        vals = {}
        for c in tp:
            v = tax[c][rank]
            if v and v != "NA":
                vals.setdefault(v, []).append(c)
        groups += [(rank, v, m) for v, m in sorted(vals.items()) if len(m) >= 2]
    for name, f in HIGHER.items():
        m = [c for c in tp if f(tax[c])]
        if len(m) >= 2:
            groups.append(("higher", name, m))
    desc = lambda c: f"{c}:{tax[c]['phylum']}/{tax[c]['class']}/{tax[c]['order']}"
    bad = 0
    with open(f"{outp}.monophyly.tsv", "w") as o:
        o.write("rank\tgroup\tn\tmonophyletic\tcore_n\tmisplaced_members\tintruders\n")
        for rank, g, m in groups:
            ms = set(m)
            anc = t.get_common_ancestor(m)
            mono = set(anc.get_leaf_names()) <= ms
            core = best_clade(t, ms)
            inside = set(core.get_leaf_names())
            outside = sorted(ms - inside)
            intr = sorted(inside - ms)
            o.write(f"{rank}\t{g}\t{len(m)}\t{int(mono)}\t{core.k}\t"
                    f"{';'.join(desc(c) for c in outside)}\t"
                    f"{';'.join(desc(c) if c in tax else c for c in intr)}\n")
            if not mono:
                bad += 1
                print(f"[{rank}] {g} (n={len(m)}): outside={[desc(c) for c in outside]} "
                      f"intruders={[desc(c) if c in tax else c for c in intr]}")
    print(f"{len(groups)} groups tested, {bad} not monophyletic; tips {len(tp)} (+{len(t) - len(tp)} untaxed)")


if __name__ == "__main__":
    main()
