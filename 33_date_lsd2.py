#!/usr/bin/env python3
"""Fossil-calibrated time tree with IQ-TREE LSD2 on a fixed rooted topology.
Calibrations: refs/calibrations_benton2015.tsv (Benton et al. 2015 Palaeontol. Electron. 18.1.1FC).
Each Benton crown node is defined by the split between two sides (taxonomy filters below, following
the node definitions in that file). A calibration is applied to MRCA(sideA + sideB) only when both
sides are sampled in the tree, so the sampled MRCA is the defined crown node (never a younger node).
Every calibration is checked on the tree rooted on the outgroups: tips under MRCA(sideA + sideB) that
belong to neither side ('intruders': a misplaced tip, or a taxonomy-filter gap) are listed in
calibrations_used.tsv with the composition of the MRCA's two child clades. LSD_INTRUDERS=skip (default)
leaves a calibration with intruders out; LSD_INTRUDERS=apply keeps it.
LSD2 takes hard bounds: interval '-soft_max:-min' (IQ-TREE date-file syntax) in Ma (tips = 0). Soft maxima are treated as
hard here; mcmctree (soft bounds) is the robustness check.
Usage: [LSD_BLFIX=1] [LSD_MEM=64G] [LSD_DRYRUN=1] python 33_date_lsd2.py TREE PARTITION_NEX|ALIGNMENT OUTDIR [OUTGROUP1,OUTGROUP2] [NTHREADS]
  With outgroups the tree is rooted on them and they are dropped from the time tree (--date-no-outgroup).
  LSD_BLFIX=1 keeps the input tree's branch lengths (use when TREE already carries ML lengths for PARTITION_NEX).
  A second argument not ending in .nex/.nexus is read as one concatenated alignment under LSD_MODEL
  (default LG+G4; the model the branch lengths were estimated with). LSD_MEM caps IQ-TREE's memory; it is
  only passed in this single-alignment mode because iqtree3 rejects -mem with partition models.
  LSD_DRYRUN=1 writes dates.txt and calibrations_used.tsv and stops before IQ-TREE.
Writes OUTDIR/{calibrations_used.tsv, dates.txt, lsd.timetree.nwk ...}."""
import csv, os, subprocess, sys
from ete3 import Tree

R = os.path.expanduser("~/deepsea_conv")
PROTO = {"Mollusca", "Annelida", "Brachiopoda", "Bryozoa", "Nemertea", "Rotifera", "Arthropoda",
         "Nematoda", "Priapulida", "Tardigrada", "Onychophora", "Platyhelminthes"}
LOPHO = {"Annelida", "Brachiopoda", "Bryozoa", "Nemertea"}       # sister side to Mollusca
ECDY = {"Arthropoda", "Nematoda", "Priapulida", "Tardigrada", "Onychophora"}
SPIR = PROTO - ECDY                                               # Benton's Lophotrochozoa includes rotifers and flatworms
DEUT = {"Chordata", "Echinodermata", "Hemichordata"}
ACULIFERA = {"Polyplacophora", "Solenogastres", "Caudofoveata"}
ERRANTIA = {"Phyllodocida", "Eunicida"}                           # Aciculata
SEDENTARIA = {"Sabellida", "Spionida", "Terebellida", "Capitellida", "Opheliida", "Scolecida"}
PALAEO_FAM = {"Oweniidae", "Magelonidae", "Chaetopteridae", "Amphinomidae"}  # branch outside Errantia+Sedentaria (refs in meta/topo_fix.tsv)
TUNIC = {"Ascidiacea", "Thaliacea", "Appendicularia"}
CYCLO = {"Myxini", "Hyperoartia"}
SARCO_CL = {"Sarcopterygii", "Dipnoi", "Actinistia", "Amphibia", "Reptilia", "Aves", "Mammalia"}
SARCO_OR = {"Ceratodontiformes", "Lepidosireniformes", "Coelacanthiformes"}  # lungfishes, coelacanths (class left blank in species_codes)
HOLOSTEI = {"Amiiformes", "Lepisosteiformes", "Semionotiformes"}  # gars are filed as Semionotiformes in species_codes
NONTEL = {"Polypteriformes", "Acipenseriformes", "Polyodontiformes"} | HOLOSTEI
ELOPO = {"Anguilliformes", "Notacanthiformes", "Elopiformes", "Albuliformes", "Saccopharyngiformes"}
OSTEOGL = {"Osteoglossiformes", "Hiodontiformes"}
OTOCEPH = {"Alepocephaliformes", "Clupeiformes", "Cypriniformes", "Siluriformes", "Characiformes",
           "Gonorynchiformes", "Gymnotiformes"}                   # Alepocephaliformes = Otomorpha (Betancur-R 2017)
NONACANTHO = NONTEL | ELOPO | OSTEOGL | OTOCEPH | {"Argentiniformes", "Osmeriformes", "Salmoniformes", "Esociformes",
              "Galaxiiformes", "Stomiiformes", "Aulopiformes", "Myctophiformes", "Ateleopodiformes"}
LAMPRI = {"Lampriformes", "Lampridiformes"}

P = lambda *ph: (lambda t: t["phylum"] in ph)
C = lambda *cl: (lambda t: t["class"] in cl)
O = lambda *od: (lambda t: t["order"] in od)
sarco = lambda t: t["phylum"] == "Chordata" and (t["class"] in SARCO_CL or t["order"] in SARCO_OR)
actinopterygii = lambda t: t["class"] in ("Actinopteri", "Cladistia")
osteich = lambda t: actinopterygii(t) or sarco(t)
gnatho = lambda t: t["class"] == "Chondrichthyes" or osteich(t)
teleost = lambda t: t["class"] == "Actinopteri" and t["order"] not in NONTEL | {""}
acantho = lambda t: t["class"] == "Actinopteri" and t["order"] not in NONACANTHO | {""}
CALS = {  # Benton node_no: (side A, side B)
    1: (P("Porifera"), lambda t: t["phylum"] not in ("Porifera", "")),
    2: (P("Cnidaria"), lambda t: t["phylum"] in PROTO | DEUT),
    3: (C("Anthozoa"), C("Hydrozoa", "Scyphozoa", "Cubozoa", "Staurozoa")),
    4: (lambda t: t["phylum"] in PROTO, lambda t: t["phylum"] in DEUT),
    5: (lambda t: t["phylum"] in SPIR, lambda t: t["phylum"] in ECDY),
    6: (P("Arthropoda"), lambda t: t["phylum"] in ECDY - {"Arthropoda"}),
    8: (C("Merostomata", "Arachnida", "Pycnogonida"), lambda t: t["phylum"] == "Arthropoda" and t["class"] not in ("Merostomata", "Arachnida", "Pycnogonida")),
    9: (P("Mollusca"), lambda t: t["phylum"] in LOPHO),
    10: (lambda t: t["phylum"] == "Mollusca" and t["class"] in ACULIFERA, lambda t: t["phylum"] == "Mollusca" and t["class"] not in ACULIFERA),
    11: (lambda t: t["phylum"] == "Annelida" and t["order"] in ERRANTIA and t["family"] not in PALAEO_FAM,
         lambda t: t["phylum"] == "Annelida" and t["order"] in SEDENTARIA and t["family"] not in PALAEO_FAM),
    12: (P("Chordata"), P("Echinodermata", "Hemichordata")),
    13: (P("Echinodermata"), P("Hemichordata")),
    14: (lambda t: t["phylum"] == "Echinodermata" and t["class"] == "Crinoidea", lambda t: t["phylum"] == "Echinodermata" and t["class"] != "Crinoidea"),
    16: (C("Leptocardii"), lambda t: t["phylum"] == "Chordata" and t["class"] != "Leptocardii"),
    17: (lambda t: t["class"] in TUNIC, lambda t: t["phylum"] == "Chordata" and t["class"] not in TUNIC | {"Leptocardii"}),
    18: (lambda t: t["class"] in CYCLO, gnatho),
    19: (C("Myxini"), C("Hyperoartia")),
    20: (C("Chondrichthyes"), osteich),
    21: (lambda t: t["class"] == "Chondrichthyes" and t["order"] == "Chimaeriformes", lambda t: t["class"] == "Chondrichthyes" and t["order"] not in ("Chimaeriformes", "")),
    22: (actinopterygii, sarco),
    25: (C("Cladistia"), C("Actinopteri")),
    27: (lambda t: t["class"] == "Actinopteri" and t["order"] in HOLOSTEI, teleost),
    28: (O("Amiiformes"), O("Lepisosteiformes", "Semionotiformes")),
    30: (lambda t: teleost(t) and t["order"] in ELOPO, lambda t: teleost(t) and t["order"] not in ELOPO),
    31: (lambda t: teleost(t) and t["order"] in OTOCEPH, lambda t: teleost(t) and t["order"] not in OTOCEPH | ELOPO | OSTEOGL),
    34: (lambda t: acantho(t) and t["order"] in LAMPRI, lambda t: acantho(t) and t["order"] not in LAMPRI),
}

def rooted(fn, og):
    t = Tree(open(fn).read().strip(), format=1, quoted_node_names=True)
    if og:
        t.set_outgroup(next(l for l in t if l.name not in og))
        t.set_outgroup(t & og[0] if len(og) == 1 else t.get_common_ancestor(og))
    return t

def check(t, a, b, tax):
    """intruders under MRCA(a + b) and the side composition of its child clades"""
    m = t.get_common_ancestor(a + b)
    sa, sb = set(a), set(b)
    intr = sorted(set(m.get_leaf_names()) - sa - sb)
    kids = "|".join(f"A{len(sa & set(c.get_leaf_names()))}B{len(sb & set(c.get_leaf_names()))}" for c in m.children)
    desc = [f"{x}({tax[x]['phylum']}/{tax[x]['class']}/{tax[x]['order']})" if x in tax else f"{x}(outgroup)" for x in intr]
    return intr, desc, kids

def main():
    tree, part, out = sys.argv[1:4]
    og = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] != "-" else ""
    nt = sys.argv[5] if len(sys.argv) > 5 else "8"
    policy = os.environ.get("LSD_INTRUDERS", "skip")
    if policy not in ("skip", "apply"):
        sys.exit("LSD_INTRUDERS must be skip or apply")
    os.makedirs(out, exist_ok=True)
    tax = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t")}
    cal = {int(r["node_no"]): r for r in csv.DictReader(open(f"{R}/refs/calibrations_benton2015.tsv"), delimiter="\t")}
    alltp = Tree(open(tree).read().strip(), format=1, quoted_node_names=True).get_leaf_names()
    tp = [t for t in alltp if t in tax]
    ogl = [o for o in og.split(",") if o in alltp] if og else []   # keep only the outgroups sampled in this tree
    og = ",".join(ogl)
    t = rooted(tree, ogl)
    rows, napp = [], 0
    with open(f"{out}/dates.txt", "w") as fd:
        for n, (fa, fb) in CALS.items():
            a = [x for x in tp if fa(tax[x])]; b = [x for x in tp if fb(tax[x])]
            if not a or not b:
                continue
            lo, hi = float(cal[n]["min_Ma"]), float(cal[n]["soft_max_Ma"])
            intr, desc, kids = check(t, a, b, tax)
            ok = not intr or policy == "apply"
            if intr:
                print(f"WARN node {n} {cal[n]['node']}: {len(intr)} intruder(s) under MRCA ({kids}): "
                      + ", ".join(desc[:10]) + (" ..." if len(desc) > 10 else "")
                      + ("; applied anyway" if ok else "; calibration skipped"), flush=True)
            if ok:
                fd.write(",".join(sorted(set(a + b))) + f" -{hi}:-{lo}\n")
                napp += 1
            rows.append((n, cal[n]["node"], len(a), len(b), lo, hi, len(intr), ",".join(desc) or "-", kids,
                         "yes" if ok else "no:intruders"))
    with open(f"{out}/calibrations_used.tsv", "w") as f:
        f.write("node_no\tnode\tn_sideA\tn_sideB\tmin_Ma\tmax_Ma\tn_intruders\tintruders\tmrca_children\tapplied\tsource\n")
        for u in rows:
            f.write("\t".join(map(str, u)) + f"\tBenton et al. 2015 PE 18.1.1FC node {u[0]}\n")
    print(f"{napp} of {len(rows)} sampled calibrations applied:",
          ", ".join(f"{u[0]}:{u[1]}" for u in rows if u[-1] == "yes"), flush=True)
    if not napp:
        sys.exit("no calibration applicable")
    if os.environ.get("LSD_DRYRUN") == "1":
        return
    if part.endswith((".nex", ".nexus")):
        data = ["-p", part]
    else:
        data = ["-s", part, "-m", os.environ.get("LSD_MODEL", "LG+G4")]
        if os.environ.get("LSD_MEM"):
            data += ["-mem", os.environ["LSD_MEM"]]
    cmd = [f"{R}/envs/dc/bin/iqtree"] + data + ["-te", tree, "--date", f"{out}/dates.txt",
           "--date-tip", "0", "--date-ci", "100",
           "-T", nt, "--prefix", f"{out}/lsd", "-redo", "-quiet"]
    if og:
        cmd += ["-o", og, "--date-no-outgroup"]  # root on the outgroups, then drop them from the time tree
    if os.environ.get("LSD_BLFIX") == "1":
        cmd += ["-blfix"]                      # keep the input tree's branch lengths (no re-optimisation)
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
