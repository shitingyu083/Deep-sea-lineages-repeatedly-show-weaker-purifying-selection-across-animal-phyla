#!/usr/bin/env python3
"""Master species table: panel_v1 + lab-local molluscs, one row per species.

source = local       : lab proteome/CDS (06_prep_local.py)
         ncbi_annot  : NCBI RefSeq/GenBank annotation (08_prep_ncbi_annot.py)
         miniprot    : unannotated assembly -> homology projection (10_*)
Output meta/species_codes.tsv
"""
import csv, os
R = os.path.expanduser("~/deepsea_conv")

LOCAL_CODES = {  # species -> code for lab proteomes kept as-is
    "Archivesica marissinica": "Ama", "Bathyacmaea lactea": "Bla",
    "Gigantidas platifrons": "Bpl", "Conchocele bisecta": "Cbi",
    "Chrysomallon squamiferum": "Csq", "Xylophaga dorsalis": "Xdo",
    "Hanleya hanleyi": "Hha", "Mytilus coruscus": "Mco", "Rapana venosa": "Rve",
    "Sinonovacula constricta": "Sco", "Sepia pharaonis": "Sph", "Teredo navalis": "Tna",
}
EXTRA = [  # local species not in panel_v1: role, phylum, class, order, family, species, acc
    ("deep", "Mollusca", "Bivalvia", "Venerida", "Vesicomyidae", "Archivesica marissinica", "GCA_014843695.1"),
    ("deep", "Mollusca", "Gastropoda", "Patellogastropoda", "Pectinodontidae", "Bathyacmaea lactea", "GCA_053477725.1"),
    ("deep", "Mollusca", "Bivalvia", "Mytilida", "Mytilidae", "Gigantidas platifrons", "GCA_040955395.1"),
    ("deep", "Mollusca", "Bivalvia", "Myida", "Xylophagaidae", "Xylophaga dorsalis", "GCA_965225065.2"),
    ("deep", "Mollusca", "Cephalopoda", "Oegopsida", "Architeuthidae", "Architeuthis dux", "GCA_006491835.1"),
    ("shallow", "Mollusca", "Polyplacophora", "Lepidopleurida", "Hanleyidae", "Hanleya hanleyi", "GCA_036873755.1"),
    ("shallow", "Mollusca", "Gastropoda", "Lepetellida", "Haliotidae", "Haliotis rubra", "GCF_003918875.1"),
    ("shallow", "Mollusca", "Gastropoda", "Patellogastropoda", "Lottiidae", "Lottia gigantea", "GCF_000327385.1"),
    ("shallow", "Mollusca", "Bivalvia", "Mytilida", "Mytilidae", "Mytilus coruscus", "GCA_011752425.2"),
    ("shallow", "Mollusca", "Gastropoda", "Neogastropoda", "Muricidae", "Rapana venosa", "GCA_028751875.1"),
    ("shallow", "Mollusca", "Bivalvia", "Adapedonta", "Pharidae", "Sinonovacula constricta", "GCA_007844125.1"),
    ("shallow", "Mollusca", "Cephalopoda", "Sepiida", "Sepiidae", "Sepia pharaonis", ""),
    ("shallow", "Mollusca", "Bivalvia", "Myida", "Teredinidae", "Teredo navalis", ""),
]
FIXED = {"Haliotis rubra": "Hru", "Lottia gigantea": "Lgi", "Architeuthis dux": "Adu",
         "Gigantopelta aegis": "Gae", "Mercenaria mercenaria": "Mme"}


def main():
    rows = list(csv.DictReader(open(f"{R}/meta/panel_v1.tsv"), delimiter="\t"))
    have = {r["species"] for r in rows}
    for role, ph, cl, od, fa, sp, acc in EXTRA:
        if sp in have:
            continue
        rows.append(dict(role=role, phylum=ph, **{"class": cl}, order=od, family=fa,
                         species=sp, accession=acc, annot="none"))
    obis = {}
    for r in csv.DictReader(open(f"{R}/meta/obis_depth_by_species.tsv"), delimiter="	"):
        obis[r["sp"]] = r
    for r in rows:
        o = obis.get(r["species"])
        if o and not r.get("depth_median"):
            r["depth_median"], r["frac_gt200"] = o["depth_median"], o["frac_gt200"]
            r["n_depth"] = o["n_depth"]
    used = set(LOCAL_CODES.values()) | set(FIXED.values())
    out = []
    for r in rows:
        sp = r["species"]
        if sp in LOCAL_CODES:
            code, src = LOCAL_CODES[sp], "local"
        else:
            code = FIXED.get(sp)
            if code is None:
                g, s = (sp.split() + ["sp"])[:2]
                base = (g[:3] + s[:3]).capitalize()
                code, k = base, 1
                while code in used:
                    k += 1; code = f"{base}{k}"
                used.add(code)
            ann = r.get("annot", "none")
            if sp in ("Haliotis rubra", "Lottia gigantea"):
                ann = "refseq"
            src = "ncbi_annot" if ann in ("refseq", "genbank") else "miniprot"
        out.append([code, sp, r.get("accession", ""), src, r["role"], r["phylum"], r["class"],
                    r["order"], r["family"], r.get("depth_median", ""), r.get("frac_gt200", "")])
    with open(f"{R}/meta/species_codes.tsv", "w") as f:
        f.write("code\tspecies\taccession\tsource\trole\tphylum\tclass\torder\tfamily\tdepth_median\tfrac_gt200\n")
        for o in out:
            f.write("\t".join(map(str, o)) + "\n")
    from collections import Counter
    print(len(out), Counter(o[3] for o in out), Counter(o[4] for o in out))


if __name__ == "__main__":
    main()
