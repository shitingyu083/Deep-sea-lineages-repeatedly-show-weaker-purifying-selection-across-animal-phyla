#!/usr/bin/env python3
"""Add supplementary species (deep arthropods / deep Aulopiformes + shallow matches)
to meta/species_codes.tsv and write meta/panel_v1b.tsv for 05_download.sh.
Depth for 'sp.' or OBIS-absent species is the genus-level OBIS summary; the basis is
recorded in meta/depth_overrides.tsv. Run from ~/deepsea_conv."""
import csv, shutil, time, statistics, collections

R = "."
ADD = [  # species, accession, source, role, phylum, class, order, family, depth_basis
    ("Bathynomus jamesi", "GCA_023014485.1", "miniprot", "deep", "Arthropoda", "Malacostraca", "Isopoda", "Cirolanidae", "genus"),
    ("Ceratothoa steindachneri", "GCA_963969495.1", "miniprot", "shallow_match", "Arthropoda", "Malacostraca", "Isopoda", "Cymothoidae", "species"),
    ("Idotea baltica", "GCA_023373965.1", "ncbi_annot", "shallow_match", "Arthropoda", "Malacostraca", "Isopoda", "Idoteidae", "literature"),
    ("Palinurus mauritanicus", "GCA_052075845.1", "miniprot", "deep", "Arthropoda", "Malacostraca", "Decapoda", "Palinuridae", "species"),
    ("Bathypterois sp.", "GCA_048564785.1", "miniprot", "deep", "Chordata", "Actinopteri", "Aulopiformes", "Ipnopidae", "genus"),
    ("Ipnops sp.", "GCA_048544905.1", "miniprot", "deep", "Chordata", "Actinopteri", "Aulopiformes", "Ipnopidae", "genus"),
    ("Parasudis fraserbrunneri", "GCA_900302295.1", "miniprot", "deep", "Chordata", "Actinopteri", "Aulopiformes", "Chlorophthalmidae", "species"),
    ("Harpadon nehereus", "GCA_051994445.1", "miniprot", "shallow_match", "Chordata", "Actinopteri", "Aulopiformes", "Synodontidae", "species"),
    ("Latropiscis purpurissatus", "GCA_988284235.1", "miniprot", "shallow_match", "Chordata", "Actinopteri", "Aulopiformes", "Aulopidae", "species"),
]
# Idotea baltica is absent from OBIS; intertidal/shallow subtidal (littoral algae beds)
LIT = {"Idotea baltica": (5.0, 0.0, "literature: littoral/sublittoral algal beds, <30 m")}

obis = {r["sp"]: r for r in csv.DictReader(open(f"{R}/meta/obis_depth_by_species.tsv"), delimiter="\t")}

def depth(sp, basis):
    if basis == "literature":
        return LIT[sp]
    if basis == "species":
        r = obis[sp]
        return float(r["depth_median"]), float(r["frac_gt200"]), f"OBIS species n_depth={r['n_depth']}"
    g = sp.split()[0]
    rows = [r for k, r in obis.items() if k.split()[0] == g and r["depth_median"] not in ("", "nan")
            and int(float(r["n_depth"])) >= 2]
    med = statistics.median(float(r["depth_median"]) for r in rows)
    frac = statistics.median(float(r["frac_gt200"]) for r in rows)
    return med, frac, f"OBIS genus {g}: median of {len(rows)} species medians"

codes_f = f"{R}/meta/species_codes.tsv"
shutil.copy(codes_f, f"{codes_f}.bak_{time.strftime('%Y%m%d_%H%M%S')}")
rows = list(csv.DictReader(open(codes_f), delimiter="\t"))
cols = list(rows[0].keys())
have_acc = {r["accession"] for r in rows}
have_code = {r["code"] for r in rows}
best = {r["accession"]: r for r in csv.DictReader(open(f"{R}/meta/species_best_assembly.tsv"), delimiter="\t")}

new, panel, ov = [], [], []
for sp, acc, src, role, ph, cl, od, fam, basis in ADD:
    if acc in have_acc:
        continue
    g, s = sp.split()[0], sp.split()[1].rstrip(".")
    code = (g[:3] + s[:3]).capitalize()
    i = 1
    while code in have_code:
        code = (g[:3] + s[:2]).capitalize() + str(i); i += 1
    have_code.add(code)
    med, frac, note = depth(sp, basis)
    new.append(dict(code=code, species=sp, accession=acc, source=src, role=role, phylum=ph, **{"class": cl},
                    order=od, family=fam, depth_median=round(med, 2), frac_gt200=round(frac, 3)))
    b = best[acc]
    panel.append([role, ph, cl, od, fam, sp, b["tax_id"], acc, b["annot"], b["level"],
                  b["contig_n50"], b["total_len"], "", med, frac])
    ov.append([code, sp, basis, round(med, 2), round(frac, 3), note])

with open(codes_f, "a") as f:
    w = csv.DictWriter(f, cols, delimiter="\t", lineterminator="\n")
    for r in new:
        w.writerow(r)
hdr = "role phylum class order family species tax_id accession annot level contig_n50 total_len n_depth depth_median frac_gt200".split()
with open(f"{R}/meta/panel_v1b.tsv", "w") as f:
    f.write("\t".join(hdr) + "\n")
    for p in panel:
        f.write("\t".join(map(str, p)) + "\n")
with open(f"{R}/meta/depth_overrides.tsv", "a") as f:
    for o in ov:
        f.write("\t".join(map(str, o)) + "\n")
for o in ov:
    print(*o, sep="\t")
print(collections.Counter(r["role"] for r in rows + new))
