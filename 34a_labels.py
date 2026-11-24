#!/usr/bin/env python3
"""One habitat-label table for every comparative layer (ASR, PGLS, RERconverge, selection summaries),
so all analyses use the same deep / shallow-background / unlabelled calls.
  lab         D / S / U exactly as 21_hog_select.labels() (deep = role 'deep'; U = eurybathic, no depth
              data, anchialine, non-marine or parasitic background; S = remaining shallow-marine species)
  lab_strict  as lab, but D only with species- or specimen-level depth evidence (else U)
  lab_1000    as lab, but D only if the median depth is > 1000 m (else U)
  asr         D / S / NA for ancestral-state reconstruction: D = deep; S = shallow-marine background
              plus non-marine species (informative 'not deep' states for the habitat history);
              NA = eurybathic, no depth data, anchialine or parasitic (habitat = host)
  asr_strict  as asr, but deep species without species/specimen depth evidence -> NA
  asr_1000    as asr, but deep species with median depth <= 1000 m -> NA
Writes meta/labels_v2.tsv (every species in species_codes.tsv, with its qc status for filtering).
Usage: python 34a_labels.py"""
import csv, importlib.util, os
from collections import Counter, defaultdict

R = os.path.expanduser("~/deepsea_conv")
spec = importlib.util.spec_from_file_location("hs", f"{R}/scripts/21_hog_select.py")
hs = importlib.util.module_from_spec(spec); spec.loader.exec_module(hs)


def clade(r):
    if r["class"] == "Actinopteri":
        return "Actinopteri"
    if r["class"] in ("Chondrichthyes", "Myxini", "Hyperoartia"):
        return "Chondrichthyes"
    if r["phylum"] in ("Annelida", "Nemertea"):
        return "Annelida"
    return r["phylum"]


def main():
    lab = hs.labels()
    sc = list(csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"))
    qc = {r["code"]: r["status"] for r in csv.DictReader(open(f"{R}/meta/qc_table.tsv"), delimiter="\t")}
    di = defaultdict(set)
    for r in csv.DictReader(open(f"{R}/meta/depth_integrated.tsv"), delimiter="\t"):
        if r["evidence_level"]:
            di[r["species"]].add(r["evidence_level"])
    bg = {r["code"]: r for r in csv.DictReader(open(f"{R}/meta/bg_eligibility.tsv"), delimiter="\t")}
    rank = {"species": 0, "specimen": 1, "locality": 2, "genus": 3}
    cols = ["code", "species", "phylum", "class", "order", "family", "clade", "role", "source", "qc",
            "lab", "lab_strict", "lab_1000", "asr", "asr_strict", "asr_1000",
            "evidence_level", "depth_median", "frac_gt200", "bg_reason"]
    n, miss = Counter(), []
    with open(f"{R}/meta/labels_v2.tsv.tmp", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in sc:
            c = r["code"]
            if c not in lab:          # not labelled by 21_hog_select -> unknown, reported below
                miss.append(c)
            L = lab.get(c, "U")
            ev = min(di.get(r["species"], set()), key=lambda e: rank.get(e, 9), default="")
            try:
                dm = float(r.get("depth_median", ""))
            except ValueError:
                dm = float("nan")
            strict = "U" if L == "D" and ev not in ("species", "specimen") else L
            l1000 = "U" if L == "D" and not dm > 1000 else L
            reason = bg.get(c, {}).get("reason", "")
            # non-marine (WoRMS isMarine != 1) background is 'not deep' for the habitat history; eurybathic,
            # no-depth, anchialine and parasitic (habitat = host) species stay unknown
            notdeep = L == "U" and r["role"] != "deep" and reason.startswith("WoRMS isMarine")
            asr = "D" if L == "D" else ("S" if L == "S" or notdeep else "NA")
            asr_s = "NA" if asr == "D" and strict != "D" else asr
            asr_k = "NA" if asr == "D" and l1000 != "D" else asr
            row = dict(r, clade=clade(r), qc=qc.get(c, "NA"), lab=L, lab_strict=strict, lab_1000=l1000,
                       asr=asr, asr_strict=asr_s, asr_1000=asr_k, evidence_level=ev or "NA",
                       bg_reason=reason or "NA")
            o.write("\t".join(str(row.get(k, "")) or "NA" for k in cols) + "\n")
            if row["qc"] == "pass":
                for k, v in (("lab", L), ("strict", strict), ("1000", l1000), ("asr", asr),
                             ("asr_strict", asr_s), ("asr_1000", asr_k)):
                    n[(k, v)] += 1
    os.replace(f"{R}/meta/labels_v2.tsv.tmp", f"{R}/meta/labels_v2.tsv")
    print("QC-pass label counts:", dict(sorted(n.items())))
    print(f"codes absent from 21_hog_select.labels() (set U): {len(miss)}", " ".join(miss[:20]))


if __name__ == "__main__":
    main()
