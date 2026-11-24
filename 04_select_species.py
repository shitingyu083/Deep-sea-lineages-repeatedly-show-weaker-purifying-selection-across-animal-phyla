#!/usr/bin/env python
"""Select species panel: deep lineages + matched shallow relatives + class backbone.

Rules (fixed before any selection-signal analysis is run):
  quality  : annotated OR contig N50 >= 20 kb
  deep     : best <=3 per family (annotation > level > N50)
  shallow  : per deep family, best <=3 shallow in same family; if none, best <=3
             in same order; class backbone = best 1 annotated shallow per class
Output: meta/panel_v1.tsv
"""
import pandas as pd

df = pd.read_csv("meta/candidates_all.tsv", sep="\t", dtype={"tax_id": str})
m = df[df.in_obis & ~df.air_breather].copy()
m["n50"] = pd.to_numeric(m.contig_n50, errors="coerce").fillna(0)
m["ok"] = (m.n50 >= 20000) | (m.annot != "none")
lv = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}
m["q"] = list(zip(m.annot.map({"refseq": 2, "genbank": 1, "none": 0}),
                  m.level.map(lv).fillna(0), m.n50))
m = m[m.ok & m.family.notna()]

def best(sub, k):
    return sub.sort_values("q", ascending=False).head(k)

deep = m[m.depth_class == "deep"]
sh = m[m.depth_class == "shallow"]
rows = []
for fam, g in deep.groupby("family"):
    rows.append(best(g, 3).assign(role="deep"))
    s = sh[sh.family == fam]
    if len(s) == 0:
        s = sh[sh.order == g.order.iloc[0]]
    if len(s):
        rows.append(best(s, 3).assign(role="shallow_match"))
bb = sh[sh.annot != "none"].groupby("class", group_keys=False).apply(lambda g: best(g, 1))
rows.append(bb.assign(role="backbone"))
p = pd.concat(rows).drop_duplicates("species")
p = p[["role", "phylum", "class", "order", "family", "species", "tax_id", "accession", "annot",
       "level", "contig_n50", "total_len", "n_depth", "depth_median", "frac_gt200"]]
p.sort_values(["phylum", "class", "order", "family", "role"]).to_csv("meta/panel_v1.tsv", sep="\t", index=False)
print(p.role.value_counts().to_dict())
print(p.groupby(["phylum", "role"]).size().unstack(fill_value=0))
print("total Gb:", round(pd.to_numeric(p.total_len).sum() / 1e9, 1))
