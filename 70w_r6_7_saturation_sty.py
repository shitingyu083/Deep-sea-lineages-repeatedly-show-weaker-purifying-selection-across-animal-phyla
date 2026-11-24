# R6-7 extraction (ADDENDUM_R6_7.md md5 3250794e, POST HOC saturation audit). Read-only. Run from ~/deepsea_conv/work.
import json, os, statistics as st
CL = ["Actinopteri_m", "Annelida", "Arthropoda", "Chondrichthyes", "Cnidaria", "Echinodermata", "Mollusca", "Porifera"]
os.makedirs("addon/R6/R6_7", exist_ok=True)
out = open("addon/R6/R6_7/r6_7_fits.tsv", "w")
out.write("clade\thog\tlabel\tp\tK\tn_br\tD_gene\tD_fg\tD_ref\tmax_br\tn_fg\n")
miss = {}
for c in CL:
    for lab, lst, suf in (("real", "hylist.txt", ""), ("null1", "nullhy.txt", ".null1")):
        hogs = [l.split()[0] for l in open(f"codon/{c}/tiers/{lst}") if l.strip()]
        for h in hogs:
            f = f"codon/{c}/sel/{h}{suf}.relax.json"
            try:
                j = json.load(open(f))
                tr = j["test results"]; p = tr["p-value"]; K = tr["relaxation or intensification parameter"]
                ba = j["branch attributes"]["0"]; te = j["tested"]["0"]
            except Exception as e:
                miss[(c, lab)] = miss.get((c, lab), 0) + 1; continue
            L = {b: v["Nucleotide GTR"] for b, v in ba.items() if isinstance(v, dict) and "Nucleotide GTR" in v}
            fg = [L[b] for b in L if te.get(b) == "Test"]; rf = [L[b] for b in L if te.get(b) == "Reference"]
            if not L or not fg:
                miss[(c, lab)] = miss.get((c, lab), 0) + 1; continue
            out.write("\t".join(map(str, [c, h, lab, p, K, len(L), st.mean(L.values()), st.mean(fg),
                                          st.mean(rf) if rf else "NA", max(L.values()), len(fg)])) + "\n")
out.close()
print("missing/unusable JSON by clade,label:", miss)
