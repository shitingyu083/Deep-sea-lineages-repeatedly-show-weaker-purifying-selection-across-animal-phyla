#!/usr/bin/env python3
"""Select near-single-copy HOGs from a finished per-clade OrthoFinder run.
Source: <clade>/of/Results_*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv
Rules (ingroup only; fixed outgroups are dropped, they are too distant for dN/dS):
  - a species with >1 gene in the HOG is kept only if its genes are species-specific
    inparalogs: they form one clade in the OG resolved gene tree restricted to the HOG's
    genes (PhyloPyPruner-style 'collapse'); the longest protein is kept. Otherwise the
    species is dropped from that HOG (paralog-safe). Needed for redundant TSA proteomes.
  - keep the HOG if retained species >= max(4, MINOCC * n_ingroup) and it has
    >= 2 deep and >= 2 shallow species (labels from species_codes.tsv)
Depth labels: deep = role 'deep'; U (unlabelled) = UNLAB set; shallow = all other.
Writes work/codon/<clade>/hogs.tsv (hog, n, nD, nS, species:gene list).
Usage: python 21_hog_select.py CLADE [MINOCC=0.6]"""
import csv, glob, gzip, os, sys

R = os.path.expanduser("~/deepsea_conv")
OUTG = {"Actinopteri": ["Mobhyp", "Lamflu"], "Chondrichthyes": ["Hypnip", "Bralan"],
        "Mollusca": ["Linana", "Galcae"], "Annelida": ["Linana", "Mme"],
        "Arthropoda": ["Pricau", "Anisim"], "Cnidaria": ["Syccil", "Ptyfla", "Halasi"],
        "Porifera": ["Mnelei", "Pocdam", "Clyhem"], "Echinodermata": ["Ptyfla", "Bralan"]}
# eurybathic / depth unknown: kept in alignments, never foreground or background
UNLAB = {"Hha", "Mco", "Xdo"}
# curated habitat overrides (meta/habitat_override.tsv): straddling 200 m or aphotic-shallow
UNLAB_ROLES = {"eury", "anchialine", "nodata"}

def bg_unlab():
    """Non-deep species that are not shallow-marine background (meta/bg_eligibility.tsv:
    WoRMS non-marine, eurybathic, no depth data, anchialine, parasite) -> U."""
    f = f"{R}/meta/bg_eligibility.tsv"
    if not os.path.exists(f):
        return set()
    return {r["code"] for r in csv.DictReader(open(f), delimiter="\t") if r["bg_label"] == "U"}

def labels():
    lab = {}
    BGU = bg_unlab()
    for r in csv.DictReader(open(f"{R}/meta/species_codes.tsv"), delimiter="\t"):
        u = r["code"] in UNLAB or r["role"] in UNLAB_ROLES or (r["role"] != "deep" and r["code"] in BGU)
        lab[r["code"]] = "U" if u else ("D" if r["role"] == "deep" else "S")
    return lab

def parse_nwk(s):
    """Minimal Newick -> (children dict, leaf names). Node ids are ints; labels/lengths ignored."""
    ch, leaves, st, cur, i, n = {}, {}, [], None, 0, 0
    root = n; ch[root] = []; cur = root; n += 1; s = s.strip().rstrip(";")
    if s.startswith("("): i = 1
    while i < len(s):
        c = s[i]
        if c == "(":
            new = n; n += 1; ch[new] = []; ch[cur].append(new); st.append(cur); cur = new; i += 1
        elif c == ",": i += 1
        elif c == ")":
            i += 1
            while i < len(s) and s[i] not in ",()": i += 1   # internal label / length
            cur = st.pop() if st else root
        else:
            j = i
            while j < len(s) and s[j] not in ",():": j += 1
            name = s[i:j]; i = j
            while i < len(s) and s[i] not in ",()": i += 1   # length
            leaf = n; n += 1; ch[leaf] = []; ch[cur].append(leaf); leaves[leaf] = name
    return ch, leaves, root

def clade_sets(ch, leaves, root):
    """leaf-name set below every node (iterative post-order)."""
    out, order, stack = {}, [], [root]
    while stack:
        v = stack.pop(); order.append(v); stack.extend(ch[v])
    for v in reversed(order):
        out[v] = {leaves[v]} if v in leaves else set().union(*(out[c] for c in ch[v]))
    return out

def inparalog_ok(tree_sets, hog_genes, sp_genes):
    """True if sp_genes form a clade (one side of a split) of the tree restricted to hog_genes."""
    if tree_sets is None: return False
    hg, sg = set(hog_genes), set(sp_genes)
    rest = hg - sg
    if not rest: return True
    for ls in tree_sets:
        x = ls & hg
        if x == sg or x == rest: return True
    return False

_LEN = {}
def plen(clade, gene):
    if clade not in _LEN:
        d = {}; base = f"{R}/work/orthofinder/{clade}"
        for fa in glob.glob(f"{base}/input/*.faa") + glob.glob(f"{base}/add_v2/*.faa") \
                + glob.glob(f"{base}/add_v3/*.faa"):
            name = None
            for line in open(fa):
                if line.startswith(">"): name = line[1:].split()[0]; d[name] = 0
                elif name: d[name] += len(line.strip())
        _LEN[clade] = d
    return _LEN[clade].get(gene, 0)

def latest_results(clade):
    # v1 run, -b passes on it (v2, fish sub77), and the panel-v3 -b pass on a v2 -b run (29d_of_add_v3.sh)
    base = f"{R}/work/orthofinder/{clade}"
    res = [r for r in glob.glob(f"{base}/of/Results_*") + glob.glob(f"{base}/add_v2/OrthoFinder/Results_*")
           + glob.glob(f"{base}/of/Results_*/WorkingDirectory/OrthoFinder/Results_*")
           + glob.glob(f"{base}/of/Results_*/WorkingDirectory/OrthoFinder/Results_*/WorkingDirectory/OrthoFinder/Results_*")
           if os.path.exists(f"{r}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")]
    return max(res, key=lambda r: os.path.getmtime(f"{r}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")) if res else None

def main():
    clade = sys.argv[1]; minocc = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6
    res = latest_results(clade)
    if not res:
        sys.exit(f"no OrthoFinder results for {clade}")
    n0 = f"{res}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
    tcache = {}
    def tree_sets(og):
        if og not in tcache:
            if len(tcache) > 5000: tcache.clear()
            fn = f"{res}/Resolved_Gene_Trees/{og}_tree.txt"
            if os.path.exists(fn):
                ch, lv, rt = parse_nwk(open(fn).read())
                tcache[og] = list(clade_sets(ch, lv, rt).values())
            else:
                tcache[og] = None
        return tcache[og]
    lab = labels()
    out = os.environ.get("HOG_OUT", f"{R}/work/codon/{clade}"); os.makedirs(out, exist_ok=True)
    with open(n0) as f:
        rd = csv.reader(f, delimiter="\t"); hdr = next(rd)
        sp_cols = [(i, h) for i, h in enumerate(hdr) if i >= 3 and h not in OUTG[clade]]
        ning = len(sp_cols); need = max(4, minocc * ning)
        kept = 0; ncoll = 0
        with open(f"{out}/hogs.tsv", "w") as fo:
            fo.write("hog\tn\tnD\tnS\tmembers\n")
            for row in rd:
                mem, multi = [], {}
                for i, sp in sp_cols:
                    g = [x for x in (row[i].split(", ") if i < len(row) and row[i] else []) if x]
                    if len(g) == 1:
                        mem.append((sp, g[0]))
                    elif len(g) > 1:
                        multi[sp] = g
                if multi:
                    hog_genes = [x for i in range(3, len(row)) if row[i] for x in row[i].split(", ")]
                    ts = tree_sets(row[1])
                    for sp, g in multi.items():
                        if inparalog_ok(ts, hog_genes, g):
                            mem.append((sp, max(g, key=lambda q: plen(clade, q)))); ncoll += 1
                nD = sum(lab.get(s) == "D" for s, _ in mem); nS = sum(lab.get(s) == "S" for s, _ in mem)
                if len(mem) >= need and nD >= 2 and nS >= 2:
                    kept += 1
                    fo.write(f"{row[0]}\t{len(mem)}\t{nD}\t{nS}\t" + ",".join(f"{s}:{g}" for s, g in mem) + "\n")
    print(f"{clade}\tingroup={ning}\tneed>={need:.1f}\tkept={kept}\tinparalog_collapsed={ncoll}\t{n0}")

if __name__ == "__main__":
    main()
