#!/usr/bin/env python
"""R2 power simulation v2 (PREREG_nature.md; DEVIATIONS: per-HOG tree-length calibration after the pilot gate). One job = clade hog scenario rep outdir.
Simulates codon data under the RELAX-null fitted Reference omega distribution on the HOG's own tree,
imposes the scenario on FG branches, applies the real gap mask, runs RELAX or BUSTED-PH, writes res.tsv."""
import sys, os, re, json, random, subprocess, hashlib, itertools

B = os.path.expanduser("~/deepsea_conv")
IQ = B + "/envs/dc/bin/iqtree"
HY = B + "/envs/dc/bin/hyphy"
SCEN = {  # name: (test, K, half, epi_p, epi_w)
    "H0r": ("relax", 1.0, 0, 0, 0), "K75": ("relax", 0.75, 0, 0, 0), "K50": ("relax", 0.5, 0, 0, 0),
    "K50h": ("relax", 0.5, 1, 0, 0), "H0b": ("bph", 1.0, 0, 0, 0),
    "P1": ("bph", 1.0, 0, 0.05, 5.0), "P2": ("bph", 1.0, 0, 0.10, 10.0)}

# ---------- newick ----------
class N:
    def __init__(s): s.ch = []; s.name = ""; s.tag = ""; s.bl = 0.0; s.up = None

def parse(t):
    t = t.strip().rstrip(";"); i = 0; root = N(); cur = root; stack = []
    tok = re.compile(r"[^(),:;\[]+")
    while i < len(t):
        c = t[i]
        if c == "(":
            n = N(); n.up = cur; cur.ch.append(n); stack.append(cur); cur = n; i += 1
        elif c == ",":
            n = N(); par = stack[-1]; n.up = par; par.ch.append(n); cur = n; i += 1
        elif c == ")":
            cur = stack.pop(); i += 1
        elif c == ":":
            m = re.match(r":([0-9.eE+-]+)", t[i:]); cur.bl = float(m.group(1)); i += len(m.group(0))
        elif c == "[":
            j = t.index("]", i); i = j + 1
        else:
            m = tok.match(t, i); s = m.group(0); i = m.end()
            mt = re.match(r"(.*?)\{(\w+)\}$", s)
            if mt: cur.name, cur.tag = mt.group(1), mt.group(2)
            else: cur.name = s
    return root

def nodes(n):
    yield n
    for c in n.ch: yield from nodes(c)

def tips(n):
    return [n.name] if not n.ch else [x for c in n.ch for x in tips(c)]

def out(n, lab):
    s = "(" + ",".join(out(c, lab) for c in n.ch) + ")" if n.ch else n.name
    return s + lab(n)

# ---------- codon model ----------
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSSRRRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODE = dict(zip(CODONS, AA))
SENSE = [c for c in CODONS if CODE[c] != "*"]
TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

def f3x4(seqs):
    cnt = [[0] * 4 for _ in range(3)]
    for s in seqs.values():
        for k in range(0, len(s) - 2, 3):
            cd = s[k:k + 3]
            if cd in CODE and CODE[cd] != "*":
                for p in range(3): cnt[p][BASES.index(cd[p])] += 1
    return [[(x + 1) / (sum(r) + 4) for x in r] for r in cnt]

def rate_parts(f, kappa):
    pi = {c: f[0][BASES.index(c[0])] * f[1][BASES.index(c[1])] * f[2][BASES.index(c[2])] for c in SENSE}
    z = sum(pi.values()); pi = {c: v / z for c, v in pi.items()}
    S = Nn = 0.0
    for i in SENSE:
        for j in SENSE:
            d = [p for p in range(3) if i[p] != j[p]]
            if len(d) != 1: continue
            p = d[0]; q = f[p][BASES.index(j[p])] * (kappa if (i[p], j[p]) in TS else 1.0)
            if CODE[i] == CODE[j]: S += pi[i] * q
            else: Nn += pi[i] * q
    return S, Nn

def read_fas(p):
    d = {}; k = None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = []
        elif k: d[k].append(l)
    return {k: "".join(v).upper() for k, v in d.items()}

def seeded(*a):
    return int(hashlib.md5("|".join(map(str, a)).encode()).hexdigest()[:8], 16)


def iqfit(aln_path, topo, pre):
    subprocess.run([IQ, "-s", aln_path, "-te", topo, "-m", "GY+F3X4", "-st", "CODON", "-T", "1", "--prefix", pre,
                    "-redo", "-quiet", "-seed", "20261028"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    txt = open(pre + ".iqtree").read()
    return (float(re.search(r"\(omega\):\s*([0-9.]+)", txt).group(1)), float(re.search(r"\(kappa\):\s*([0-9.]+)", txt).group(1)),
            float(re.search(r"Total tree length \(sum of branch lengths\):\s*([0-9.]+)", txt).group(1)))

def setup(C, hog, fitd):
    W = f"{B}/work/codon/{C}"  # trees + RELAX fits of the analysed labelling (Actinopteri_m has its own)
    WA = f"{B}/work/codon/{'Actinopteri' if C.startswith('Actinopteri') else C}"  # alignments
    tr = parse(open(f"{W}/trees/{hog}.hyphy.nwk").read()); aln = read_fas(f"{WA}/aln/{hog}.fas")
    names = tips(tr); L = min(len(aln[n]) for n in names) // 3
    os.makedirs(fitd, exist_ok=True); pre = f"{fitd}/fit"
    if not os.path.exists(pre + ".iqtree"):
        open(f"{fitd}/topo.nwk", "w").write(out(tr, lambda n: "") + ";")
        iqfit(f"{WA}/aln/{hog}.fas", f"{fitd}/topo.nwk", pre)
    txt = open(pre + ".iqtree").read()
    w_hat = float(re.search(r"\(omega\):\s*([0-9.]+)", txt).group(1)); kap = float(re.search(r"\(kappa\):\s*([0-9.]+)", txt).group(1))
    fitL = float(re.search(r"Total tree length \(sum of branch lengths\):\s*([0-9.]+)", txt).group(1))
    ft = parse(open(pre + ".treefile").read()); allt = frozenset(names); split = {}
    for n in nodes(ft):
        if n.up is None: continue
        s = frozenset(tips(n)); split[s] = n.bl; split[allt - s] = n.bl
    for n in nodes(tr):
        if n.up is None: continue
        n.bl = split.get(frozenset(tips(n)), 0.0)
        if n.up is tr and len(tr.ch) == 2: n.bl /= 2.0
    j = json.load(open(next(p for p in (f"{W}/sel/{hog}.relax.json", f"{W}/sel_rescue/{hog}.relax.json") if os.path.exists(p))))
    rd = j["fits"]["RELAX null"]["Rate Distributions"]["Reference"]
    cls = [(float(v["omega"]), float(v["proportion"])) for _, v in sorted(rd.items())]
    f = f3x4({n: aln[n] for n in names}); S, N1 = rate_parts(f, kap)
    return dict(W=W, tr=tr, aln=aln, names=names, L=L, w_hat=w_hat, kap=kap, fitL=fitL, cls=cls, S=S, N1=N1,
                fstr=",".join(f"{x:.6f}" for r in f for x in r), fitd=fitd)

def simulate(c, sc_cls, hit, ep, seed, scale, path, wd):
    os.makedirs(wd, exist_ok=True); tr, names, L = c["tr"], c["names"], c["L"]
    S, N1, w_hat, kap, fstr = c["S"], c["N1"], c["w_hat"], c["kap"], c["fstr"]
    nc = [int(round(p * L)) for _, _, p in sc_cls]; nc[0] += L - sum(nc); cols = []
    for ci, ((wb, wf, p), n_c) in enumerate(zip(sc_cls, nc)):
        if n_c <= 0: continue
        def lab(n, wb=wb, wf=wf):
            if n.up is None: return ""
            w = wf if (n.tag == "FG" if ep else id(n) in hit) else wb
            return f":{n.bl * scale * (S + w * N1) / (S + w_hat * N1):.8f}[&model=GY{{{max(w,1e-4):.6f},{kap:.4f}}}+F3X4{{{fstr}}}]"
        pre = f"{wd}/c{ci}"; open(pre + ".nwk", "w").write(out(tr, lab) + ";")
        subprocess.run([IQ, "--alisim", pre, "-t", pre + ".nwk", "-m", f"GY{{{max(wb,1e-4):.6f},{kap:.4f}}}+F3X4{{{fstr}}}",
                        "--length", str(n_c * 3), "-seed", str(seed + ci), "-af", "fasta", "-redo"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        sq = read_fas(pre + ".fa")
        for k_ in range(n_c): cols.append({t: sq[t][3 * k_:3 * k_ + 3] for t in names})
        for x in (".fa", ".nwk", ".nwk.log"):
            if os.path.exists(pre + x): os.remove(pre + x)
    random.Random(seed).shuffle(cols)
    with open(path, "w") as fo:
        for t in names:
            s = c["aln"][t]
            fo.write(f">{t}\n" + "".join("---" if any(ch in "-N?" for ch in s[3*k:3*k+3]) or len(s[3*k:3*k+3]) < 3 else col[t]
                                         for k, col in enumerate(cols)) + "\n")

def calibrate(c):
    """Scale the fitted lengths until an H0 simulation refitted by the same single-omega GY model reproduces the fitted
    tree length (within 5 %, max 4 rounds): the 3-class mixture saturates the high-omega class, so unscaled simulations
    come out 17-26 % short (pilot, DEVIATIONS_addon.md)."""
    fd = c["fitd"]; sf = f"{fd}/scale.txt"
    if os.path.exists(sf): return float(open(sf).read().split()[0])
    scale = 1.0; logl = []
    sc_cls = [(w, w, p) for w, p in c["cls"]]
    for it in range(4):
        simulate(c, sc_cls, set(), 0, 20261028 + it, scale, f"{fd}/cal.fas", f"{fd}/calwd")
        _, _, simL = iqfit(f"{fd}/cal.fas", f"{fd}/topo.nwk", f"{fd}/cal{it}")
        logl.append(f"{it}\t{scale:.4f}\t{simL:.4f}\t{c['fitL']:.4f}")
        if abs(simL / c["fitL"] - 1) < 0.05: break
        scale *= c["fitL"] / simL
    open(f"{fd}/cal_log.tsv", "w").write("iter\tscale\tsimL\tfitL\n" + "\n".join(logl) + "\n")
    open(sf, "w").write(f"{scale:.6f}\t{logl[-1].split()[2]}\t{c['fitL']}\n")
    return scale

def main():
    if sys.argv[1] == "calib":                # calib CLADE HOG FITDIR
        c = setup(sys.argv[2], sys.argv[3], sys.argv[4]); calibrate(c); return
    C, hog, sc, rep, od = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
    test, K, half, ep, ew = SCEN[sc]
    os.makedirs(od, exist_ok=True)
    if os.path.exists(f"{od}/res.tsv"): return
    c = setup(C, hog, f"{od}/../_fit"); scale = calibrate(c)
    seed = seeded(C, hog, sc, rep); rng = random.Random(seed); tr = c["tr"]; cls = c["cls"]
    fg = [n for n in nodes(tr) if n.up is not None and n.tag == "FG"]; hit = set(id(n) for n in fg)
    if half:
        k = sorted(fg, key=lambda n: ",".join(sorted(tips(n)))); rng.shuffle(k); hit = set(id(n) for n in k[: max(1, len(k) // 2)])
    sc_cls = ([(w, w, p * (1 - ep)) for w, p in cls] + [(cls[1][0], ew, ep)]) if ep > 0 else [(w, w ** K, p) for w, p in cls]
    simulate(c, sc_cls, hit, ep, seed, scale, f"{od}/sim.fas", f"{od}/wd")
    tf = f"{c['W']}/trees/{hog}.hyphy.nwk"; os.chdir(od)
    if test == "relax":
        cmd = [HY, "CPU=1", "relax", "--alignment", "sim.fas", "--tree", tf, "--test", "FG", "--reference", "BG", "--models", "Minimal", "--output", "o.json"]
    else:
        cmd = [HY, "CPU=1", "busted-ph", "--alignment", "sim.fas", "--tree", tf, "--branches", "FG", "--srv", "No", "--output", "o.json"]
    subprocess.run(cmd, stdout=open("hy.log", "w"), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    r = json.load(open("o.json")); t_ = r.get("test results", {})
    vals = ([t_.get("p-value"), t_.get("relaxation or intensification parameter"), "", ""] if test == "relax" else
            [t_.get("p-value"), "", r.get("Background selection test results", {}).get("p-value"),
             r.get("Comparative selection test results", {}).get("p-value")])
    with open("res.tsv", "w") as fo:
        fo.write("clade\thog\tscen\trep\tp\tK\tbg_p\tdiff_p\tL\tw_hat\tkappa\tscale\tseed\tn_fg\tn_hit\n")
        fo.write("\t".join(map(str, [C, hog, sc, rep] + vals + [c["L"], c["w_hat"], c["kap"], scale, seed, len(fg), len(hit)])) + "\n")
    if os.path.exists("sim.fas"): os.system("gzip -f sim.fas")

if __name__ == "__main__":
    main()
