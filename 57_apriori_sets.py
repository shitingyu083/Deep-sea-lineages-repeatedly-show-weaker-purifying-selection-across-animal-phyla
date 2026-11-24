#!/usr/bin/env python3
"""A5 (PREREG_addon.md): a-priori pathway sets with positive and negative controls, on 47's HOG annotation.
Sets: a HOG is in a set when its 47 annotation (work/func/<C>/hog_annot.tsv; >= 50 % of the annotated members) carries
any of the set's KEGG maps (kegg_pathway) or KOs (kegg_ko, 'ko:' prefix dropped). P1, P2 positive controls, N1
negative control, S1-S9 candidates (SETS below).
Indicators per HOG and label set (work/codon/<C>/sel_results_v3.tsv; the paired frames of 47 / 25b): relax_relaxed
(RELAX p < 0.05 and K < 1), relax_intens (p < 0.05 and K > 1), bph_strict (BUSTED-PH FG p < 0.05 and diff p < 0.05)
on the nullhy HOGs with a real and a null1 value; cs_qual (>= 1 qualified CSUBST pair) on the hylist HOGs with a real
and >= 1 null value.
Per clade, set and indicator (>= 3 set HOGs in the frame): S = sum of d = I_real - mean I_null; null = 47's
exchangeable permutation (per HOG one of its 1 + K sets drawn as real; perm_term copied from 47_func_enrich.py),
10,000 draws, seed (20261025, crc32(clade), indicator index); z = (S - mean S*) / sd S*, one-sided p (S* >= S);
clades with sd S* = 0 carry no information and are dropped. Across clades as
48: weighted Stouffer Z = sum(w z) / sqrt(sum w^2), w = sqrt(n HOGs). p_up / p_two are exact (Monte Carlo) p of Z
under the joint permutation null (draw i of each clade's independent stream = joint draw i, the same z
standardisation), because sparse sets give discrete z for which the normal tail is anti-conservative; the normal
p are reported alongside.
Families (Holm over the testable members): positive P1, P2 x relax_relaxed (p_up, over 2); negative N1 x every
indicator (p_two, over 4; the control passes when every Holm p >= 0.05, a 5 % false-fail rate under the null);
candidates S1-S9 x {bph_strict, cs_qual} (p_up, over 18); the rest descriptive.
Usage: python 57_apriori_sets.py [--clades A,B,...] [--func DIR] [--codon DIR] [--out DIR] [--nperm 10000]
       [--sub norescue]   (reads CODON/<C>/norescue/sel_results_v3.tsv, writes OUT/norescue)
Writes OUT/set_members.tsv, OUT/perclade.tsv, OUT/combined.tsv, OUT/report.txt, OUT/A5.DONE."""
import csv, math, os, sys, zlib
from collections import defaultdict
import numpy as np

nsf = lambda x: 0.5 * math.erfc(x / math.sqrt(2))   # standard normal upper tail

R = os.environ.get("DSC_ROOT", os.path.expanduser("~/deepsea_conv"))
A = sys.argv
opt = lambda k, d: A[A.index(k) + 1] if k in A else d
FUNC = opt("--func", f"{R}/work/func")
COD = opt("--codon", f"{R}/work/codon")
OUT = opt("--out", f"{R}/work/addon/A5")
CL = opt("--clades", "Actinopteri_m,Chondrichthyes,Mollusca,Annelida,Arthropoda,Cnidaria,Porifera,Echinodermata").split(",")
NPERM = int(opt("--nperm", "10000"))
SEED, CHUNK, MINH = 20261025, 1000, 3
NA = float("nan")
TESTS = ("relax_relaxed", "relax_intens", "bph_strict", "cs_qual")
SETS = {  # name: (label, KEGG maps, KOs)
    "P1": ("phototransduction", {"map04744", "map04745"}, set()),
    "P2": ("circadian rhythm / entrainment", {"map04710", "map04711", "map04713"}, set()),
    "N1": ("ribosome", {"map03010"}, set()),
    "S1": ("DNA repair", {"map03410", "map03420", "map03430", "map03440", "map03450"}, set()),
    "S2": ("unsaturated fatty-acid biosynthesis", {"map01040"}, set()),
    "S3": ("fatty-acid elongation", {"map00062"}, set()),
    "S4": ("proteasome", {"map03050"}, set()),
    "S5": ("ER protein processing + Hsp70 / Hsp90", {"map04141"}, {"K03283", "K04079"}),
    "S6": ("oxidative phosphorylation", {"map00190"}, set()),
    "S7": ("actin cytoskeleton + actin / tubulin", {"map04810"}, {"K05692", "K07374", "K07375"}),
    "S8": ("Na+/K+-ATPase", set(), {"K01539", "K01540"}),
    "S9": ("TMAO synthesis (FMO)", set(), {"K00485"}),
}
FAM = {("P1", "relax_relaxed"): "positive", ("P2", "relax_relaxed"): "positive"}
FAM.update({("N1", t): "negative" for t in TESTS})
FAM.update({(f"S{i}", t): "candidate" for i in range(1, 10) for t in ("bph_strict", "cs_qual")})
csv.field_size_limit(1 << 30)


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NA


def sel_indicators(fn):
    """= 47 sel_indicators with RELAX split by K: test -> hog -> [real, null...]"""
    by = defaultdict(dict)
    with open(fn) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by[r["hog"]][r["set"]] = r
    out = {t: {} for t in TESTS}
    for hog, s in by.items():
        real = s.get("real")
        if real is None:
            continue
        n1 = s.get("null1")
        if n1 is not None and real.get("in_nullhy") == "1":
            def hy(r):
                rp, rk, bp, dp = fnum(r["relax_p"]), fnum(r["relax_K"]), fnum(r["bph_p"]), fnum(r["bph_diff_p"])
                ok = rp == rp and rk == rk
                return (float(rp < 0.05 and rk < 1) if ok else NA, float(rp < 0.05 and rk > 1) if ok else NA,
                        NA if (bp != bp or dp != dp) else float(bp < 0.05 and dp < 0.05))
            a, b = hy(real), hy(n1)
            for i, t in enumerate(TESTS[:3]):
                if a[i] == a[i] and b[i] == b[i]:
                    out[t][hog] = [a[i], b[i]]
        if real.get("in_hylist") == "1":
            v = [fnum(r["cs_nqual"]) for k, r in sorted(s.items(), key=lambda kv: (kv[0] != "real", kv[0]))
                 if k == "real" or k.startswith("null")]
            ok = [x for x in v if x == x]
            if v[0] == v[0] and len(ok) >= 2:
                out["cs_qual"][hog] = [float(x >= 1) for x in ok]
    return out


def perm_term(groups, terms, rng, keep=None):
    """= 47_func_enrich.perm_term. groups: hog -> indicator list; terms: term -> list of hogs -> (S, mean, sd, p).
    57 addition: if keep is a dict, keep[term] = the NPERM permutation values S* (same draws, same summaries)."""
    hogs = sorted(groups)
    hi = {h: i for i, h in enumerate(hogs)}
    D = []
    for h in hogs:
        g = np.asarray(groups[h], float); n = len(g)
        D.append(g - (g.sum() - g) / (n - 1))
    obs_d = np.array([d[0] for d in D])
    tl = sorted(terms)
    M = np.zeros((len(tl), len(hogs)), np.float32)
    for k, t in enumerate(tl):
        M[k, [hi[h] for h in terms[t]]] = 1
    S = M @ obs_d.astype(np.float32)
    var = np.array([d.var() for d in D])
    live = np.where(var > 0)[0]
    Ml = M[:, live]
    const = M[:, np.setdiff1d(np.arange(len(hogs)), live)] @ obs_d[np.setdiff1d(np.arange(len(hogs)), live)].astype(np.float32)
    ge = np.zeros(len(tl)); s1 = np.zeros(len(tl)); s2 = np.zeros(len(tl))
    Dl = [D[i] for i in live]
    sizes = np.array([len(d) for d in Dl])
    Dpad = np.zeros((len(Dl), sizes.max() if len(Dl) else 1), np.float32)
    for i, d in enumerate(Dl):
        Dpad[i, :len(d)] = d
    done = 0; chunks = []
    while done < NPERM:
        b = min(CHUNK, NPERM - done)
        pick = (rng.random((b, len(Dl))) * sizes).astype(int)
        X = Dpad[np.arange(len(Dl))[None, :], pick]
        Ss = X @ Ml.T + const[None, :]
        ge += (Ss >= S[None, :] - 1e-6).sum(0); s1 += Ss.sum(0); s2 += (Ss.astype(np.float64) ** 2).sum(0)
        if keep is not None:
            chunks.append(Ss.astype(np.float64))
        done += b
    mu = s1 / NPERM; sd = np.sqrt(np.maximum(s2 / NPERM - mu ** 2, 0))
    p = (1 + ge) / (1 + NPERM)
    if keep is not None:
        allS = np.vstack(chunks)
        for k, t in enumerate(tl):
            keep[t] = allS[:, k]
    return {t: (float(S[k]), float(mu[k]), float(sd[k]), float(p[k])) for k, t in enumerate(tl)}


def holm(p):
    p = np.asarray(p, float); m = len(p); q = np.full(m, NA)
    ok = np.where(p == p)[0]; k = len(ok)
    if not k:
        return q
    o = ok[np.argsort(p[ok])]; r = np.minimum(1, p[o] * (k - np.arange(k)))
    q[o] = np.maximum.accumulate(r)
    return q


def set_members(fn):
    """hog_annot.tsv -> set -> {hog: (og, name, matched maps/KOs)}"""
    out = {s: {} for s in SETS}
    with open(fn) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            kp = set(x for x in r["kegg_pathway"].split(",") if x)
            ko = set(x.split(":", 1)[-1] for x in r["kegg_ko"].split(",") if x)
            for s, (_, maps, kos) in SETS.items():
                hit = sorted(kp & maps) + sorted(ko & kos)
                if hit:
                    out[s][r["hog"]] = (r["og"], r["name"], ",".join(hit))
    return out


def main():
    sub = opt("--sub", "")
    out = f"{OUT}/{sub}" if sub else OUT
    os.makedirs(out, exist_ok=True)
    pc, rep, mem, ZD = [], [], [], {}
    for c in CL:
        fa = f"{FUNC}/{c}/hog_annot.tsv"
        fs = f"{COD}/{c}/{sub + '/' if sub else ''}sel_results_v3.tsv"
        if not (os.path.exists(fa) and os.path.exists(fs)):
            rep.append(f"{c}: missing {fa if not os.path.exists(fa) else fs}, skipped"); continue
        ar = f"{FUNC}/{c}/annot_report.txt"
        rep.append(f"{c}: " + (open(ar).readline().strip() if os.path.exists(ar) else "annot_report.txt missing"))
        sm = set_members(fa)
        ind = sel_indicators(fs)
        for s in SETS:
            for h, (og, nm, hit) in sorted(sm[s].items()):
                mem.append([c, s, h, og, nm, hit] + [int(h in ind[t]) for t in TESTS])
        for ti, t in enumerate(TESTS):
            g = ind[t]
            terms = {s: sorted(h for h in sm[s] if h in g) for s in SETS}
            terms = {s: v for s, v in terms.items() if len(v) >= MINH}
            for s in SETS:
                if s not in terms:
                    pc.append([c, s, t, sum(h in g for h in sm[s]), NA, NA, NA, NA, NA, NA, NA, "lt3"])
            if not terms:
                continue
            u = sorted(set(h for v in terms.values() for h in v))
            rng = np.random.default_rng([SEED, zlib.crc32(c.encode()), ti])
            kd = {}
            res = perm_term({h: g[h] for h in u}, terms, rng, keep=kd)
            for s, v in terms.items():
                S, mu, sd, p = res[s]
                nreal = int(sum(g[h][0] for h in v)); nnull = float(sum(np.mean(g[h][1:]) for h in v))
                z = (S - mu) / sd if sd > 0 else NA
                pc.append([c, s, t, len(v), nreal, nnull, S, mu, sd, z, p, "ok" if sd > 0 else "sd0"])
                if sd > 0:
                    ZD[(c, s, t)] = (kd[s] - mu) / sd   # permutation z values of this clade
        print(f"{c}: " + ", ".join(f"{t} {len(ind[t])}" for t in TESTS) + " HOGs in frame; set HOGs "
              + " ".join(f"{s}:{len(sm[s])}" for s in SETS), flush=True)
    comb = []
    for s in SETS:
        for t in TESTS:
            rows = [r for r in pc if r[1] == s and r[2] == t and r[11] == "ok"]
            fam = FAM.get((s, t), "descriptive")
            if rows:
                w = np.sqrt([r[3] for r in rows]); z = np.array([r[9] for r in rows]); den = np.sqrt((w ** 2).sum())
                Z = float((w * z).sum() / den)
                pu, p2 = nsf(Z), 2 * nsf(abs(Z))
                # exact null of Z: the clade streams are independent, so draw i of every clade is one joint draw
                Zs = sum(wc * ZD[(r[0], s, t)] for wc, r in zip(w, rows)) / den
                ppu = (1 + int((Zs >= Z - 1e-9).sum())) / (1 + NPERM)
                ppl = (1 + int((Zs <= Z + 1e-9).sum())) / (1 + NPERM)
                pp2 = min(1.0, 2 * min(ppu, ppl))
            else:
                Z = pu = p2 = ppu = pp2 = NA
            comb.append(dict(set=s, label=SETS[s][0], test=t, family=fam, n_clades=len(rows),
                             n_clades_p05=sum(r[10] < 0.05 for r in rows), n_hogs=sum(r[3] for r in rows),
                             n_real=sum(r[4] for r in rows), n_null=sum(r[5] for r in rows), Z=Z, p_up=ppu,
                             p_two=pp2, p_up_normal=pu, p_two_normal=p2, p_adj=NA, call=""))
    for fam, pk in (("positive", "p_up"), ("candidate", "p_up"), ("negative", "p_two")):
        ix = [i for i, r in enumerate(comb) if r["family"] == fam]
        q = holm([comb[i][pk] for i in ix])
        for i, qq in zip(ix, q):
            r = comb[i]; r["p_adj"] = float(qq)
            if fam == "negative":   # Holm over the 4 indicators, two-sided: 5 % false-fail rate under the null
                r["call"] = "not testable" if qq != qq else ("FAIL (departs from null)" if qq < 0.05 else "pass")
            else:
                r["call"] = "not testable" if qq != qq else ("enriched" if qq < 0.05 else "n.s.")
    neg = [r for r in comb if r["family"] == "negative"]
    with open(f"{out}/set_members.tsv", "w") as f:
        f.write("clade\tset\thog\tog\tname\tmatched\t" + "\t".join(f"in_{t}" for t in TESTS) + "\n")
        for m in mem:
            f.write("\t".join(map(str, m)) + "\n")
    with open(f"{out}/perclade.tsv", "w") as f:
        f.write("clade\tset\ttest\tn_hogs\tn_real\tn_null\tS\tS_null_mean\tS_null_sd\tz\tperm_p_up\tstatus\n")
        for r in pc:
            f.write("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    keys = list(comb[0])
    with open(f"{out}/combined.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in comb:
            f.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    with open(f"{out}/report.txt", "w") as f:
        f.write(f"A5 a-priori sets | sel {sub or 'main (rescued)'} | NPERM {NPERM} | seed {SEED} | min {MINH} HOGs per clade\n")
        f.write("annotation (47 annot_report, first line):\n" + "\n".join("  " + x for x in rep) + "\n")
        f.write("p = joint-permutation p of the weighted Stouffer Z (normal approximation in brackets)\n")
        for fam in ("positive", "negative", "candidate"):
            f.write(f"{fam} family ({'two-sided' if fam == 'negative' else 'one-sided'} p, Holm within family):\n")
            for r in comb:
                if r["family"] == fam:
                    f.write(f"  {r['set']} {r['label']} x {r['test']}: clades {r['n_clades']} (p<0.05 {r['n_clades_p05']}), "
                            f"HOGs {r['n_hogs']}, real {r['n_real']} vs null {r['n_null']:.1f}, Z {r['Z']:.3g}, "
                            f"p_up {r['p_up']:.3g} [{r['p_up_normal']:.3g}], p_two {r['p_two']:.3g} [{r['p_two_normal']:.3g}], "
                            f"Holm {r['p_adj']:.3g} -> {r['call']}\n")
        nt = [r for r in neg if r["call"] != "not testable"]
        f.write(f"negative control N1: {'passes' if nt and all(r['call'] == 'pass' for r in nt) else 'FAILS' if nt else 'not testable'}"
                f" ({sum(r['call'] == 'pass' for r in nt)}/{len(nt)} indicators with Holm p >= 0.05; unadjusted p_two < 0.05: "
                f"{sum(r['p_two'] < 0.05 for r in nt)})\n")
    open(f"{out}/A5.DONE", "w").write("done\n")
    print(open(f"{out}/report.txt").read(), flush=True)


if __name__ == "__main__":
    main()
