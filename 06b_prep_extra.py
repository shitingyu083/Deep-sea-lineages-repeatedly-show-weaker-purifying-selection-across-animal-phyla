#!/usr/bin/env python3
"""Standardise author-released (figshare/Zenodo) proteomes+CDS for panel v2.

Same pairing/verification as 06_prep_local.py (pep<->cds by order, else by ID;
CDS translation must agree with the peptide), plus a generic isoform collapse
(longest isoform per gene for IDs ending in .tN / -RA / -TN / _iN / .mRNA N).
Input table meta/extra_local.tsv: code species role pep cds  (paths relative to
data/extra/). Outputs data/proteomes/<CODE>.faa, data/cds/<CODE>.fna,
meta/idmap/<CODE>.tsv and appends meta/extra_prep_report.tsv.
Usage: 06b_prep_extra.py [CODE ...]   (default: all rows)
"""
import importlib.util, os, re, sys, csv

R = os.path.expanduser('~/deepsea_conv')
spec = importlib.util.spec_from_file_location('pl', f'{R}/scripts/06_prep_local.py')
pl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pl)

ISO = re.compile(r'(\.t\d+|-R[A-Z]{1,2}|-T\d+|_i\d+|\.mRNA\d*|\.p\d+)$')


def gene_key(pid):
    return ISO.sub('', pid)


def rescue(pep, tx):
    q = pep.replace('*', '').upper()
    if len(q) < 30:
        return None
    for f in (0, 1, 2):
        t = pl.translate(tx[f:])
        i = t.find(q)
        if i >= 0:
            a = f + 3 * i
            b = a + 3 * len(q)
            if pl.translate(tx[b:b + 3]) == '*':
                b += 3
            return tx[a:b]
    return None


def prep(code, sp, role, pf, cf, rep):
    pep, cds = pl.read_fa(f'{R}/data/extra/{pf}'), pl.read_fa(f'{R}/data/extra/{cf}')
    pep = [(h, q.replace('.', '*')) for h, q in pep]  # some releases use '.' as stop
    pairs, mode = [], 'order'
    if len(pep) == len(cds):
        pairs = list(zip(pep, cds))
        if sum(pl.agree(p[1], c[1]) for p, c in pairs[:500]) < 0.9 * min(500, len(pairs)):
            pairs = []
    if not pairs:
        mode = 'id'
        cmap = {pl.norm_id(h, code): (h, s) for h, s in cds}
        for h, s in pep:
            k = pl.norm_id(h, code)
            if k not in cmap:
                k = re.sub(r'(\.\d+)+$', '', k)
            if k in cmap:
                pairs.append(((h, s), cmap[k]))
    good = []
    for p, c in pairs:
        if pl.agree(p[1], c[1]):
            good.append((p, c))
        else:  # MAKER-style 'transcript' FASTA carries UTRs: locate the CDS in-frame
            sub = rescue(p[1], c[1])
            if sub:
                good.append((p, (c[0], sub)))
    best = {}
    for (ph, ps), (ch, cs) in good:
        pid = pl.norm_id(ph, code)
        gm = re.search(r'\bgene=(\S+)', ph)  # explicit gene attribute wins
        g, L = (gm.group(1) if gm else gene_key(pid)), len(ps.replace('*', ''))
        if L >= 30 and (g not in best or L > best[g][0]):
            best[g] = (L, pid, ps, cs)
    with open(f'{R}/data/proteomes/{code}.faa', 'w') as fa, open(f'{R}/data/cds/{code}.fna', 'w') as fn, \
            open(f'{R}/meta/idmap/{code}.tsv', 'w') as mp:
        for n, g in enumerate(sorted(best), 1):
            L, pid, ps, cs = best[g]
            nid, c = f'{code}_{n:06d}', cs.upper()
            if pl.translate(c[-3:]) == '*':
                c = c[:-3]
            fa.write(f'>{nid}\n{ps.replace("*", "")}\n')
            fn.write(f'>{nid}\n{c}\n')
            mp.write(f'{nid}\t{pid}\t{g}\n')
    row = [code, sp, role, len(pep), len(cds), mode, len(good), len(best)]
    rep.write('\t'.join(map(str, row)) + '\n')
    print(*row, flush=True)


def main():
    want = set(sys.argv[1:])
    rows = list(csv.DictReader(open(f'{R}/meta/extra_local.tsv'), delimiter='\t'))
    new = not os.path.exists(f'{R}/meta/extra_prep_report.tsv')
    with open(f'{R}/meta/extra_prep_report.tsv', 'a') as rep:
        if new:
            rep.write('code\tspecies\trole\tn_pep\tn_cds\tpairing\tn_agree\tn_genes_kept\n')
        for r in rows:
            if not want or r['code'] in want:
                prep(r['code'], r['species'], r['role'], r['pep'], r['cds'], rep)


if __name__ == '__main__':
    main()
