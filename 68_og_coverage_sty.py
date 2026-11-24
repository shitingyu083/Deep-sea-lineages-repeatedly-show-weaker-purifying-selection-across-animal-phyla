#!/usr/bin/env python3
"""Post hoc (2026-10-29): number of clades in which each eggNOG OG was testable
(relax test, HyPhy frame) -> distribution of coverage for the recurrence frame.
Read-only on work/func_v2; writes work/addon/R4_direction/og_coverage_*.tsv"""
import csv, collections, os
CL = ['Chondrichthyes', 'Porifera', 'Arthropoda', 'Echinodermata', 'Cnidaria',
      'Annelida', 'Mollusca', 'Actinopteri_m']
B = 'work/func_v2'
OUT = 'work/addon/R4_direction'
os.makedirs(OUT, exist_ok=True)
summ = []
for test in ['relax', 'bph_strict', 'cs_qual']:
    cov = collections.defaultdict(set)
    for c in CL:
        og = {}
        with open(f'{B}/{c}/hog_annot.tsv') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                og[r['hog']] = r['og']
        with open(f'{B}/{c}/hog_sel.tsv') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if r['test'] != test:
                    continue
                o = og.get(r['hog'], '')
                if o.startswith('Metazoa'):   # same key rule as 48_func_combine.py
                    cov[o].add(c)
    dist = collections.Counter(len(v) for v in cov.values())
    for k in range(1, 9):
        summ.append((test, k, dist.get(k, 0), sum(dist.get(j, 0) for j in range(k, 9))))
    print(test, 'n_og', len(cov), dict(sorted(dist.items())))
with open(f'{OUT}/og_coverage_summary.tsv', 'w') as f:
    f.write('test\tn_clades_tested\tn_og\tn_og_ge\n')
    for r in summ:
        f.write('\t'.join(map(str, r)) + '\n')
