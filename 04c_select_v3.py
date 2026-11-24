#!/usr/bin/env python3
"""Panel v3 = panel v2 + (i) deep genomes previously excluded by contig N50 < 2 kb (the
compleasm/protein QC gate, not contiguity, now decides; several are chromosome-scale
scaffolds, e.g. Paraliparis) and (ii) backbone taxa needed for fossil calibrations
(Benton et al. 2015 nodes 25/27/34 + Sarcopterygii). Shallow controls keep N50 >= 2 kb.
Writes meta/panel_v3.tsv; panel v2 selection is unchanged for existing families.

Panel v2: every deep (>200 m) metazoan genome + family-matched shallow controls.

Deep set  : all genome species with final_class == deep in meta/depth_integrated.tsv
            (occurrence, sample-depth and congener layers; see 02d). No family cap.
Shallow   : per deep family, up to K shallow species of the same family; if the
            family has < 2 shallow genomes, top up from the same order to 3.
            Ranked by annotation (refseq > genbank > none), assembly level,
            contig N50.
Quality gate is applied later on compleasm (C >= 60 %, M <= 25 %), not on N50;
here we only drop assemblies with contig N50 < 2 kb (unassemblable for miniprot).
Species already coded in panel v1 (meta/species_codes.tsv) keep their code.

Usage: 04b_select_v2.py [K=5]  -> meta/panel_v2.tsv
"""
import csv, sys, os, collections

M = os.path.expanduser('~/deepsea_conv/meta/')
K = int(sys.argv[1]) if len(sys.argv) > 1 else 5
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))

DI = {r['species']: r for r in rd(M + 'depth_integrated.tsv') if r['set'] == 'genome'}
CA = {}
for r in rd(M + 'candidates_all.tsv'):
    CA.setdefault(r['species'], r)
V1 = {r['species']: r for r in rd(M + 'species_codes.tsv')}

ANN = {'refseq': 0, 'genbank': 1}
LEV = {'Chromosome': 0, 'Complete Genome': 0, 'Scaffold': 1, 'Contig': 2}


def n50(r):
    try:
        return float(r['contig_n50'] or 0)
    except ValueError:
        return 0.0


def rank(r):
    return (ANN.get(r['annot'], 2), LEV.get(r['level'], 3), -n50(r))


def ok(sp):
    r = CA.get(sp)
    return r is not None and (n50(r) >= 2000 or r['annot'] in ANN)


deep = [s for s, d in DI.items() if d['final_class'] == 'deep' and s in CA]
# shallow controls must be water-breathing marine taxa (same exclusions as 02d)
SKIP = {'Insecta', 'Arachnida', 'Mammalia', 'Aves', 'Reptilia', 'Amphibia', 'Collembola',
        'Diplopoda', 'Chilopoda', 'Clitellata'}
shal = [s for s, d in DI.items() if d['final_class'] == 'shallow' and ok(s)
        and CA[s]['class'] not in SKIP and CA[s]['phylum']]
byfam, byord = collections.defaultdict(list), collections.defaultdict(list)
for s in shal:
    byfam[CA[s]['family']].append(s)
    byord[CA[s]['order']].append(s)

sel = {s: ('deep', 'all_deep') for s in deep}
byclass = collections.defaultdict(list)
for s in shal:
    byclass[CA[s]['class']].append(s)
for fam in sorted({CA[s]['family'] for s in deep}):
    # an empty family is not a taxon: never match on it, fall back to order/class below
    pool = sorted(byfam.get(fam, []), key=lambda s: rank(CA[s])) if fam else []
    pick = pool[:K]
    if len(pick) < 2:
        d0 = CA[[s for s in deep if CA[s]['family'] == fam][0]]
        o, c = d0['order'], d0['class']
        src = byord.get(o, []) if o else byclass.get(c, []) if c else []
        extra = [s for s in sorted(src, key=lambda s: rank(CA[s])) if s not in pick]
        pick += extra[:3 - len(pick)]
    for s in pick:
        sel.setdefault(s, ('shallow_match', f'family:{fam}' if fam and CA[s]['family'] == fam else f'order_for:{fam or d0["class"]}'))
# v3 backbone additions (role backbone; not part of deep/shallow contrasts)
BACKBONE_V3 = ['Latimeria chalumnae', 'Protopterus annectens', 'Polypterus senegalus',
               'Erpetoichthys calabaricus', 'Lepisosteus oculatus', 'Amia ocellicauda',
               'Acipenser ruthenus', 'Lampris incognitus', 'Regalecus russelii']
for s in BACKBONE_V3:
    if s in CA and s not in sel and s not in V1:
        sel[s] = ('backbone', 'backbone_v3')
# v1 backbone/outgroups are kept as is
for s, r in V1.items():
    if s not in sel:
        sel[s] = (r['role'], 'v1_kept')

cols = ['species', 'code', 'role', 'reason', 'in_v1', 'accession', 'annot', 'level', 'contig_n50',
        'total_len', 'phylum', 'class', 'order', 'family', 'depth_src', 'sample_max_m']
with open(M + 'panel_v3.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(cols)
    for s, (role, why) in sorted(sel.items(), key=lambda x: (CA.get(x[0], {}).get('phylum', ''), x[0])):
        r, d, v = CA.get(s, {}), DI.get(s, {}), V1.get(s)
        if v and not r.get('phylum'):
            r = dict(r, **{k: v[k] for k in ('phylum', 'class', 'order', 'family')})
        w.writerow([s, v['code'] if v else '', role, why, 'yes' if v else 'no',
                    v['accession'] if v else r.get('accession', ''), r.get('annot', ''), r.get('level', ''),
                    r.get('contig_n50', ''), r.get('total_len', ''), r.get('phylum', ''), r.get('class', ''),
                    r.get('order', ''), r.get('family', ''), d.get('final_src', ''), d.get('sample_max_m', '')])
c = collections.Counter((CA.get(s, {}).get('phylum', '?'), role, 'v1' if s in V1 else 'new') for s, (role, _) in sel.items())
for k in sorted(c):
    print(*k, c[k], sep='\t')
print('total', len(sel), 'new', sum(s not in V1 for s in sel))
