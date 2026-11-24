#!/usr/bin/env python3
"""[v3 copy: reads meta/panel_v3.tsv, writes depth_basis_v3.tsv and data/ncbi/acc_list_v3.tsv]
Register the panel-v2 species (meta/panel_v2.tsv, in_v1 == no) in meta/species_codes.tsv.

code   : author-released gene sets (meta/extra_local.tsv) keep their prep code;
         otherwise Gen3+sp3 capitalised, collisions -> Gen3+sp2+N (as 12_add_species.py).
source : local (extra_local, proteome already in data/proteomes/), ncbi_annot (RefSeq /
         GenBank annotation), miniprot otherwise.
depth  : depth_median / frac_gt200 from the best available layer, recorded per species in
         meta/depth_basis_v2.tsv (obis species > 02b expansion > sample depth > congeners).
Extra local species absent from panel_v2 (no NCBI genome, e.g. Alvinella pompejana)
are added from EXTRA_ONLY with explicit taxonomy.
Backs up species_codes.tsv first. Run from ~/deepsea_conv. Idempotent (skips species present).
"""
import csv, os, shutil, time, collections

R = os.path.expanduser('~/deepsea_conv')
M = R + '/meta/'
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))

EXTRA_ONLY = {  # species: (phylum, class, order, family, role, accession, depth_basis_note)
    'Alvinella pompejana': ('Annelida', 'Polychaeta', 'Terebellida', 'Alvinellidae', 'deep', '',
                            'hydrothermal vents, East Pacific Rise ~2500 m (Zenodo 11241339 assembly)'),
}

codes_f = M + 'species_codes.tsv'
rows = rd(codes_f)
cols = list(rows[0].keys())
have_sp = {r['species'] for r in rows}
have_code = {r['code'] for r in rows}
extra = {r['species']: r['code'] for r in rd(M + 'extra_local.tsv')}
obis = {r['sp']: r for r in rd(M + 'obis_depth_by_species.tsv')}
dv2 = {r['species']: r for r in rd(M + 'depth_v2.tsv')}
di = {r['species']: r for r in rd(M + 'depth_integrated.tsv') if r['set'] == 'genome'}


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def depth(sp, role):
    o = obis.get(sp)
    if o and (fnum(o['n_depth']) or 0) > 0 and fnum(o['depth_median']) is not None:
        return fnum(o['depth_median']), fnum(o['frac_gt200']), f"OBIS species n_depth={int(fnum(o['n_depth']))}"
    v = dv2.get(sp)
    if v and fnum(v['depth_median']) is not None and fnum(v['frac_gt200']) is not None:
        return fnum(v['depth_median']), fnum(v['frac_gt200']), f"02b {v['source']} n_depth={v['n_depth']}"
    d = di.get(sp, {})
    s = fnum(d.get('sample_max_m'))
    if s is not None:
        return s, 1.0 if s > 200 else 0.0, f"sample depth {d.get('sample_src', '')}"
    g = fnum(d.get('genus_median_m'))
    if g is not None:
        return g, 1.0 if role == 'deep' else 0.0, f"congeners n={d.get('genus_n', '')} ({d.get('final_src', '')})"
    return None, None, 'none'


def mk_code(sp):
    g, s = sp.split()[0], (sp.split() + ['sp'])[1].rstrip('.')
    code = (g[:3] + s[:3]).capitalize()
    i = 1
    while code in have_code:
        code = (g[:3] + s[:2]).capitalize() + str(i)
        i += 1
    return code


todo = [(r['species'], r['role'], r['accession'], r['annot'], r['phylum'], r['class'], r['order'], r['family'])
        for r in rd(M + 'panel_v3.tsv') if r['in_v1'] == 'no']


ANN = {t[0]: t[3] for t in todo}
new, basis = [], []
for sp, role, acc, annot, ph, cl, od, fam in todo:
    if sp in have_sp:
        continue
    if sp in extra:
        code, src = extra[sp], 'local'
        assert os.path.getsize(f'{R}/data/proteomes/{code}.faa') > 0, code
    else:
        code = mk_code(sp)
        src = 'ncbi_annot' if annot in ('refseq', 'genbank') else 'miniprot'
    assert code not in have_code or sp in extra, code
    have_code.add(code)
    have_sp.add(sp)
    med, frac, note = depth(sp, role)
    if sp in EXTRA_ONLY:
        med, frac, note = 2500.0, 1.0, 'literature: ' + EXTRA_ONLY[sp][6]
    new.append(dict(code=code, species=sp, accession=acc, source=src, role=role, phylum=ph,
                    **{'class': cl}, order=od, family=fam,
                    depth_median='' if med is None else round(med, 2),
                    frac_gt200='' if frac is None else round(frac, 3)))
    basis.append([code, sp, role, src, note])

shutil.copy(codes_f, f"{R}/backups/species_codes.tsv.{time.strftime('%Y%m%d_%H%M%S')}")
with open(codes_f, 'a') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n')
    for r in new:
        w.writerow(r)
bf = M + 'depth_basis_v3.tsv'
first = not os.path.exists(bf)
with open(bf, 'a') as f:
    if first:
        f.write('code\tspecies\trole\tsource\tdepth_basis\n')
    for b in basis:
        f.write('\t'.join(map(str, b)) + '\n')
with open(R + '/data/ncbi/acc_list_v3.tsv', 'w') as f:
    for r in new:
        if r['source'] != 'local':
            f.write(r["accession"] + chr(9) + (ANN[r["species"]] or "none") + chr(10))
print('added', len(new))
print(collections.Counter((r['source'], r['role']) for r in new))
print(collections.Counter(b[4].split()[0] for b in basis))
