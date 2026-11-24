#!/usr/bin/env python3
"""Transcriptome layer: deep (>200 m) species that have no genome but a TSA assembly.

Selection : depth_tsa.tsv depth_class_v2 == deep, species not already in species_codes.tsv;
            the TSA master with the most contigs (meta/tsa_masters.tsv slen) is used.
Download  : NCBI wgs_aux  <PREFIX4>01.1.fsa_nt.gz
ORFs      : longest ORF per contig on both strands, >= 100 codons (TransDecoder.LongOrfs
            default); complete (ATG..stop) ORFs preferred over 5'/3'-partial ones of < 1.2x
            length. Redundancy: mmseqs2 easy-cluster at 95 % identity / 90 % coverage,
            longest member kept (isoform/allele collapse).
Output    : data/proteomes/CODE.faa, data/cds/CODE.fna, meta/idmap/CODE.tsv,
            species_codes.tsv rows with source = tsa (backup first), meta/tsa_prep_report.tsv.
Usage     : 28_tsa_prep.py [threads]      (run from ~/deepsea_conv)
"""
import csv, gzip, os, re, shutil, subprocess, sys, time, urllib.request

R = os.path.expanduser('~/deepsea_conv')
M = R + '/meta/'
T = sys.argv[1] if len(sys.argv) > 1 else '16'
E = R + '/envs/dc/bin'
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
MINAA = 100

bases = 'TCAG'
aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: aas[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases)
         for k, c in enumerate(bases)}
RC = str.maketrans('ACGTN', 'TGCAN')


def tr(s):
    return ''.join(CODON.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def best_orf(seq):
    """Return (cds, pep, complete) of the longest ORF on either strand, or None."""
    seq = seq.upper()
    cands = []
    for strand in (seq, seq.translate(RC)[::-1]):
        for f in range(3):
            p = tr(strand[f:])
            pos = 0
            for seg in p.split('*'):
                end = pos + len(seg)          # index of stop codon (or end) in aa
                has_stop = end < len(p)
                m = seg.find('M')
                if m >= 0 and len(seg) - m >= MINAA:
                    a, comp = pos + m, has_stop
                    cands.append((len(seg) - m, comp, f, a, end, has_stop, strand))
                if len(seg) >= MINAA and (pos == 0 or not has_stop):  # 5' or 3' partial
                    cands.append((len(seg), False, f, pos, end, has_stop, strand))
                pos = end + 1
    if not cands:
        return None
    longest = max(c[0] for c in cands)
    comp = [c for c in cands if c[1] and c[0] * 1.2 >= longest]
    L, c_ok, f, a, end, has_stop, strand = max(comp or cands, key=lambda c: c[0])
    cds = strand[f + 3 * a: f + 3 * end + (3 if has_stop else 0)]
    pep = tr(cds).rstrip('*')
    if 'X' in pep and pep.count('X') > len(pep) * 0.05:
        return None
    return cds, pep, c_ok


def read_fa(path):
    h, buf = None, []
    with gzip.open(path, 'rt') as f:
        for l in f:
            if l.startswith('>'):
                if h:
                    yield h, ''.join(buf)
                h, buf = l[1:].split()[0], []
            else:
                buf.append(l.strip())
    if h:
        yield h, ''.join(buf)


def mk_code(sp, have):
    g, s = sp.split()[0], (sp.split() + ['sp'])[1].rstrip('.')
    code, i = (g[:3] + s[:3]).capitalize(), 1
    while code in have:
        code, i = (g[:3] + s[:2]).capitalize() + str(i), i + 1
    return code


codes = rd(M + 'species_codes.tsv')
have_sp, have_code = {r['species'] for r in codes}, {r['code'] for r in codes}
masters = {}
for r in rd(M + 'tsa_masters.tsv'):
    if r['organism'] not in masters or int(r['slen']) > int(masters[r['organism']]['slen']):
        masters[r['organism']] = r
sel = [r for r in rd(M + 'depth_tsa.tsv') if r['depth_class_v2'] == 'deep'
       and r['species'] not in have_sp and r['valid_name'] not in have_sp and r['species'] in masters]
# TSA_V3=1: also take TSA species called deep only by the integrated layers (02d: congeners /
# sample depth), which depth_tsa.tsv (occurrence + range only) leaves as nodata/eury.
DI_TSA = {r['species']: r for r in rd(M + 'depth_integrated.tsv') if r['set'] == 'tsa'}
if os.environ.get('TSA_V3'):
    got = {r['species'] for r in sel}
    sel += [r for r in rd(M + 'depth_tsa.tsv') if r['species'] not in got
            and DI_TSA.get(r['species'], {}).get('final_class') == 'deep'
            and r['species'] not in have_sp and r['valid_name'] not in have_sp and r['species'] in masters]


def ncbi_lineage(taxid):
    """rank -> name from NCBI taxonomy (same source as the genome panel)."""
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml'
    # curl, not urllib: urllib from this server stalls on NCBI's IPv6 endpoints
    x = subprocess.run(['curl', '-s', '-4', '--max-time', '60', url], capture_output=True, text=True).stdout
    time.sleep(0.4)
    d = dict((rk, nm) for nm, rk in re.findall(r'<Taxon>\s*<TaxId>\d+</TaxId>\s*<ScientificName>([^<]+)'
                                              r'</ScientificName>\s*<Rank>([^<]+)</Rank>', x))
    d['_lineage'] = (re.search(r'<Lineage>([^<]*)</Lineage>', x) or [None, ''])[1]
    return d


# WoRMS name matching can hit homonyms (e.g. Episinus = spider genus in NCBI): keep a TSA
# species only if NCBI places its taxid in the WoRMS phylum; taxonomy columns come from NCBI.
LIN, ok_sel = {}, []
for r in sel:
    lin = ncbi_lineage(masters[r['species']]['taxid'])
    if r['w_phylum'] not in lin['_lineage'] or 'Metazoa' not in lin['_lineage']:
        print('HOMONYM_SKIP', r['species'], r['w_phylum'], '| NCBI:', lin['_lineage'][-80:], flush=True)
        continue
    LIN[r['species']] = lin
    ok_sel.append(r)
sel = ok_sel
print('TSA deep species to add:', len(sel), flush=True)

os.makedirs(R + '/data/tsa', exist_ok=True)
rep_f = M + 'tsa_prep_report.tsv'
rep = open(rep_f, 'a')
if os.path.getsize(rep_f) == 0:
    rep.write('code\tspecies\ttsa\tn_contigs\tn_orf\tn_complete\tn_kept\n')
new = []
for r in sel:
    sp, m = r['species'], masters[r['species']]
    acc = m['accession'][:4]
    code = mk_code(sp, have_code)
    have_code.add(code)
    fa = f'{R}/data/tsa/tsa.{acc}.fsa_nt.gz'
    if not os.path.exists(fa):
        # GenBank TSA layout: genbank/tsa/<A>/tsa.<ABCD>.<vol>.fsa_nt.gz (vol 1..n); gz members concatenate
        vol = 1
        with open(fa + '.part', 'wb') as out:
            while True:
                url = f'https://ftp.ncbi.nlm.nih.gov/genbank/tsa/{acc[0]}/tsa.{acc}.{vol}.fsa_nt.gz'
                x = subprocess.run(['curl', '-s', '-4', '-f', '--retry', '3', url], capture_output=True)
                if x.returncode:
                    break
                out.write(x.stdout)
                vol += 1
        if vol == 1:
            print('DLFAIL', sp, url, flush=True)
            os.remove(fa + '.part')
            continue
        os.rename(fa + '.part', fa)
    w = f'{R}/work/tsa/{code}'
    os.makedirs(w, exist_ok=True)
    n = n_orf = n_comp = 0
    orfs = {}
    for h, s in read_fa(fa):
        n += 1
        o = best_orf(s)
        if o:
            n_orf += 1
            n_comp += o[2]
            orfs[h] = o
    with open(f'{w}/orf.faa', 'w') as f:
        for h, (c, p, ok) in orfs.items():
            f.write(f'>{h}\n{p}\n')
    subprocess.run([f'{E}/mmseqs', 'easy-cluster', f'{w}/orf.faa', f'{w}/clu', f'{w}/tmp',
                    '--min-seq-id', '0.95', '-c', '0.9', '--cov-mode', '1', '--threads', T],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    reps = {}
    for l in open(f'{w}/clu_cluster.tsv'):
        rp, mem = l.split()
        if rp not in reps or len(orfs[mem][1]) > len(orfs[reps[rp]][1]):
            reps[rp] = mem
    keep = sorted(set(reps.values()), key=lambda h: (-len(orfs[h][1]), h))
    with open(f'{R}/data/proteomes/{code}.faa', 'w') as fp, open(f'{R}/data/cds/{code}.fna', 'w') as fc, \
            open(f'{M}idmap/{code}.tsv', 'w') as fm:
        for i, h in enumerate(keep, 1):
            nid = f'{code}_{i:06d}'
            c, p, ok = orfs[h]
            if c[-3:] in ('TAA', 'TAG', 'TGA'):
                c = c[:-3]
            fp.write(f'>{nid}\n{p}\n')
            fc.write(f'>{nid}\n{c}\n')
            fm.write(f'{nid}\t{h}\t{h}\n')
    shutil.rmtree(f'{w}/tmp', ignore_errors=True)
    row = [code, sp, m['accession'], n, n_orf, n_comp, len(keep)]
    rep.write('\t'.join(map(str, row)) + '\n')
    rep.flush()
    print(*row, flush=True)
    med = r['depth_median'] or DI_TSA.get(sp, {}).get('genus_median_m') or DI_TSA.get(sp, {}).get('sample_max_m') or (str((float(r['range_sh']) + float(r['range_dp'])) / 2)
                                if r['range_sh'] and r['range_dp'] else '')
    new.append(dict(code=code, species=sp, accession=m['accession'], source='tsa', role='deep',
                    phylum=LIN[sp].get('phylum', r['w_phylum']), **{'class': LIN[sp].get('class', r['w_class'])},
                    order=LIN[sp].get('order', r['w_order']), family=LIN[sp].get('family', r['w_family']), depth_median=med, frac_gt200=r['frac_gt200'] or '1.0'))

if new:
    shutil.copy(M + 'species_codes.tsv', f"{R}/backups/species_codes.tsv.{time.strftime('%Y%m%d_%H%M%S')}")
    with open(M + 'species_codes.tsv', 'a') as f:
        wr = csv.DictWriter(f, list(codes[0].keys()), delimiter='\t', lineterminator='\n')
        for x in new:
            wr.writerow(x)
print('added', len(new))
