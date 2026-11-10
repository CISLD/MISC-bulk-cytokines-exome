#!/usr/bin/env python3
"""Tier annotated panel variants per patient and test recurrence.

Tiers (all require qual_pass: DP >= 10, GQ >= 20, het alt fraction 0.2-0.8):
  A  ClinVar 20260615 P/LP (not conflicting, not a swapped-allele lift record), gnomAD AF < 0.05
  B  LoF, gnomAD exome AF < 0.001 and in-house case frequency <= 0.5%
  C  biallelic candidate: homozygous/hemizygous, or >= 2 heterozygous (phase unknown) in one gene;
     LoF or missense/inframe, gnomAD AF < 0.01 and in-house case frequency <= 1%
  D  rare missense/inframe, gnomAD AF < 0.001 and in-house <= 0.5%, in the MIS-C candidate or
     inborn-errors-of-immunity categories of Table S9

usage: wes_tiers.py <annotated.tsv> <out_dir> <n_inhouse_cases>
"""
import csv, collections, sys

ann_path, out_dir, n_cases = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = list(csv.DictReader(open(ann_path), delimiter='\t'))
codes = sorted({r['code'] for r in rows})


def f(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


PLP = {'P', 'LP', 'P/LP'}
IMMUNE_CATS = ('MIS-C', 'Primary Immunodeficiency')
for r in rows:
    r['gaf'] = f(r['gnomad_exome_af'])
    r['inhouse_frq'] = f(r['inhouse_cases']) / n_cases
    r['ok'] = r['qual_pass'] == 'True'
    r['coding'] = r['variant_class'] in ('LoF', 'missense/inframe')
    r['plp'] = r['clinvar_class'] in PLP and r['clinvar_swapped'] != 'true'
    r['rare_dom'] = r['gaf'] < 0.001 and r['inhouse_frq'] <= 0.005
    r['rare_rec'] = r['gaf'] < 0.01 and r['inhouse_frq'] <= 0.01
    r['tiers'] = []

for r in rows:
    if r['ok'] and r['plp'] and r['gaf'] < 0.05:
        r['tiers'].append('A')
    if r['ok'] and r['variant_class'] == 'LoF' and r['rare_dom']:
        r['tiers'].append('B')
    if r['ok'] and r['coding'] and r['rare_dom'] and any(c in r['panel_category'] for c in IMMUNE_CATS) \
            and r['variant_class'] != 'LoF':
        r['tiers'].append('D')

by_gene = collections.defaultdict(list)
for r in rows:
    if r['ok'] and r['coding'] and r['rare_rec']:
        by_gene[(r['code'], r['gene'])].append(r)
for (code, gene), vs in by_gene.items():
    if any(v['zygosity'] == 'hom' for v in vs) or sum(v['zygosity'] == 'het' for v in vs) >= 2:
        for v in vs:
            v['tiers'].append('C')

cols = ['code', 'tiers', 'gene', 'panel_category', 'chrom', 'pos', 'ref', 'alt', 'zygosity', 'dp', 'gq',
        'alt_fraction', 'effect', 'variant_class', 'transcript', 'hgvs_c', 'hgvs_p', 'gnomad_exome_af',
        'inhouse_cases', 'inhouse_hom', 'clinvar_id', 'clinvar_class', 'clinvar_sig', 'clinvar_stars',
        'clinvar_swapped', 'clinvar_bcbio']
with open(f'{out_dir}/candidates_by_patient.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(cols)
    for r in rows:
        if r['tiers']:
            r['tiers'] = ','.join(sorted(set(r['tiers'])))
            w.writerow([r[c] for c in cols])

cand = [r for r in rows if r['tiers']]
with open(f'{out_dir}/summary_by_patient.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['code', 'panel_variants_pass', 'rare_coding_AF<0.01', 'A_clinvar_PLP', 'B_rare_LoF',
                'C_biallelic_variants', 'C_biallelic_genes', 'D_rare_missense_immune'])
    for c in codes:
        pv = [r for r in rows if r['code'] == c and r['ok']]
        cc = [r for r in cand if r['code'] == c]
        w.writerow([c, len(pv), sum(r['coding'] and r['rare_rec'] for r in pv),
                    sum('A' in r['tiers'] for r in cc), sum('B' in r['tiers'] for r in cc),
                    sum('C' in r['tiers'] for r in cc), len({r['gene'] for r in cc if 'C' in r['tiers']}),
                    sum('D' in r['tiers'] for r in cc)])

# recurrence: the manuscript filter (exonic, MAF < 0.01, present in >= 5 patients) and gene-level sharing
site = collections.defaultdict(set)
gene = collections.defaultdict(set)
for r in rows:
    if r['ok'] and r['coding'] and r['rare_rec']:
        site[(r['gene'], r['chrom'], r['pos'], r['ref'], r['alt'], r['hgvs_p'], r['gnomad_exome_af'],
              r['inhouse_cases'])].add(r['code'])
    if r['ok'] and r['coding'] and r['rare_dom']:
        gene[r['gene']].add(r['code'])
with open(f'{out_dir}/recurrence_sites_AF001.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['gene', 'chrom', 'pos', 'ref', 'alt', 'hgvs_p', 'gnomad_exome_af', 'inhouse_cases',
                'n_patients'])
    for k, v in sorted(site.items(), key=lambda kv: -len(kv[1])):
        if len(v) >= 2:
            w.writerow(list(k) + [len(v)])
with open(f'{out_dir}/recurrence_genes_AF0001.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['gene', 'n_patients', 'patients'])
    for g, v in sorted(gene.items(), key=lambda kv: -len(kv[1])):
        if len(v) >= 2:
            w.writerow([g, len(v), ','.join(sorted(v))])

n5 = sum(len(v) >= 5 for v in site.values())
print(f'patients: {len(codes)}; candidate rows: {len(cand)}; '
      f'rare (AF<0.01) coding sites in >=5 patients: {n5}')
