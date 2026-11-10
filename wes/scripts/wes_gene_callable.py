#!/usr/bin/env python3
"""Fraction of coding bases (RefSeq, longest NM_ transcript per gene) inside each sample's bcbio
callable.bed (CALLABLE intervals), for a list of key genes.
usage: wes_gene_callable.py <sample_key.tsv> <run_final_dir> <hg19_refGene.txt> <out.tsv> GENE ..."""
import bisect, collections, csv, glob, os, sys
key, final, refgene, out, *genes = sys.argv[1:]
want = set(genes); tx = {}
for f in csv.reader(open(refgene), delimiter='\t'):
    name, chrom, g = f[1], f[2].replace('chr', ''), f[12]
    if g not in want or not name.startswith('NM_') or '_' in chrom:
        continue
    cs, ce = int(f[6]), int(f[7])
    ex = [(max(int(s), cs), min(int(e), ce)) for s, e in zip(f[9].rstrip(',').split(','), f[10].rstrip(',').split(','))]
    ex = [(s, e) for s, e in ex if e > s]
    L = sum(e - s for s, e in ex)
    if g not in tx or L > tx[g][2]:
        tx[g] = (chrom, ex, L)
rows = []
for r in csv.DictReader(open(key), delimiter='\t'):
    bed = glob.glob(os.path.join(final, r['folder'], '*callable.bed'))[0]
    iv = collections.defaultdict(list)
    for l in open(bed):
        c, s, e, st = l.rstrip('\n').split('\t')[:4]
        if st == 'CALLABLE':
            iv[c].append((int(s), int(e)))
    for c in iv:  # sort and merge overlapping intervals
        m = []
        for s0, e0 in sorted(iv[c]):
            if m and s0 <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e0))
            else:
                m.append((s0, e0))
        iv[c] = m
    starts = {c: [s for s, _ in v] for c, v in iv.items()}
    for g in genes:
        if g not in tx:
            rows.append([r['code'], g, 'NA', 'NA']); continue
        chrom, ex, L = tx[g]; cov = 0
        for s, e in ex:
            v = iv.get(chrom, []); i = max(0, bisect.bisect_right(starts.get(chrom, []), s) - 1)
            while i < len(v) and v[i][0] < e:
                cov += max(0, min(e, v[i][1]) - max(s, v[i][0])); i += 1
        rows.append([r['code'], g, L, f'{cov / L:.3f}'])
with open(out, 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['code', 'gene', 'coding_bp', 'fraction_callable']); w.writerows(rows)
