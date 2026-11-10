"""Link WES codes to the flare (M) and remission (K) bulk RNA-seq FASTQs of the same child by
surname + first initial (WES folder vs FASTQ name). The PAXgene codes in the single_cell tab are
not unique (same remission code on several rows), so they are not used.
usage: rna_key.py <wes/sample_key.tsv> <wes/patient_key.tsv> <fastq_dir> <out.tsv>
Prints codes only."""
import csv, os, sys, unicodedata

def norm(s):
    return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode().upper()

wkey, pkey, fq, dst = sys.argv[1:5]
fastqs = sorted(f for f in os.listdir(fq) if f.endswith('_R1.fastq.gz'))
parsed = []
for f in fastqs:
    t = norm(f.replace('_UMI_R1.fastq.gz', '')).split('_')[0].split('-')
    cond = {'M': 'flare', 'K': 'remission'}.get(t[-1], '?')
    parsed.append((f, set(t[2:-1]), cond, '-'.join(t[:2])))
wes = {r['code']: r for r in csv.DictReader(open(wkey), delimiter='\t')}
rows = []
for r in csv.DictReader(open(pkey), delimiter='\t'):
    w = norm(wes[r['wes_code']]['folder']).split('-')   # lab, yy, initial, surname, wes, sex
    initial, surname = w[2], w[3]
    for cond in ('flare', 'remission'):
        hits = [p for p in parsed if p[2] == cond and surname in p[1] and initial in p[1]]
        if not hits:   # double surnames: FASTQ carries surname + initial of the second surname
            hits = [p for p in parsed if p[2] == cond and surname in p[1]]
        rows.append([r['wes_code'], r['file_code'], r['manuscript_code'], cond,
                     hits[0][3] if len(hits) == 1 else '', len(hits),
                     os.path.join(fq, hits[0][0]) if len(hits) == 1 else ''])
with open(dst, 'w') as h:
    h.write('wes_code\tfile_code\tmanuscript_code\tcondition\tpax_code\tn_match\tr1\n')
    for o in rows:
        h.write('\t'.join(map(str, o)) + '\n')
for o in rows:
    print(*o[:6], sep='\t')
