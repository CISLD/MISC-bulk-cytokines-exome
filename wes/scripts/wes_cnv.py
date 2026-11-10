#!/usr/bin/env python3
"""CNV calls (WesFlow consensus, AnnotSV full-mode rows) overlapping Table S9 genes, per patient.
The AnnotSV header is shifted by one in these files: 'SV_type' holds the copy number and 'user#3'
the DEL/DUP type; AnnotSV ranking was not run (ACMG_class is NA), so support = number of callers.
usage: wes_cnv.py <sample_key.tsv> <CNV-WesFlow dir> <panel_genes.tsv> <out_dir>"""
import csv, glob, os, sys, collections
key, cnvdir, panel_path, out_dir = sys.argv[1:5]
csv.field_size_limit(10**9)
panel = {r[0] for r in csv.reader(open(panel_path), delimiter='\t')}
rows, summ = [], []
for k in csv.DictReader(open(key), delimiter='\t'):
    f = glob.glob(os.path.join(cnvdir, '*', 'results', k['folder'] + '*', '*.anno.tsv'))
    if not f:
        summ.append([k['code'], 'no CNV file', '', '']); continue
    n = 0
    for r in csv.DictReader(open(f[0]), delimiter='\t'):
        if r['Annotation_mode'] != 'full':
            continue
        n += 1
        callers = r['user#1'].split(';'); typ = r['user#3'].split(';')[0]
        genes = [g for g in r['Gene_name'].split(';') if g in panel]
        if genes:
            rows.append([k['code'], r['SV_chrom'], r['SV_start'], r['SV_end'], typ, r['SV_type'].split(';')[0],
                         r['user#1'], len(callers), ';'.join(genes),
                         r.get('B_loss_AFmax', '') if typ == 'DEL' else r.get('B_gain_AFmax', '')])
    mine = [x for x in rows if x[0] == k['code']]
    summ.append([k['code'], n, len(mine), sum(x[7] >= 2 for x in mine)])
with open(f'{out_dir}/cnv_panel_overlaps.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['code', 'chrom', 'start', 'end', 'type', 'copy_number', 'callers', 'n_callers',
                'panel_genes', 'benign_AFmax'])
    w.writerows(rows)
with open(f'{out_dir}/cnv_summary.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['code', 'cnv_calls', 'calls_overlapping_panel_genes', 'of_which_>=2_callers'])
    w.writerows(summ)
