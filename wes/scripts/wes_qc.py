#!/usr/bin/env python3
"""Per-sample QC from the bcbio run: mean target coverage, VerifyBamID FREEMIX, peddy sex check.
Output carries only the de-identified code.

usage: wes_qc.py <sample_key.tsv> <run_dir> <out.tsv>
"""
import csv, glob, os, re, sys

key_path, run, out_path = sys.argv[1:4]
final = glob.glob(os.path.join(run, '*', 'final'))[0]
summary = glob.glob(os.path.join(final, '2024-*', 'project-summary.yaml'))[0]

# project-summary.yaml: per-sample metrics blocks; read Name and Avg_coverage pairs without a yaml lib
cov, name, avg = {}, None, None
for line in open(summary):
    m = re.match(r'\s+Avg_coverage:\s*([\d.]+)', line)
    if m:
        avg = float(m.group(1))
    m = re.match(r'\s+Name:\s*(\S+)', line)
    if m:
        name = m.group(1)
        if avg is not None:
            cov[name] = avg
        avg = None

rows = []
for r in csv.DictReader(open(key_path), delimiter='\t'):
    folder = r['folder']
    d = os.path.join(final, folder, 'qc')
    freemix = ''
    for f in glob.glob(os.path.join(d, 'contamination', '*.selfSM')):
        rec = list(csv.DictReader(open(f), delimiter='\t'))
        if rec:
            freemix = rec[0].get('FREEMIX(alpha)', '')
    pred = ''
    for f in glob.glob(os.path.join(d, 'peddy', '*sex_check.csv')):
        rec = list(csv.DictReader(open(f)))
        if rec:
            pred = rec[0]['predicted_sex']
    reported = {'F': 'female', 'M': 'male'}.get(folder.rsplit('-', 1)[-1], '')
    rows.append([r['code'], cov.get(folder, ''), freemix, reported, pred,
                 'ok' if reported == pred else 'MISMATCH'])

with open(out_path, 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['code', 'mean_target_coverage', 'freemix', 'sex_label', 'sex_predicted', 'sex_check'])
    w.writerows(rows)
print(open(out_path).read())
