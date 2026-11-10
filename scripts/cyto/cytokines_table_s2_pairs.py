#!/usr/bin/env python3
"""Revised Table S2, 8 complete pairs only (unpaired flare samples excluded), from
cytokines_primary_LODsqrt2.tsv (+ the drop-S6 sensitivity run).
usage: cytokines_table_s2.py <primary.tsv> <sensitivity_primary.tsv> <out.tsv>"""
import csv, math, sys
prim, sens, out = sys.argv[1:4]
S = {r['analyte']: r for r in csv.DictReader(open(sens), delimiter='\t')}
def g(x, d=1):
    try: v = float(x)
    except ValueError: return 'NA'
    return f'{v:,.{d}f}'
def p(x):
    try: v = float(x)
    except ValueError: return 'NA'
    return f'{v:.3f}' if v >= 0.001 else f'{v:.1e}'
rows = list(csv.DictReader(open(prim), delimiter='\t'))
rows.sort(key=lambda r: (r['testable'] != 'TRUE', float(r['p_wilcoxon_exact']), r['analyte']))
with open(out, 'w') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['Analyte', 'LLOQ at 2x (pg/mL)', 'Flare n (% < LLOQ)', 'Remission n (% < LLOQ)',
                'Flare median [IQR]', 'Remission median [IQR]', 'Pairs (n)',
                'Pairs up / down', 'Median paired fold change (95% CI)', 'P (exact paired Wilcoxon)',
                'P adj. (BH)', 'P adj. without pair S6', 'Note'])
    for r in rows:
        note = []
        if r['testable'] != 'TRUE':
            note.append(f"not quantifiable: above upper limit in {r['samples_above_uloq']} of "
                        f"{int(r['n_flare']) + int(r['n_remission'])} samples")
        if r['heavily_censored'] == 'TRUE':
            note.append('≥50% below LLOQ in one group; fold change depends on LLOQ substitution')
        fc = lambda x: g(2 ** float(x), 2) if x not in ('NA', '') else 'NA'
        w.writerow([r['analyte'], g(r['lloq_2x'], 2) if r['lloq_2x'] != 'NA' else '—',
                    f"{r['n_flare']} ({g(r['pct_below_flare'])})",
                    f"{r['n_remission']} ({g(r['pct_below_remission'])})",
                    f"{g(r['median_flare_paired'])} [{g(r['q1_flare_paired'])}–{g(r['q3_flare_paired'])}]",
                    f"{g(r['median_remission'])} [{g(r['q1_remission'])}–{g(r['q3_remission'])}]",
                    r['n_pairs'], f"{r['pairs_up']} / {r['pairs_down']}",
                    f"{fc(r['median_log2FC'])} ({fc(r['log2FC_ci_low'])}–{fc(r['log2FC_ci_high'])})",
                    p(r['p_wilcoxon_exact']) if r['testable'] == 'TRUE' else 'not tested',
                    p(r['p_adj_BH']), p(S[r['analyte']]['p_adj_BH']), '; '.join(note)])
print(open(out).read())
