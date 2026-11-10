#!/usr/bin/env python3
"""Trace every value in the 'grouped by cytokine' sheet (the input of Table S2) back to the
raw wells: which matrix (DIL 2x, DIL 10x, serum) and which summary (mean of both wells,
first well, second well) it equals.

usage: cytokines_reconcile.py <Legendplex.xlsx> <wells_long.tsv> <out.tsv>
"""
import csv, re, sys, collections
from cytokines_parse import read_sheets

ALIAS = {  # grouped-sheet header -> analyte name in the raw sheets
    'IL-6': 'IL-6', 'CCL2': 'MCP-1 (CCL2)', 'G-CSF': 'G-CSF', 'IFN-a2': 'IFN-a2', 'IL-2': 'IL-2',
    'CXCL10': 'CXCL10 (IP-10)', 'RANTES': 'RANTES (CCL5)', 'INF-y': 'IFN-y', 'IL-7': 'IL-7',
    'IL-1RA': 'IL-1RA', 'CXCL8': 'CXCL8 (IL-8)', 'TNF-a': 'TNF-a', 'MIP1/CCL3': 'MIP-1α (CCL3)',
    'IL10': 'IL-10', 'IL-13': 'IL-13', 'GM-CSF': 'GM-CSF', 'IL-1β': 'IL-1β', 'IL-5': 'IL-5',
    'IL-2Ra': 'sCD25 (IL-2Ra)', 'IL-4': 'IL-4', 'VEGF': 'VEGF', 'IL-17A': 'IL-17A', 'IL-18': 'IL-18',
    'CCL4': 'CCL4 (MIP-1β)', 'IL-15': 'IL-15', 'IL-12p70': 'IL-12p70', 'APRIL': 'APRIL',
}
SAMPLE = re.compile(r'^S\d+(-\d+)?$')


def num(s):
    try:
        return float(str(s).replace(',', '.'))
    except ValueError:
        return None


def main():
    xlsx, wells_path, dst = sys.argv[1:4]
    rows = read_sheets(xlsx)['grouped by cytokine']
    wells = collections.defaultdict(dict)
    for r in csv.DictReader(open(wells_path), delimiter='\t'):
        wells[(r['analyte'], r['sample'], r['matrix'])][int(r['well'])] = r

    grouped = []  # (analyte, side, sample, cell)
    header = None
    for r in rows:
        names = {i: v for i, v in enumerate(r) if v in ALIAS}
        if names:
            header = names
            continue
        if header is None:
            continue
        for i, name in header.items():
            if i + 1 < len(r) and SAMPLE.match(r[i] if i < len(r) else ''):
                grouped.append((ALIAS[name], 'flare', r[i], r[i + 1]))
            if i + 3 < len(r) and SAMPLE.match(r[i + 2]):
                grouped.append((ALIAS[name], 'remission', r[i + 2], r[i + 3]))

    out = csv.writer(open(dst, 'w'), delimiter='\t', lineterminator='\n')
    out.writerow(['analyte', 'side', 'sample', 'grouped_value', 'source'])
    for analyte, side, sample, cell in grouped:
        g = num(cell)
        src = []
        for m in ('DIL2X', 'DIL10X', 'SERUM'):
            w = wells.get((analyte, sample, m), {})
            v = {k: num(x['value']) for k, x in w.items() if x['flag'] == 'ok'}
            cands = {}
            if 1 in v and 2 in v:
                cands['mean'] = (v[1] + v[2]) / 2
            for k in v:
                cands[f'well{k}'] = v[k]
            if len(v) == 1:
                cands['single_detected_well'] = list(v.values())[0]
            for lab, c in cands.items():
                if g is not None and abs(c - g) <= max(1e-3, 1e-4 * abs(g)):
                    src.append(f'{m}:{lab}')
        if g is None:
            source = f'text:{cell}'
        elif not src:
            source = 'no_match'
        elif any(s.startswith('DIL2X:mean') for s in src):
            source = 'DIL2X:mean'
        else:
            source = src[0] if len({s.split(':')[0] for s in src}) == 1 else '|'.join(sorted(set(src)))
        out.writerow([analyte, side, sample, cell, source])


if __name__ == '__main__':
    main()
