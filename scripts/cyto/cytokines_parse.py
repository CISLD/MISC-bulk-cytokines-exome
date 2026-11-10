#!/usr/bin/env python3
"""Parse the LEGENDplex export (sheets '14plex' and '13plex') into one row per well.

The export has flare samples on the left and their paired remission sample on the same row
on the right. Each sample block is: a label row (sample ID and/or matrix such as 'DIL 2x',
'DIL 10X', 'serum'), a second well on the next row, then an 'average' row.

Output columns: panel, analyte, side, sample, pair_flare, matrix, well, raw, value, flag, limit
  flag  = ok | below (raw '<x') | above (raw '>x') | missing
  limit = the detection limit printed in a '<x' or '>x' cell (pg/mL, dilution-corrected)

usage: cytokines_parse.py <Legendplex.xlsx> <out.tsv>
"""
import re, sys, zipfile
from xml.etree import ElementTree as ET

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
RNS = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
SAMPLE = re.compile(r'^S\d+(-\d+)?$')
MATRIX = re.compile(r'^(dil\s*\d+x|serum)$', re.I)


def read_sheets(path):
    z = zipfile.ZipFile(path)
    ss = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).iter(NS + 'si'):
            ss.append(''.join(t.text or '' for t in si.iter(NS + 't')))
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    out = {}
    for sh in wb.iter(NS + 'sheet'):
        tgt = rels[sh.get(RNS + 'id')].lstrip('/')
        tgt = tgt if tgt.startswith('xl/') else 'xl/' + tgt
        rows = []
        for row in ET.fromstring(z.read(tgt)).iter(NS + 'row'):
            cells = {}
            for c in row.iter(NS + 'c'):
                ref = re.match(r'[A-Z]+', c.get('r')).group()
                col = 0
                for ch in ref:
                    col = col * 26 + ord(ch) - 64
                v, t = c.find(NS + 'v'), c.get('t')
                if t == 's' and v is not None:
                    val = ss[int(v.text)]
                elif t == 'inlineStr':
                    val = ''.join(x.text or '' for x in c.iter(NS + 't'))
                else:
                    val = v.text if v is not None else ''
                cells[col - 1] = (val or '').strip()
            n = max(cells) + 1 if cells else 0
            rows.append([cells.get(i, '') for i in range(n)])
        out[sh.get('name')] = rows
    return out


def analyte_name(h):
    """'MCP-1 (CCL2) (A4)' -> 'MCP-1 (CCL2)'; drops the bead-region suffix."""
    return re.sub(r'\s*\([AB]\d+\)\s*$', '', h).strip()


def parse_value(s):
    s = s.strip()
    if not s:
        return '', 'missing', ''
    m = re.match(r'^([<>])\s*([\d.,]+)$', s)
    if m:
        lim = float(m.group(2).replace(',', '.'))
        return '', 'below' if m.group(1) == '<' else 'above', lim
    try:
        return float(s.replace(',', '.')), 'ok', ''
    except ValueError:
        return '', 'missing', ''


def parse_panel(rows, panel):
    hdr = rows[1]
    cols = [i for i, h in enumerate(hdr) if re.search(r'\([AB]\d+\)\s*$', h)]
    # two analyte blocks: flare (left) and remission (right)
    split = next(i for i in range(1, len(cols)) if cols[i] != cols[i - 1] + 1)
    flare_cols, rem_cols = cols[:split], cols[split:]
    rem_id_col = rem_cols[0] - 1
    out = []
    state = {'flare': [None, None, 0], 'remission': [None, None, 0]}
    for r in rows[2:]:
        r = r + [''] * (max(cols) + 2 - len(r))
        for side, id_col, acols in (('flare', 0, flare_cols), ('remission', rem_id_col, rem_cols)):
            st = state[side]
            lab0 = r[id_col]
            lab1 = r[1] if side == 'flare' else ''
            if side == 'flare' and SAMPLE.match(lab0):
                st[:] = [lab0, lab1 or 'DIL 2x', 0]
            elif side == 'flare' and MATRIX.match(lab1):
                st[1], st[2] = lab1, 0
            elif side == 'remission' and SAMPLE.match(lab0):
                st[:] = [lab0, 'DIL 2x', 0]
            elif side == 'remission' and MATRIX.match(lab0):
                st[1], st[2] = lab0, 0
            if (lab1 if side == 'flare' else lab0).lower().startswith('average'):
                continue
            vals = [r[c] for c in acols]
            if st[0] is None or not any(vals):
                continue
            st[2] += 1
            if st[2] > 2:  # every block has two wells; a third row is an unlabelled average
                continue
            for c, v in zip(acols, vals):
                value, flag, lim = parse_value(v)
                out.append([panel, analyte_name(hdr[c]), side, st[0], '', re.sub(r'(?i)dil\s*', 'DIL ', st[1]).upper().replace('DIL ', 'DIL'),
                            st[2], v, value, flag, lim])
    # pairing: remission sample on the same block row as its flare sample
    pair = {}
    for r in rows[2:]:
        r = r + [''] * (rem_id_col + 1 - len(r))
        if SAMPLE.match(r[0]) and SAMPLE.match(r[rem_id_col]):
            pair[r[rem_id_col]] = r[0]
            pair[r[0]] = r[0]
    for o in out:
        o[4] = pair.get(o[3], '')
    return out


def main():
    src, dst = sys.argv[1:3]
    sheets = read_sheets(src)
    rows = parse_panel(sheets['14plex'], 'panel1_14plex') + parse_panel(sheets['13plex'], 'panel2_13plex')
    with open(dst, 'w') as fh:
        fh.write('\t'.join(['panel', 'analyte', 'side', 'sample', 'pair_flare', 'matrix', 'well',
                            'raw', 'value', 'flag', 'limit']) + '\n')
        for o in rows:
            fh.write('\t'.join('' if x is None else str(x) for x in o) + '\n')
    print(f'{len(rows)} wells written to {dst}')


if __name__ == '__main__':
    main()
