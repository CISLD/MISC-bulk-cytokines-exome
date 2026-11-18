#!/usr/bin/env python3
"""Update the numbers in the BioRender Figure 1 PDF without redrawing it.

v21 (2026-11-18): copy of figure1_biorender_numbers.py that also relabels "Single Cell Multiome" as
"Single-Nucleus Multiome" (the assay profiles nuclei). Same input (original submitted BioRender PDF).

Edits the page content stream in place, using glyphs of the embedded Roboto subset, so artwork, font
and layout are unchanged: cytokine cohort N=11 -> N=8, and the analysed pairs per assay.

usage: figure1_biorender_numbers.py <in.pdf> <out.pdf>
"""
import re
import sys
import zlib

src, dst = sys.argv[1:3]
b = open(src, "rb").read()
objs = {int(m.group(1)): m.group(2) for m in re.finditer(rb"(\d+) 0 obj(.*?)endobj", b, re.S)}

def stream(k):
    h, rest = objs[k].split(b"stream", 1)
    s = rest.lstrip(b"\r\n").rsplit(b"endstream", 1)[0]
    return h, (zlib.decompress(s) if b"FlateDecode" in h else s)

# Regular font (F2): ToUnicode is one bfrange array, CID i -> character
cm = stream(203)[1].decode()
arr = re.search(r"<0000> <[0-9a-f]+> \[(.*)\]", cm).group(1)
glyphs = [bytes.fromhex(x.replace(" ", "")).decode("utf-16-be") for x in re.findall(r"<([0-9a-f ]+)>", arr)]
cid = {g: i for i, g in enumerate(glyphs) if len(g) == 1}

def enc(text):
    missing = [c for c in text if c not in cid]
    if missing:
        raise SystemExit(f"glyphs not in embedded subset: {missing}")
    return "<" + "".join(f"{cid[c]:04x}" for c in text) + ">"

def dec(line):
    return "".join(glyphs[int(h[j:j + 4], 16)] for h in re.findall(r"<([0-9a-f]+)>", line) for j in range(0, len(h), 4))

hdr, content = stream(5)
lines = content.decode("latin1").split("\n")
labels = [i for i, l in enumerate(lines) if l.endswith("TJ")]
text = {i: dec(lines[i]) for i in labels}

def find(t, nth=0):
    hits = [i for i in labels if text[i] == t]
    return hits[nth]

def append(i, extra):
    assert lines[i].endswith(" 0] TJ")
    lines[i] = lines[i][:-len(" 0] TJ")] + " " + enc(extra) + " 0] TJ"

i = find("(N=11)")
lines[i] = "[" + enc("(N=8)") + " 0] TJ"
append(find("Cytokine profiling"), " (8 pairs)")
append(find("Bulk RNA sequencing", 0), " (9 pairs)")
append(find("Bulk RNA sequencing", 1), " (9 pairs)")
append(find("Whole exome sequencing"), " (10 patients)")

# Multiome: second line "(10 pairs)" in the same text block, one line below
i = find("Single Cell Multiome ATAC + Gene Expression")
lines[i] = "[" + enc("Single-Nucleus Multiome ATAC +") + " 0] TJ"
text[i] = "Single-Nucleus Multiome ATAC +"
tm = i - 1
assert lines[tm].endswith("Tm")
lines.insert(i + 1, "1 0 0 -1 0 32 Tm")
lines.insert(i + 2, "[" + enc("Gene Expression (10 pairs)") + " 0] TJ")

new = "\n".join(lines).encode("latin1")
objs[5] = (b"\n<< /Length " + str(len(new)).encode() + b" >>\nstream\n" + new + b"\nendstream\n")

# Rebuild the file with a fresh xref
trailer = re.search(rb"trailer\s*<<(.*?)>>\s*startxref", b, re.S).group(1)
trailer = re.sub(rb"/Size \d+", b"", trailer)
trailer = re.sub(rb"/Prev \d+", b"", trailer)
out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
offs = {}
for k in sorted(objs):
    offs[k] = len(out)
    out += f"{k} 0 obj".encode() + objs[k] + b"endobj\n"
n = max(objs) + 1
x = len(out)
out += f"xref\n0 {n}\n0000000000 65535 f \n".encode()
for k in range(1, n):
    out += (f"{offs[k]:010d} 00000 n \n" if k in offs else "0000000000 65535 f \n").encode()
out += b"trailer\n<< /Size " + str(n).encode() + b" " + trailer.strip() + b" >>\nstartxref\n" + str(x).encode() + b"\n%%EOF\n"
open(dst, "wb").write(out)
print("edited:", ", ".join(sorted({"N=11->N=8", "pairs added"})))
