#!/usr/bin/env python3
"""Figure 8, v8: figure8_v6.py with a new bottom panel for naive CD4+ and
CD8+ T cells (RNA: paired pseudobulk; chromatin: STAT1/STAT4 chromVAR motif deviation, pre-declared extension in
results/chromvar_naive_extension), plus Treg and naive B context notes. The effector part is unchanged, drawn at the
same scale in the upper axes.

Values: paired pseudobulk log2FC flare vs remission (results/pseudobulk_followup/
key_gene_panel_pseudobulk.tsv). Stars from genome-wide FDR; where that is >= 0.05 or not computed
(low expression), from the panel FDR, marked †.

usage: figure8_v8.py <key_gene_panel_pseudobulk.tsv> <motif_deviation_flare_vs_remission.tsv> <out prefix>
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch

src, motif_src, prefix = sys.argv[1:4]
MOT = {}
with open(motif_src) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        MOT[(r["cell_type"], r["feature"])] = r
PB = {}
with open(src) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        PB[(r["cell_type"], r["gene"])] = r

def stars(f):
    return "***" if f < 0.001 else "**" if f < 0.01 else "*" if f < 0.05 else ""

def val(ct, gene, label=None):
    r = PB[(ct, gene)]
    fc = float(r["logFC"])
    gw = r["FDR_genomewide"]
    s = stars(float(gw)) if gw not in ("NA", "") else ""
    if not s:
        s = stars(float(r["FDR_panel"]))
        s = s + "†" if s else " (ns)"
    return f"{label or gene} {fc:+.2f}{s}"

TXT, GREY, FADE = "#222222", "#666666", "#aaaaaa"
RED, BLUE, GREEN, ORANGE = "#b03a2e", "#1f4e79", "#2e7d32", "#b9770e"
# Same scale as v6 (100 x 123 units on 8.5 x 10.4 in): upper axes = v6 drawing above y 14; lower axes = naive panel,
# summary box and legend (LOW units high).
U = 10.4 / 123
TOP0, LOW = 14, 49
H = (123 - TOP0 + LOW) * U
fig = plt.figure(figsize=(8.5, H))
ax = fig.add_axes([0, LOW * U / H, 1, (123 - TOP0) * U / H])
ax.set_xlim(0, 100)
ax.set_ylim(TOP0, 123)
ax.axis("off")
axb = fig.add_axes([0, 0, 1, LOW * U / H])
axb.set_xlim(0, 100)
axb.set_ylim(0, LOW)
axb.axis("off")

def box(x, y, w, h, text, fc="white", ec=TXT, ls="-", fs=7, color=TXT, weight="normal"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.3,rounding_size=0.8",
                                fc=fc, ec=ec, ls=ls, lw=0.9, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=color, weight=weight,
            linespacing=1.25, zorder=4)

def arrow(x1, y1, x2, y2, color=TXT, lw=1.1, rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9, color=color,
                                 lw=lw, ls=(0, (3, 2)), connectionstyle=f"arc3,rad={rad}", zorder=2))

def cell(x, y, w, h, fc, ec):
    """Cell body and nucleus; returns the y of the nucleus top."""
    ax.add_patch(Ellipse((x, y), w, h, fc=fc, ec=ec, lw=1.2, zorder=1))
    ax.add_patch(Ellipse((x, y - h * 0.217), w * 0.85, h * 0.5, fc="white", ec=ec, lw=0.8, zorder=1))

def header(x, ytop, name, ec, sub=None, w=20):
    box(x, ytop + 4.5, w, 3.2, name, ec=ec, fs=8.5, weight="bold")
    if sub:
        ax.text(x, ytop + 1.0, sub, ha="center", va="center", fontsize=6.3, color=BLUE, style="italic",
                linespacing=1.2)

def signalling(x, y, ct, stats, il12_box=None):
    """Receptor line (IL6ST, IL-12R), STAT boxes, dashed receptor→STAT arrows."""
    ax.text(x - 6.5, y + 11, val(ct, "IL6ST") + " ↑", fontsize=6.3, color=GREEN, ha="center")
    ax.text(x + 8, y + 11, "IL-12R", fontsize=6.3, color=GREY, ha="center")
    xs = [x - 13, x - 4] if len(stats) == 2 else [x - 9]
    for sx, s in zip(xs, stats):
        box(sx, y + 6, 7.5, 2.4, s, ec=GREEN, fs=7)
        arrow(x - 6.5, y + 10.2, sx, y + 7.6, color=GREEN)
    if il12_box:
        box(x + 9, y + 6, 15, 2.4, il12_box, ec=BLUE, fs=6.5, color=BLUE, weight="bold")
    else:
        box(x + 9, y + 6, 7.5, 2.4, "STAT4", ec=GREY, fs=7)
    arrow(x + 8, y + 10.2, x + 9, y + 7.6, color=GREY)

# Title


# Upstream signals
ax.text(50, 120.5, "Upstream signals (plasma, 8 flare–remission pairs)", ha="center", fontsize=7.5,
        color=GREY, style="italic")
box(27, 116, 13, 4, "IL-6\nplasma ↑", fc="#f3e5f5", ec="#6a1b9a", fs=7.5, weight="bold")
box(42, 116, 13, 4, "IL-18\nplasma ↑", fc="#f3e5f5", ec="#6a1b9a", fs=7.5, weight="bold")
box(58, 116, 15, 4, "IL-27\nhypothetical, not measured", ec=FADE, ls="--", fs=6.5, color=FADE)
box(75, 116, 15, 4, "IL-12p70 / IL-15\nnot raised in plasma", ec=GREY, fs=6.5, color=GREY)
ax.plot([27, 27, 42, 42], [113.6, 111.5, 111.5, 113.6], color=TXT, lw=0.9)
ax.text(34.5, 112.5, "source unspecified", ha="center", va="center", fontsize=6, color=GREY, style="italic")

# Cells
cx, cy = 22, 78          # CD8 TEM
gx, gy = 78, 72          # γδ T
nx, ny = 50, 37          # NK
cell(cx, cy, 38, 30, "#eef6ee", GREEN)
cell(gx, gy, 38, 30, "#fbeeee", RED)
cell(nx, ny, 38, 28, "#eaf2fb", BLUE)

# Ligands to cells (inferred)
arrow(34.5, 111.5, cx + 1, cy + 22.2, rad=0.1)
arrow(34.5, 111.5, gx - 9, gy + 20.4, rad=-0.2)
arrow(34.5, 111.5, nx, ny + 19.8, rad=-0.05)
arrow(58, 113.6, gx - 4, gy + 20.4, color=FADE)

# CD8 TEM
header(cx, cy + 15, "CD8⁺ TEM", GREEN, "T-bet-associated effector chromatin (flare and remission)", w=20)
signalling(cx, cy, "CD8 TEM", ["STAT1"])
ax.text(cx, cy + 2.6, val("CD8 TEM", "STAT4") + "  ·  STAT-family motifs ↑ (10/10)", ha="center",
        fontsize=5.6, color=GREY)
box(cx, cy - 1.8, 17, 2.2, "TBX21  RNA unchanged", fc="#fff8e1", ec=ORANGE, fs=6.3)
for i, g in enumerate(["BATF", "STAT1", "PRDM1", "IRF4"]):
    ax.text(cx - 12, cy - 5.0 - i * 2.0, val("CD8 TEM", g), fontsize=6.3)
ax.text(cx + 4, cy - 6.0, val("CD8 TEM", "IFNG") + " ↑", fontsize=6.3, color=RED, weight="bold")
ax.text(cx + 4, cy - 8.2, val("CD8 TEM", "IFNG-AS1") + " ↓", fontsize=6.3, color=BLUE, weight="bold")
arrow(cx - 1, cy - 7.5, cx + 3.5, cy - 6.6, color=GREEN)
box(cx, cy - 21, 34, 7.5,
    "Effector and checkpoint output\n"
    f'{val("CD8 TEM", "CTLA4")}, {val("CD8 TEM", "LAG3")}\n'
    f'{val("CD8 TEM", "HAVCR2")}, {val("CD8 TEM", "PDCD1")}\n'
    f'{val("CD8 TEM", "GZMB")}, {val("CD8 TEM", "PRF1")}\n'
    "TOX and NR4A2/NR4A3 not increased", ec=BLUE, fs=6.2)

# γδ T
header(gx, gy + 15, "γδ T cell", RED, w=18)
signalling(gx, gy, "gdT", ["STAT1", "STAT3"])
ax.text(gx, gy + 2.6, val("gdT", "STAT4") + "  ·  STAT-family motifs ↑ (10/10)", ha="center",
        fontsize=5.6, color=GREY)
for i, g in enumerate(["STAT1", "STAT3", "PRDM1", "BATF", "JUNB"]):
    lab = val("gdT", g) + (" (low expr.)" if g == "BATF" else "")
    ax.text(gx - 11.5, gy - 2.2 - i * 1.95, lab, fontsize=6.0)
ax.text(gx + 3, gy - 5.0, val("gdT", "IFNG") + " ↑", fontsize=6.3, color=RED, weight="bold")
ax.text(gx + 3, gy - 7.2, val("gdT", "IFNG-AS1") + " ↓", fontsize=6.3, color=BLUE, weight="bold")
box(gx, gy - 19, 30, 2.8, "NFKB1 ↓ / KLF3 ↑ module associations (exploratory)", ec=FADE, ls="--",
    fs=6.0, color=GREY)

# NK
header(nx, ny + 14, "NK cell", BLUE, w=14)
signalling(nx, ny, "NK", ["STAT1"], il12_box=val("NK", "STAT4") + " ↓")
box(nx, ny + 2.4, 33, 2.2, "T-box motif accessibility ↑ (8/9), with most other motif families",
    fc="#fff8e1", ec=ORANGE, fs=5.8)
for i, g in enumerate(["STAT1", "PRDM1", "BATF"]):
    ax.text(nx - 12, ny - 4.0 - i * 2.1, val("NK", g), fontsize=6.3)
ax.text(nx + 3.5, ny - 5.0, val("NK", "IFNG") + " ↑", fontsize=6.3, color=RED, weight="bold")
ax.text(nx + 3.5, ny - 7.1, val("NK", "IFNG-AS1"), fontsize=6.3, color=GREY)
box(nx, ny - 19.5, 30, 4.2,
    "Cytotoxic output\n"
    f'{val("NK", "GZMB")}, {val("NK", "PRF1")}, {val("NK", "NKG7")}', ec=BLUE, fs=6.2)

# Monocyte and secreted IFN-γ (inferred)
mx, my = 86, 101
ax.add_patch(Ellipse((mx, my), 22, 11, fc="#fdf2e9", ec=ORANGE, lw=1.2, zorder=1))
ax.text(mx, my + 3.0, "CD14⁺ monocyte", ha="center", fontsize=8, weight="bold")
ax.text(mx, my + 0.6, "IFN-γ-responsive", ha="center", fontsize=6.8)
ax.text(mx, my - 1.6, val("CD14 Mono", "STAT1") + " ↑", ha="center", fontsize=6.8)
ax.text(mx, my - 3.6, "PD-L1 ↑ (cell-level only)", ha="center", fontsize=6.3, color=GREY)
px, py = 53.5, 75
ax.add_patch(Ellipse((px, py), 8, 4, fc="#fdecea", ec=RED, lw=1.0, zorder=3))
ax.text(px, py, "IFN-γ", ha="center", va="center", fontsize=7.5, color=RED, weight="bold", zorder=4)
arrow(cx + 15, cy - 4, px - 4, py - 0.8, color=RED, rad=-0.15)
arrow(gx - 17, gy - 1, px + 4, py - 0.3, color=RED, rad=0.15)
arrow(nx + 8, ny + 13.5, px + 1.5, py - 2.1, color=RED, rad=0.25)
arrow(px + 1.5, py + 2.0, mx - 11, my - 1, color=RED, rad=-0.3)
arrow(mx - 11, my + 1, px - 1.5, py + 2.0, color=FADE, rad=0.35)
ax.text(60, 95, "feedback?", fontsize=6.3, color=FADE, ha="center", style="italic")

# ---------------------------------------------------------------- naive T cells (lower axes)
ax = axb


def motif(ct, tf):
    m = MOT[(ct, tf)]
    s = stars(float(m["FDR"])) or " (ns)"
    return f"{tf} {float(m['delta_z']):+.2f}{s}", f"{m['patients_flare_higher']}/{m['n_pairs']}"


ax.plot([4, 96], [48.3, 48.3], color=FADE, lw=0.6)
ax.text(50, 46.2, "Naive T cells: the interferon response is seen in both RNA and chromatin",
        ha="center", va="center", fontsize=8.5, weight="bold", color=TXT)
NAIVE = [("CD4 Naive", "Naive CD4⁺ T cell", 26), ("CD8 Naive", "Naive CD8⁺ T cell", 74)]
yc = 30
for ct, name, x in NAIVE:
    cell(x, yc, 40, 24, "#f1eef8", "#5e35b1")
    ax.text(x, yc + 8.6, name, ha="center", va="center", fontsize=8.5, weight="bold", zorder=4)
    (m1, n1), (m4, n4) = motif(ct, "STAT1"), motif(ct, "STAT4")
    box(x, yc + 4.4, 33, 3.4, f"STAT-family motif accessibility ↑\nSTAT1 motif {m1.split()[1]}, STAT4 motif {m4.split()[1]}"
        f" (higher in {n1} patients)", fc="#fff8e1", ec=ORANGE, fs=5.8)
    for i, g in enumerate(["STAT1", "GBP1", "SOCS1"]):
        ax.text(x - 13, yc - 1.9 - i * 2.0, val(ct, g) + " ↑", fontsize=6.3)
    for i, g in enumerate(["SOCS3", "IL6ST"]):
        ax.text(x + 1.5, yc - 1.9 - i * 2.0, val(ct, g) + " ↑", fontsize=6.3)
    ax.text(x + 1.5, yc - 5.9, val(ct, "IFNG"), fontsize=6.3, color=GREY)
    ax.text(x, yc - 8.6, val(ct, "STAT4", "STAT4 RNA"), ha="center", fontsize=5.8, color=GREY)
# IFN-γ (raised in plasma) acting on naive T cells: inferred
ax.add_patch(Ellipse((50, 38.5), 9, 4.4, fc="#fdecea", ec=RED, lw=1.0, zorder=3))
ax.text(50, 38.5, "IFN-γ", ha="center", va="center", fontsize=7.5, color=RED, weight="bold", zorder=4)
ax.text(50, 35.0, "plasma ↑\n(7/7 pairs)", ha="center", va="center", fontsize=5.8, color=GREY, linespacing=1.15)
arrow(45.5, 38.5, 40.5, 37.0, color=RED, rad=0.1)
arrow(54.5, 38.5, 59.5, 37.0, color=RED, rad=-0.1)
# Context: regulatory T and naive B cells
tm1, tn = motif("Treg", "STAT1")
tm4, _ = motif("Treg", "STAT4")
box(26, 15.6, 40, 3.4, f"Regulatory T cells: same pattern\n{val('Treg', 'STAT1')} RNA; STAT1 motif {tm1.split()[1]}, "
    f"STAT4 motif {tm4.split()[1]} ({tn})", ec="#5e35b1", fs=5.8)
bm1, bn = motif("B naive", "STAT1")
bm4, _ = motif("B naive", "STAT4")
box(74, 15.6, 40, 3.4, f"Naive B cells: RNA only\n{val('B naive', 'STAT1')}, {val('B naive', 'GBP1')} RNA; "
    f"STAT motifs not significant", ec=GREY, fs=5.8, color=GREY)

# Bottom box and legend
box(50, 9.9, 88, 4.6,
    "STAT1 and IL6ST ↑ in CD8⁺ TEM, NK, γδ T and naive T cells · STAT-family motif accessibility ↑ in effector and naive T cells\n"
    "IFNG-AS1 ↓ in CD8⁺ TEM and γδ T cells · upstream signal(s) unresolved",
    fc="#fdf2e9", ec=RED, fs=7.0, weight="bold")
ax.text(5, 5.2, "↑ / ↓ with value: measured change", fontsize=6.3, va="center")
ax.plot([33, 39], [5.2, 5.2], color=TXT, lw=1.1, ls=(0, (3, 2)))
ax.text(40, 5.2, "arrows: inferred / not measured (receptor signalling; STAT phosphorylation not measured)",
        fontsize=6.3, va="center")
ax.text(5, 3.1, "RNA: log₂FC flare vs remission, paired pseudobulk; * FDR < 0.05, ** < 0.01, *** < 0.001 "
        "(genome-wide); † FDR within the targeted (post hoc) gene panel; ns, not significant",
        fontsize=5.8, va="center", color=GREY)
ax.text(5, 1.2, "Motifs: chromVAR deviation z difference, flare vs remission (ATAC pseudobulk); stars as above, FDR within cell type; "
        "STAT1, STAT3, STAT4 and STAT5 motifs are near-identical",
        fontsize=5.8, va="center", color=GREY)

fig.savefig(prefix + ".pdf")
fig.savefig(prefix + ".png", dpi=400)
