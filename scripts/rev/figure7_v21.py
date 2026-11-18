#!/usr/bin/env python3
"""Figure 7 v21 (2026-11-18): the effector module in public cohorts, now with 95% confidence intervals.

Effect estimates with 95% confidence intervals (moderated-t SE and df written by public_comparators.R).

usage: figure7_v21.py <results/public_comparators/comparator_module_score.tsv> <out dir>
"""
import csv
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
FLARE = "#eb6834"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False})
with open(src) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
nice = {"MIS.C.acute": "MIS-C, acute", "MIS.C.3m.": "MIS-C, 3 months", "COVID.19.acute": "Paediatric COVID-19, acute",
        "Kawasaki": "Kawasaki disease", "DefiniteViral": "Viral infection", "DefiniteBacterial": "Bacterial infection",
        "Inflammatory": "Non-infectious inflammation", "SIRS": "SIRS", "Sepsis": "Sepsis", "SepticShock": "Septic shock",
        "MIS.C": "MIS-C", "KD": "Kawasaki disease", "MIS.C_MYO": "MIS-C with myocarditis"}
items = []
for r in rows:
    case = r["contrast"].split(" vs ")[0]
    items.append((r["dataset"], nice.get(case, case), r["composition_adjusted"] == "TRUE", float(r["delta"]), float(r["P"]),
                  r["n_cases"], r["n_controls"], float(r["CI95_low"]), float(r["CI95_high"])))
order_ds = ["GSE225220 whole blood", "GSE73461 whole blood", "GSE66099 whole blood", "GSE167028 PBMC"]
labels, ylist = [], []
fig, ax = plt.subplots(figsize=(7.2, 5.4))
fig.subplots_adjust(left=0.36, right=0.97, top=0.93, bottom=0.10)
y = 0
XP = 2.85
for ds in order_ds:
    cases = []
    for it in items:
        if it[0] == ds and it[1] not in cases:
            cases.append(it[1])
    labels.append(ds); ylist.append(y)
    y -= 1
    for c in cases:
        first = [it for it in items if it[0] == ds and it[1] == c][0]
        for adj in (False, True):
            hit = [it for it in items if it[0] == ds and it[1] == c and it[2] == adj]
            if not hit:
                continue
            _, _, _, d, p, nc, nk, lo, hi = hit[0]
            yy = y + (0.17 if not adj else -0.17)
            ax.plot([lo, hi], [yy, yy], color=FLARE if adj else MUTED, lw=1.0, zorder=1, solid_capstyle="butt")
            ax.scatter([d], [yy], s=28, zorder=3, color=FLARE if adj else "white", ec=FLARE, lw=1.2)
            ax.text(XP, yy, "p<0.0001" if p < 0.0001 else f"p = {p:.2g}", fontsize=6, va="center", color=INK2)
        labels.append(f"{c} (n = {first[5]} vs {first[6]})"); ylist.append(y)
        y -= 1
    y -= 0.4
ax.axvline(0, color=MUTED, lw=0.6)
ax.set_ylim(y + 0.3, 0.6)
ax.set_yticks(ylist)
ax.set_yticklabels(labels, fontsize=6.6)
for t, l in zip(ax.get_yticklabels(), labels):
    if l in order_ds:
        t.set_fontweight("bold"); t.set_color(INK); t.set_fontsize(7.2)
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_xlim(-1.1, 3.3)
ax.set_xticks([-1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5])
ax.set_xlabel("Module score difference vs controls (mean z of BATF, STAT1, PRDM1, IRF4), 95% CI")
ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="-", color=MUTED, mfc="white", mec=FLARE, mew=1.2, ms=5, lw=1.0, label="unadjusted"),
                   plt.Line2D([], [], marker="o", ls="-", color=FLARE, ms=5, lw=1.0, label="adjusted for blood composition")],
          loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.5)
for ext, kw in (("pdf", {}), ("png", {"dpi": 400})):
    fig.savefig(os.path.join(OUT, f"Figure7_public_cohorts.{ext}"), **kw)
print("Figure 7 written to", OUT)
