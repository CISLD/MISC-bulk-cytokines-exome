#!/usr/bin/env python3
"""Revised Figures 3, 4, 5, 6A and 7 for the MIS-C manuscript, from the patient-level results.

usage: make_figures_v5.py <results dir> <wes/patient_key.tsv> <out dir> [figure ...]
  figures: 3 4 5 6A 7 (default: all whose inputs exist)
Colours: reference palette of the dataviz method (blue #2a78d6 = remission, orange #eb6834 = flare;
diverging blue - grey #f0efec - red for fold changes). Stars: genome-wide FDR; † = panel FDR only.
"""
import csv
import itertools
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

RES, KEY, OUT = sys.argv[1:4]
WANT = set(sys.argv[4:]) or {"3", "4", "5", "6A", "7"}
os.makedirs(OUT, exist_ok=True)

FLARE, REM = "#eb6834", "#2a78d6"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec", "#e8826f", "#b3261e"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False})

def read(*p):
    with open(os.path.join(*p)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

key = {r["file_code"]: r["manuscript_code"] for r in read(KEY) if r["manuscript_code"] not in ("-", "")}
LABEL = {"gdT": "γδ T", "CD4 Naive": "CD4 naive", "CD8 Naive": "CD8 naive", "NK_CD56bright": "NK CD56bright",
         "CD14 Mono": "CD14 mono", "CD16 Mono": "CD16 mono", "NK Proliferating": "NK prolif.",
         "CD8 Proliferating": "CD8 prolif.", "B naive": "B naive"}
lab = lambda c: LABEL.get(c, c)

def stars(f):
    return "***" if f < 0.001 else "**" if f < 0.01 else "*" if f < 0.05 else ""

def star_rule(gw, panel):
    s = stars(float(gw)) if gw not in ("NA", "") else ""
    if not s and panel not in ("NA", ""):
        s = stars(float(panel)) + ("†" if stars(float(panel)) else "")
    return s

def signed_rank_exact(d):
    """Two-sided exact Wilcoxon signed-rank P (zeros dropped, average ranks for ties)."""
    d = [x for x in d if x != 0]
    n = len(d)
    if n == 0:
        return 1.0
    a = sorted(range(n), key=lambda i: abs(d[i]))
    r = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(d[a[j + 1]]) == abs(d[a[i]]):
            j += 1
        for k in range(i, j + 1):
            r[a[k]] = (i + j) / 2 + 1
        i = j + 1
    obs = sum(r[i] for i in range(n) if d[i] > 0)
    mid = sum(r) / 2
    dist = [sum(r[i] for i in range(n) if s[i]) for s in itertools.product([0, 1], repeat=n)]
    return sum(1 for w in dist if abs(w - mid) >= abs(obs - mid) - 1e-9) / len(dist)

def paired_panel(ax, pairs, title, ylab, note=None, logy=False):
    """pairs: list of (patient, remission value, flare value)."""
    for p, r, f in pairs:
        ax.plot([0, 1], [r, f], color="#b9b8b3", lw=0.8, zorder=1)
    ax.scatter([0] * len(pairs), [p[1] for p in pairs], s=16, color=REM, ec="white", lw=0.5, zorder=3)
    ax.scatter([1] * len(pairs), [p[2] for p in pairs], s=16, color=FLARE, ec="white", lw=0.5, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Remission", "Flare"])
    ax.set_xlim(-0.4, 1.4)
    if logy:
        ax.set_yscale("log")
    ax.set_title(title, fontsize=7.8, color=INK, pad=3)
    ax.set_ylabel(ylab, fontsize=7)
    if note:
        ax.text(0.5, 1.0, note, transform=ax.transAxes, ha="center", va="bottom", fontsize=6.4, color=INK2)
    ax.grid(axis="y", color=GRID, lw=0.5)

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=400)
    plt.close(fig)
    print("wrote", name)

# ---------------------------------------------------------------- Figure 3
def fig3():
    counts = read(RES, "a6_composition", "nuclei_counts_per_subset_per_sample.tsv")
    tests = {(r["set"], r["cell_type"]): r for r in read(RES, "a6_composition", "composition_tests_recomputed.tsv")}
    e20 = {r["cell_type"]: r for r in read(RES, "a6_composition", "composition_sensitivity_each20.tsv")}
    samples = [c for c in counts[0] if c != "cell_type"]
    tot = {s: sum(int(r[s]) for r in counts) for s in samples}
    pts = sorted({s.split(".", 1)[1] for s in samples}, key=lambda p: int(key.get(p, "P99")[1:]))
    prop = {r["cell_type"]: {s: int(r[s]) / tot[s] for s in samples} for r in counts}
    cnt = {r["cell_type"]: {s: int(r[s]) for s in samples} for r in counts}
    show = ["NK Proliferating", "B naive", "MAIT", "CD8 TEM", "CD8 Proliferating", "cDC2", "CD4 CTL", "gdT",
            "CD4 TCM", "Plasmablast"]
    um_path = os.path.join(RES, "figure_data", "umap.tsv")
    has_umap = os.path.exists(um_path)
    fig = plt.figure(figsize=(7.2, 9.0 if has_umap else 5.6))
    top = 0.60 if has_umap else 0.97
    if has_umap:
        um = read(RES, "figure_data", "umap.tsv")
        types = sorted({r["cell_type"] for r in um}, key=lambda t: -sum(1 for r in um if r["cell_type"] == t))
        LINEAGE = [("CD4⁺ T", "#2a78d6", ("CD4", "Treg")), ("CD8⁺ T", "#eb6834", ("CD8",)),
                   ("Unconventional T", "#1baf7a", ("MAIT", "gdT", "dnT")), ("NK / ILC", "#eda100", ("NK", "ILC")),
                   ("B / plasmablast", "#8e5bd0", ("B ", "Plasmablast")), ("Myeloid / DC", "#d6457a", ("Mono", "DC", "ASDC")),
                   ("Other", "#9a9994", ())]
        def lineage(t):
            for name, c, keys in LINEAGE:
                if any(k in t + " " for k in keys):
                    return name, c
            return "Other", "#9a9994"
        col = {t: lineage(t)[1] for t in types}
        import statistics
        cen = {t: (statistics.median(float(r["umap1"]) for r in um if r["cell_type"] == t),
                   statistics.median(float(r["umap2"]) for r in um if r["cell_type"] == t)) for t in types}
        big = [t for t in types if sum(1 for r in um if r["cell_type"] == t) >= 800
               and t not in ("Eryth", "Treg", "B intermediate", "B memory")]
        for k, cond in enumerate(["Remission", "Flare"]):
            ax = fig.add_axes([0.05 + k * 0.45, 0.66, 0.40, 0.31])
            phase = lambda c: "Flare" if c.upper() in ("MISC", "FLARE") else "Remission"
            sub = [r for r in um if phase(r["condition"]) == cond]
            import random
            random.seed(1)
            n = min(len(sub), 37000)
            sub = random.sample(sub, n)
            ax.scatter([float(r["umap1"]) for r in sub], [float(r["umap2"]) for r in sub], s=0.15, lw=0,
                       c=[col[r["cell_type"]] for r in sub], rasterized=True)
            ax.set_title(f"{cond} ({n:,} nuclei shown)", fontsize=8, color=INK)
            for t in big:
                ax.text(cen[t][0], cen[t][1], "B cells" if t == "B naive" else lab(t), fontsize=5.2, ha="center", va="center", color=INK,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if k == 0:
                fig.text(0.01, 0.975, "A", fontsize=12, weight="bold", va="top")
        lax = fig.add_axes([0.05, 0.615, 0.9, 0.04]); lax.axis("off")
        h = [plt.Line2D([], [], marker="o", ls="", color=c, ms=5, label=name) for name, c, _ in LINEAGE]
        lax.legend(handles=h, ncol=7, fontsize=6.4, loc="center", handletextpad=0.2, columnspacing=1.0)
    # paired proportions (10 populations) and CD4/CD8 ratio on a 2 x 6 grid
    letter = "B" if has_umap else "A"
    ys, hh = ((0.35, 0.05), 0.18) if has_umap else ((0.59, 0.09), 0.27)
    fig.text(0.01, ys[0] + hh + 0.045, letter, fontsize=12, weight="bold", va="top")
    def slot(k):
        row, colm = divmod(k, 6)
        return fig.add_axes([0.065 + colm * 0.158, ys[row], 0.115, hh])
    for k, ct in enumerate(show):
        ax = slot(k)
        pairs = [(p, 100 * prop[ct][f"remission.{p}"], 100 * prop[ct][f"flare.{p}"]) for p in pts]
        t_all, t20 = tests[("all", ct)], e20.get(ct)
        note = f"FDR {float(t_all['fdr_exact_30']):.3f}\nboth ≥20 nuclei: " + (
            f"{float(t20['fdr_each20_30']):.2f} ({t20['n_pairs_each20']} pairs)" if t20 and t20["testable_each20"] == "True" else "n.t.")
        paired_panel(ax, pairs, lab(ct), "% of PBMC nuclei" if k % 6 == 0 else "")
        ax.text(0.5, -0.36 if has_umap else -0.22, note, transform=ax.transAxes, ha="center", va="top", fontsize=5.6, color=INK2)
    cd4 = ["CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating", "Treg"]
    cd8 = ["CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating"]
    ratio = lambda s: sum(cnt[c][s] for c in cd4 if c in cnt) / max(1, sum(cnt[c][s] for c in cd8 if c in cnt))
    pairs = [(p, ratio(f"remission.{p}"), ratio(f"flare.{p}")) for p in pts]
    d = [f - r for _, r, f in pairs]
    pval = signed_rank_exact(d)
    up = sum(1 for x in d if x > 0)
    ax = slot(10)
    paired_panel(ax, pairs, "CD4⁺/CD8⁺ T ratio", "")
    ax.text(0.5, -0.36 if has_umap else -0.22, f"{up}/10 higher\nexact P = {pval:.3f}", transform=ax.transAxes,
            ha="center", va="top", fontsize=5.6, color=INK2)
    fig.text(0.065 + 5 * 0.158 - 0.01, ys[1] + hh * 0.5,
             "Lines connect samples\nof the same patient.\nExact paired Wilcoxon\ntest, BH over 30\npopulations; second\n"
             "line: only pairs with\n≥20 nuclei in both\nsamples (n.t., < 4).", fontsize=5.8, color=INK2, va="center")
    save(fig, "Figure3_composition")
    return up, pval

# ---------------------------------------------------------------- Figure 4
def fig4():
    cam = read(RES, "pseudobulk", "camera_hallmark_all.tsv")
    sets = [("HALLMARK_INTERFERON_GAMMA_RESPONSE", "IFN-γ response"), ("HALLMARK_INTERFERON_ALPHA_RESPONSE", "IFN-α response"),
            ("HALLMARK_MYC_TARGETS_V1", "MYC targets V1"), ("HALLMARK_MTORC1_SIGNALING", "mTORC1 signalling"),
            ("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Oxidative phosphorylation"), ("HALLMARK_E2F_TARGETS", "E2F targets"),
            ("HALLMARK_G2M_CHECKPOINT", "G2M checkpoint"), ("HALLMARK_IL6_JAK_STAT3_SIGNALING", "IL-6/JAK/STAT3"),
            ("HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response"), ("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF-α via NF-κB"),
            ("HALLMARK_TGF_BETA_SIGNALING", "TGF-β signalling")]
    order = ["CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg", "CD8 Naive", "CD8 TCM", "CD8 TEM", "MAIT", "gdT", "dnT",
             "NK", "NK_CD56bright", "B naive", "B intermediate", "B memory", "Plasmablast", "CD14 Mono", "CD16 Mono", "pDC"]
    idx = {(r["set"], r["cell_type"]): r for r in cam}
    order = [c for c in order if any((s, c) in idx for s, _ in sets)]
    M = [[0.0] * len(order) for _ in sets]
    for i, (s, _) in enumerate(sets):
        for j, c in enumerate(order):
            r = idx.get((s, c))
            if r:
                v = min(10.0, -math.log10(max(float(r["FDR"]), 1e-300)))
                M[i][j] = v if r["Direction"] == "Up" else -v
    kg = read(RES, "pseudobulk_followup", "key_gene_panel_pseudobulk.tsv")
    genes = ["STAT1", "PRDM1", "BATF", "IRF4", "IL6ST", "STAT3", "STAT4", "TBX21", "IFNG", "IFNG-AS1", "GZMB", "PRF1",
             "CTLA4", "LAG3", "HAVCR2", "PDCD1", "TOX", "NR4A2"]
    cts = ["CD8 TEM", "NK", "gdT", "CD4 Naive", "CD8 Naive", "CD14 Mono"]
    gidx = {(r["cell_type"], r["gene"]): r for r in kg}
    fig = plt.figure(figsize=(7.2, 8.4))
    ax = fig.add_axes([0.24, 0.60, 0.66, 0.36])
    im = ax.imshow(M, cmap=DIV, norm=TwoSlopeNorm(0, -10, 10), aspect="auto")
    for i in range(len(sets)):
        for j, c in enumerate(order):
            r = idx.get((sets[i][0], c))
            if r and float(r["FDR"]) < 0.05:
                ax.text(j, i, "•", ha="center", va="center", fontsize=6, color=INK)
            if not r:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc="white", ec=GRID, hatch="////", lw=0))
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([lab(c) for c in order], rotation=50, ha="right")
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels([n for _, n in sets])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("signed −log₁₀ FDR (camera; capped at 10)", fontsize=6.5)
    cb.outline.set_visible(False)
    fig.text(0.01, 0.975, "A", fontsize=12, weight="bold", va="top")
    fig.text(0.24, 0.965, "Hallmark gene sets, flare vs remission (paired pseudobulk; • FDR < 0.05)", fontsize=7.5, color=INK2)
    ax2 = fig.add_axes([0.24, 0.05, 0.50, 0.44])
    L = [[float(gidx[(c, g)]["logFC"]) if (c, g) in gidx else float("nan") for c in cts] for g in genes]
    im2 = ax2.imshow(L, cmap=DIV, norm=TwoSlopeNorm(0, -3, 3), aspect="auto")
    for i, g in enumerate(genes):
        for j, c in enumerate(cts):
            r = gidx.get((c, g))
            if r:
                s = star_rule(r["FDR_genomewide"], r["FDR_panel"])
                v = float(r["logFC"])
                ax2.text(j, i, f"{v:+.2f}{s}", ha="center", va="center", fontsize=5.6,
                         color="white" if abs(v) > 2.0 else INK)
            else:
                ax2.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc="white", ec=GRID, hatch="////", lw=0))
    ax2.set_xticks(range(len(cts)))
    ax2.set_xticklabels([lab(c) for c in cts], rotation=30, ha="right")
    ax2.set_yticks(range(len(genes)))
    ax2.set_yticklabels(genes, style="italic")
    ax2.tick_params(length=0)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    cb2 = fig.colorbar(im2, ax=ax2, fraction=0.04, pad=0.03)
    cb2.set_label("log₂ fold change, flare vs remission (capped ±3)", fontsize=6.5)
    cb2.outline.set_visible(False)
    fig.text(0.01, 0.515, "B", fontsize=12, weight="bold", va="top")
    fig.text(0.24, 0.505, "Key genes (paired pseudobulk; * genome-wide FDR, † panel FDR only; *<0.05 **<0.01 ***<0.001)",
             fontsize=7.0, color=INK2)
    save(fig, "Figure4_programmes")

# ---------------------------------------------------------------- Figure 5
def fig5():
    kg = read(RES, "pseudobulk_followup", "key_gene_panel_pseudobulk.tsv")
    order = ["CD8 TEM", "gdT", "NK", "NK_CD56bright", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg", "CD8 Naive", "CD8 TCM",
             "MAIT", "B naive", "B memory", "CD14 Mono"]
    gidx = {(r["cell_type"], r["gene"]): r for r in kg}
    order = [c for c in order if (c, "IFNG-AS1") in gidx or (c, "IFNG") in gidx]
    fig = plt.figure(figsize=(7.2, 8.0))
    ax = fig.add_axes([0.10, 0.62, 0.86, 0.33])
    x = range(len(order))
    for k, (g, colr, dx) in enumerate([("IFNG", FLARE, -0.2), ("IFNG-AS1", REM, 0.2)]):
        for j, c in enumerate(order):
            r = gidx.get((c, g))
            if not r:
                continue
            v = float(r["logFC"]); s = star_rule(r["FDR_genomewide"], r["FDR_panel"])
            ax.bar(j + dx, v, width=0.36, color=colr, edgecolor="white", lw=0.6, zorder=2,
                   alpha=1.0 if s else 0.45)
            if s:
                ax.text(j + dx, v + (0.08 if v >= 0 else -0.08), s, ha="center", va="bottom" if v >= 0 else "top",
                        fontsize=5.5, color=INK)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xticks(list(x))
    ax.set_xticklabels([lab(c) for c in order], rotation=35, ha="right")
    ax.set_ylabel("log₂ fold change, flare vs remission")
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=FLARE, label="IFNG"),
                       plt.Rectangle((0, 0), 1, 1, color=REM, label="IFNG-AS1"),
                       plt.Rectangle((0, 0), 1, 1, color="#bbbbbb", label="not significant (faded)")],
              loc="lower left", fontsize=6.5, ncol=3)
    fig.text(0.01, 0.975, "A", fontsize=12, weight="bold", va="top")
    fig.text(0.10, 0.965, "Paired pseudobulk by population (* genome-wide FDR, † panel FDR only)", fontsize=7.5, color=INK2)
    pb_path = os.path.join(RES, "figure_data", "pb_ifng.tsv")
    if os.path.exists(pb_path):
        pb = read(RES, "figure_data", "pb_ifng.tsv")
        fig.text(0.01, 0.50, "B", fontsize=12, weight="bold", va="top")
        for k, ct in enumerate(["CD8 TEM", "gdT", "NK"]):
            for m, g in enumerate(["IFNG", "IFNG-AS1"]):
                ax = fig.add_axes([0.08 + k * 0.32, 0.29 - m * 0.0 + (0 if m == 0 else -0.0), 0.12, 0.17])
                ax.set_position([0.08 + k * 0.32 + m * 0.15, 0.29, 0.11, 0.17])
                rows = [r for r in pb if r["cell_type"] == ct and r["gene"] == g and int(r["n_cells"]) >= 10]
                byp = {}
                for r in rows:
                    byp.setdefault(r["patient"], {})["flare" if r["condition"].upper() in ("MISC", "FLARE") else "rem"] = float(r["log2cpm"])
                pairs = [(p, v["rem"], v["flare"]) for p, v in byp.items() if "rem" in v and "flare" in v]
                gr = gidx.get((ct, g))
                note = f"{float(gr['logFC']):+.2f}{star_rule(gr['FDR_genomewide'], gr['FDR_panel'])}" if gr else ""
                paired_panel(ax, pairs, f"{lab(ct)}: {g}", "log₂ CPM" if (k == 0 and m == 0) else "")
                ax.text(0.5, -0.32, note, transform=ax.transAxes, ha="center", fontsize=6, color=INK2)
    bk_path = os.path.join(RES, "figure_data", "bulk_ifng.tsv")
    if os.path.exists(bk_path):
        bk = read(RES, "figure_data", "bulk_ifng.tsv")
        fig.text(0.01, 0.20, "C", fontsize=12, weight="bold", va="top")
        stats = {("main", "IFNG"): "+0.73 (padj 0.19)", ("main", "IFNG-AS1"): "−2.14 (padj 2×10⁻³²)",
                 ("validation", "IFNG"): "+0.09 (padj 0.91)", ("validation", "IFNG-AS1"): "−1.52 (padj 4×10⁻⁷)"}
        for k, (coh, g) in enumerate(itertools.product(["main", "validation"], ["IFNG", "IFNG-AS1"])):
            ax = fig.add_axes([0.08 + k * 0.235, 0.06, 0.15, 0.12])
            byp = {}
            for r in bk:
                if r["cohort"] == coh and r["gene"] == g:
                    byp.setdefault(r["pair"], {})[r["condition"]] = float(r["log2norm"])
            pairs = [(p, v["remission"], v["flare"]) for p, v in byp.items() if len(v) == 2]
            paired_panel(ax, pairs, f"{'Main' if coh == 'main' else 'Validation'}: {g}",
                         "log₂ normalised count" if k == 0 else "")
            ax.text(0.5, -0.42, stats[(coh, g)], transform=ax.transAxes, ha="center", fontsize=5.8, color=INK2)
    save(fig, "Figure5_IFNG_IFNGAS1")

# ---------------------------------------------------------------- Figure 6A
def fig6a():
    p = os.path.join(RES, "figure_data", "chromvar_stat.tsv")
    if not os.path.exists(p):
        print("Figure 6A: waiting for", p); return
    cv = read(RES, "figure_data", "chromvar_stat.tsv")
    # tests from the pre-declared extension run, which reproduces results/gdt_chromvar exactly for the original populations
    st = {(r["cell_type"], r["feature"]): r for r in read(RES, "chromvar_naive_extension", "motif_deviation_flare_vs_remission.tsv")}
    condcol = "cond" if "cond" in cv[0] else "condition"
    cts = ("CD8 TEM", "gdT", "CD4 TCM", "CD4 Naive", "CD8 Naive", "Treg", "B naive")
    fig = plt.figure(figsize=(7.2, 5.6))
    for i, m in enumerate(("STAT1", "STAT4")):
        for k, ct in enumerate(cts):
            ax = fig.add_axes([0.075 + k * 0.133, 0.62 - i * 0.48, 0.085, 0.25])
            byp = {}
            for r in cv:
                if r["cell_type"] == ct and r["motif"] == m:
                    byp.setdefault(r["patient"], {})["f" if r[condcol].lower().startswith(("fl", "mis")) else "r"] = float(r["z"])
            pairs = [(p, v["r"], v["f"]) for p, v in byp.items() if len(v) == 2]
            s = st.get((ct, m))
            paired_panel(ax, pairs, f"{lab(ct)}\n{m} motif", "chromVAR deviation z" if k == 0 else "")
            ax.tick_params(axis="x", labelsize=5.6)
            ax.tick_params(axis="y", labelsize=5.8)
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
            if s:
                ax.text(0.5, -0.24, f"FDR {float(s['FDR']):.3f}\n{s['patients_flare_higher']}/{s['n_pairs']} up",
                        transform=ax.transAxes, ha="center", va="top", fontsize=5.8, color=INK2)
    fig.text(0.01, 0.97, "A", fontsize=12, weight="bold", va="top")
    save(fig, "Figure6A_STAT_motifs")

# ---------------------------------------------------------------- Figure 7
def fig7():
    rows = read(RES, "public_comparators", "comparator_module_score.tsv")
    nice = {"MIS.C.acute": "MIS-C, acute", "MIS.C.3m.": "MIS-C, 3 months", "COVID.19.acute": "Paediatric COVID-19, acute",
            "Kawasaki": "Kawasaki disease", "DefiniteViral": "Viral infection", "DefiniteBacterial": "Bacterial infection",
            "Inflammatory": "Non-infectious inflammation", "SIRS": "SIRS", "Sepsis": "Sepsis", "SepticShock": "Septic shock",
            "MIS.C": "MIS-C", "KD": "Kawasaki disease", "MIS.C_MYO": "MIS-C with myocarditis"}
    items = []
    for r in rows:
        case = r["contrast"].split(" vs ")[0]
        items.append((r["dataset"], nice.get(case, case), r["composition_adjusted"] == "TRUE", float(r["delta"]),
                      float(r["P"]), r["n_cases"], r["n_controls"]))
    order_ds = ["GSE225220 whole blood", "GSE73461 whole blood", "GSE66099 whole blood", "GSE167028 PBMC"]
    labels, ylist = [], []
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    fig.subplots_adjust(left=0.36, right=0.97, top=0.93, bottom=0.10)
    y = 0
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
                _, _, _, d, p, nc, nk = hit[0]
                yy = y + (0.17 if not adj else -0.17)
                ax.plot([0, d], [yy, yy], color="#cfcec9", lw=1.0, zorder=1)
                ax.scatter([d], [yy], s=28, zorder=3, color=FLARE if adj else "white", ec=FLARE, lw=1.2)
                ax.text(1.75, yy, "p<0.0001" if p < 0.0001 else f"p = {p:.2g}", fontsize=6, va="center", color=INK2)
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
    ax.set_xlim(-0.7, 2.05)
    ax.set_xticks([-0.5, 0, 0.5, 1.0, 1.5])
    ax.set_xlabel("Module score difference vs controls (mean z of BATF, STAT1, PRDM1, IRF4)")
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc="white", mec=FLARE, mew=1.2, ms=5, label="unadjusted"),
                       plt.Line2D([], [], marker="o", ls="", color=FLARE, ms=5, label="adjusted for blood composition")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.5)
    save(fig, "Figure7_public_cohorts")

if "3" in WANT:
    print("CD4/CD8 ratio:", fig3())
if "4" in WANT:
    fig4()
if "5" in WANT:
    fig5()
if "6A" in WANT:
    fig6a()
if "7" in WANT:
    fig7()
