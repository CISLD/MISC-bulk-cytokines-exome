# MIS-C flare vs remission: whole-blood RNA-seq, cytokines, exome/HLA, public cohorts and figures

Analysis code for the article *Interferon-linked transcription and chromatin accessibility extend to
naive T cells in multisystem inflammatory syndrome in children: a paired multi-omic cohort study*
(Jenko Bizjan et al., submitted). The single-nucleus analyses are in the companion
repository [LabenaBioinf/10X_multiome-anlysis](https://github.com/LabenaBioinf/10X_multiome-anlysis)
(pipeline in `scripts/`, article analyses in `article_analyses/`).

**No data are included.** Whole-blood RNA-seq data are in GEO (GSE301727). Exome sequencing data cannot
be made publicly available under Slovenian law and are available from the corresponding author on
request. De-identified sample keys are not part of the repository. Paths are written as `/path/to/...`;
set `WORKDIR` in each job file.

## Contents

| Article item | Script (job file) | What it does |
|---|---|---|
| Table S7, Figure 4A, C | `scripts/bulk/bulk_rerun.R` (`.sb`) | whole-blood RNA-seq, main and validation cohorts (nine pairs each): DESeq2 `~ Patient_ID + Condition`, with and without adjustment for estimated blood composition |
| Figure 5C inputs | `scripts/bulk/figure_data_bulk.R` | per-sample normalised counts |
| Figure 2, Table S2 | `scripts/cyto/cytokines_parse.py`, `cytokines_reconcile.py`, `cytokines.sb`, `cytokines_analysis.R`, `cytokines_pairs.sb`, `cytokines_table_s2_pairs.py`, `cytokines_figure2_pairs.R` | LEGENDplex export parsed per well; plasma 1:2, duplicate wells; eight flare–remission pairs; values below LLOQ set to LLOQ/√2; exact paired Wilcoxon, BH over 26 analytes |
| Exome review, Table S9 | `wes/scripts/` (see below) | per-patient variant review in MIS-C/inborn-error genes, CNV overlap |
| Table S12 | `wes/scripts/wes_hla.sb`, `hla_confirm/scripts/` | HLA class I typing with OptiType 1.3.5 (exome) and T1K 1.0.9 (exome and RNA); concordance |
| Figure 7, Table S11 | `scripts/rev/public_comparators.R` (`.sb`) | effector module (BATF, STAT1, PRDM1, IRF4) in GSE225220, GSE73461, GSE66099 and GSE167028; limma, with composition covariates for whole blood only |
| Figures 2–8, Supplementary Figure S4 | `scripts/figures.sb` | runs `figures_article_v6.R` (Figures 2–5), `make_figures_v6.py` (6A, 7), `figr_figure6_v8.R` (6B–C, S4), `figure6_composite_v6.py` and `figure8_v8.py` |

### Exome steps (`wes/scripts/`)

Input: bcbio-nextgen GRCh37 output per sample, WesFlow CNV output (AnnotSV tables), the Table S9 gene list
and a de-identified `sample_key.tsv` (not included).

| Step | Script | What it does |
|---|---|---|
| 1 | `wes_scan.sb` | per sample: normalise, PASS non-reference genotypes, pairwise genotype discordance; runs steps 2–5 |
| 2 | `wes_qc.py` | coverage, contamination (FREEMIX), sex check |
| 3 | `wes_panel_variants.py` | coding/splice variants in panel genes, DP ≥ 10, GQ ≥ 20, het alt fraction 0.2–0.8 |
| 4 | `wes_annotate.sql` | in-house variant frequencies and ClinVar |
| 5 | `wes_tiers.py` | tiers: ClinVar P/LP, rare LoF, biallelic candidates, rare missense in MIS-C/IEI genes |
| 6 | `wes_gene_callable.py` | callable coding bases per key gene |
| 7 | `wes_cnv.py` | CNV calls overlapping panel genes |
| 8 | `wes_hla.sb` + `optitype_config.ini` | MHC-region and unmapped reads → OptiType (DNA mode) |

## Software

SLURM batch files run each step in an Apptainer container (`apptainer exec`): DESeq2 1.46.0
(Bioconductor container), an R container with edgeR, limma, fgsea and ggplot2, OptiType 1.3.5, T1K 1.0.9,
bcftools/samtools 1.21, and Python 3 with matplotlib 3.5.1 for figures. Versions are given in the
article's Methods.

## Licence and citation

MIT licence (see `LICENSE`). Please cite the article (see `CITATION.cff`).
