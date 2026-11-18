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
| Whole-blood composition | `scripts/bulk/lm22_nnls_composition.R` (`.sb`) | neutrophil and lymphocyte fractions per sample: LM22 signature matrix (Newman et al. 2015) fitted to Salmon gene TPM by non-negative least squares (Lawson–Hanson, written out in the script; not CIBERSORT). `LM22.txt` is not included: download it from the CIBERSORTx website under its terms |
| Table S7, Figure 4C (whole blood) | `scripts/bulk/bulk_rerun.R` (`.sb`) | whole-blood RNA-seq, main and validation cohorts (nine pairs each): DESeq2 `~ Patient_ID + Condition`, and with the LM22 neutrophil and lymphocyte fractions as covariates (`~ Patient_ID + neut + lymph + Condition`) |
| Figure 4 | `scripts/rev/figure4_camera_v11.R` (run from `scripts/figures.sb`) | Hallmark gene sets with camera: whole blood on paired edgeR models `~ patient + condition`, genes filtered on the flare and remission groups (main and validation cohorts), next to the pseudobulk camera results of the companion repository; writes the source data `camera_hallmark_bulk_pseudobulk.tsv`. Panel C: key genes, whole-blood DESeq2 and paired pseudobulk with genome-wide FDR (genes below the primary expression filter take the genome-wide FDR of the alternative full-design filter, triangles, or are not tested) |
| Figure 5C inputs | `scripts/bulk/figure_data_bulk.R` | per-sample normalised counts |
| Figure 2, Table S2 | `scripts/cyto/cytokines_parse.py`, `cytokines_reconcile.py`, `cytokines.sb`, `cytokines_analysis.R`, `cytokines_pairs.sb`, `cytokines_table_s2_pairs.py`, `cytokines_figure2_pairs.R` | LEGENDplex export parsed per well; plasma 1:2, duplicate wells; eight flare–remission pairs; values below LLOQ set to LLOQ/√2; exact paired Wilcoxon, BH over 26 analytes |
| Exome review, Table S9 | `wes/scripts/` (see below) | per-patient variant review in MIS-C/inborn-error genes, CNV overlap |
| Table S12 | `wes/scripts/wes_hla.sb`, `hla_confirm/scripts/` | HLA class I typing with OptiType 1.3.5 (exome) and T1K 1.0.9 (exome and RNA); concordance |
| Figure 7, Table S11 | `scripts/rev/public_comparators.R` (`.sb`) | effector module (BATF, STAT1, PRDM1, IRF4) in GSE225220, GSE73461, GSE66099 and GSE167028; limma, with composition covariates for whole blood only; writes estimates, p values, moderated-t SE, df and 95% confidence intervals |
| Figures 1–8, Supplementary Figures S4–S5 | `scripts/figures.sb` | runs `figures_article_v11.R` (Figures 2, 5), `figure3_v8.R` (3), `figure4_camera_v11.R` (4), `make_figure6_v11.py` (6, STAT motifs), `figure7_v21.py` (7, with 95% confidence intervals), `figure8_v21.py` (8), `figr_figure6_v16.R` (S4 FigR heatmaps, S5 TBX21 targets) and `figure1_v21.py` (Figure 1 is a BioRender drawing; the script only edits labels and the analysed pairs in its PDF, which is not included). Key-gene significance in Figures 4C, 5 and 8 is the genome-wide FDR only. The targeted gene panel only selects the genes shown |

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
(Bioconductor container), an R container with edgeR, limma, msigdbr and ggplot2, OptiType 1.3.5, T1K 1.0.9,
bcftools/samtools 1.21, and Python 3 with matplotlib 3.5.1 for figures. Versions are given in the
article's Methods.

## Licence and citation

MIT licence (see `LICENSE`). Please cite the article (see `CITATION.cff`). The code used in the article is release v1.1, archived on Zenodo (all versions: https://doi.org/10.5281/zenodo.24120588).
