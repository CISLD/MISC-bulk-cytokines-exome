# Article versions of Figures 2-5: same data and plot types as figures_original_style_v5.R; no titles in the graphs,
# uppercase panel letters, no box outlines, British spelling, 300 dpi. Options are not rebuilt.
# Figures 2-5 in the format of the submitted figures (plot types, colours, titles and labels as in
# /path/to/MISC/MISC/R_scripts: cytokines.Rmd, UMAP_batch.R, Cell_type_porpotions.R, gsea_hallmark_final.R,
# Fig_4C_DotPlot_v3.R, expressions_all_patients_IFNG_IFNG-AS1.R), redrawn from the corrected data:
#   Fig 2  8 flare-remission pairs, plasma 1:2, LLOQ/sqrt2; exact paired Wilcoxon, BH over 26 analytes
#   Fig 3  Harmony UMAP (umap_harmony, as submitted); proportions per sample; exact paired Wilcoxon, BH over 30
#   Fig 4  fgsea (as submitted) on corrected rankings: bulk DESeq2 Wald statistic (validation pairing fixed),
#          pseudobulk signed sqrt(F) of the paired edgeR QL test; key genes from bulk DESeq2 and paired pseudobulk
#   Fig 5  paired pseudobulk IFNG / IFNG-AS1; per-patient pseudobulk CPM; whole-blood normalised counts
# Main versions (option A) are written to <out>; alternatives to <out>/options.
# usage: Rscript figures_original_style_v5.R <project dir> <harmony-UMAP seurat .rds> <out dir> <R lib with fgsea>
for (l in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(Sys.setlocale("LC_CTYPE", l))) break
args <- commandArgs(trailingOnly = TRUE)
P <- args[1]; OBJ <- args[2]; OUT <- args[3]; .libPaths(c(args[4], .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(gghalves); library(dplyr); library(tidyr)
  library(edgeR); library(limma); library(msigdbr); library(fgsea)
})
OPT <- file.path(OUT, "options"); dir.create(OPT, recursive = TRUE, showWarnings = FALSE)
R <- file.path(P, "results")
rd <- function(...) read.delim(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
exact_paired <- function(d) {
  d <- d[abs(d) > 1e-12]; n <- length(d); if (!n) return(1)
  r <- rank(abs(d)); V <- sum(r[d > 0]); s <- as.matrix(expand.grid(rep(list(c(0, 1)), n)))
  mean(abs(s %*% r - sum(r) / 2) >= abs(V - sum(r) / 2) - 1e-9)
}
save2 <- function(p, file, w, h, dpi = 150) {
  if (grepl("/options/", file)) return(invisible(NULL))
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 300)
}

## ======================================================================= Figure 2
sv <- rd(R, "cytokines_pairs_only", "sample_values_DIL2X.tsv")
st <- rd(R, "cytokines_pairs_only", "cytokines_primary_LODsqrt2.tsv")
lab <- c("IL-6" = "IL-6", "IFN-y" = "IFN-γ", "VEGF" = "VEGF", "MCP-1 (CCL2)" = "CCL2", "IL-7" = "IL-7",
         "GM-CSF" = "GM-CSF", "IL-17A" = "IL-17A", "G-CSF" = "G-CSF", "IL-1RA" = "IL-1RA", "IL-1β" = "IL-1β",
         "IL-18" = "IL-18", "IFN-a2" = "IFN-α2", "CXCL8 (IL-8)" = "CXCL8", "IL-5" = "IL-5", "CCL4 (MIP-1β)" = "CCL4",
         "IL-2" = "IL-2", "TNF-a" = "TNF-α", "sCD25 (IL-2Ra)" = "sCD25", "IL-15" = "IL-15", "CXCL10 (IP-10)" = "CXCL10",
         "MIP-1α (CCL3)" = "CCL3", "IL-12p70" = "IL-12p70", "RANTES (CCL5)" = "RANTES", "IL-10" = "IL-10",
         "APRIL" = "APRIL")
padj <- setNames(st$p_adj_BH, st$analyte)
violin2 <- function(a, logy = FALSE) {
  df <- sv[sv$analyte == a, ]
  df$Condition <- factor(ifelse(df$side == "flare", "Flare", "Remission"), levels = c("Flare", "Remission"))
  df$pacienti <- df$pair_flare
  if (logy) df$value <- log10(df$value)
  pa <- padj[a]
  ttl <- if (is.na(pa) || a == "RANTES (CCL5)") lab[[a]] else sprintf("%s\n(adj. p = %s)", lab[[a]], formatC(pa, digits = 2, format = "fg", flag = "#"))
  bw <- diff(range(df$value)) / 25 + 1e-9
  ggplot(df[df$Condition == "Flare", ], aes(x = Condition, y = value, fill = Condition)) +
    geom_half_violin(side = "l") +
    geom_half_dotplot(stackdir = "up", binaxis = "y", dotsize = 0.6, binwidth = bw) +
    geom_half_violin(data = df[df$Condition == "Remission", ], aes(x = Condition, y = value, fill = Condition), side = "r") +
    geom_half_dotplot(data = df[df$Condition == "Remission", ], aes(x = Condition, y = value, fill = Condition),
                      stackdir = "down", binaxis = "y", dotsize = 0.6, binwidth = bw) +
    geom_line(data = df, mapping = aes(group = pacienti), alpha = 0.3) +
    scale_fill_manual(values = c("Flare" = "#E41A1C", "Remission" = "#377EB8")) +
    labs(title = ttl, x = "Condition", y = "Concentration (pg/mL)") +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7, angle = 90), axis.title.y = element_text(size = 8),
          panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey90", linetype = "dotted"))
}
grid2 <- function(analytes, logy = FALSE, ncol = 7) {
  pl <- lapply(analytes, violin2, logy = logy); n <- length(pl)
  for (i in seq_along(pl)) if (i > n - ncol) pl[[i]] <- pl[[i]] +
    theme(axis.text.x = element_text(size = 7, angle = 90), axis.title.x = element_text(size = 8)) + xlab("Condition")
  wrap_plots(pl, ncol = ncol) + plot_layout(axis_titles = "collect")
}
all2 <- intersect(names(lab), unique(sv$analyte))
save2(grid2(all2), file.path(OUT, "Figure2_cytokines"), 14, 10.5)
save2(grid2(all2, logy = TRUE), file.path(OPT, "Figure2_optB_log10"), 14, 10.5)
sig2 <- all2[all2 %in% st$analyte[!is.na(st$p_adj_BH) & st$p_adj_BH < 0.05] | all2 %in% c("IL-12p70", "IFN-a2", "TNF-a", "IL-2", "IL-15")]
save2(grid2(sig2, ncol = 6), file.path(OPT, "Figure2_optC_selected"), 12, 8)
cat("Figure 2: options A (all, linear), B (log10), C (12 raised + 5 key unchanged)\n")

## ======================================================================= Figure 3
nice_ct <- function(x) gsub("^gdT$", "γδ T", gsub("_", " ", x))
counts <- rd(R, "a6_composition", "nuclei_counts_per_subset_per_sample.tsv"); counts$cell_type <- nice_ct(counts$cell_type)
tests <- rd(R, "a6_composition", "composition_tests_recomputed.tsv"); tests <- tests[tests$set == "all", ]; tests$cell_type <- nice_ct(tests$cell_type)
long <- counts %>% pivot_longer(-cell_type, names_to = "s", values_to = "n") %>%
  separate(s, c("cond", "patient"), sep = "\\.") %>% group_by(cond, patient) %>% mutate(Proportion = n / sum(n)) %>% ungroup() %>%
  mutate(Group = factor(ifelse(cond == "flare", "MIS-C flare", "Remission"), levels = c("MIS-C flare", "Remission")),
         predicted.id = cell_type)
lev <- sort(unique(long$predicted.id)); long$predicted.id <- factor(long$predicted.id, levels = lev)
sig <- tests[tests$fdr_exact_30 < 0.05, c("cell_type", "fdr_exact_30")]
ymax <- long %>% group_by(predicted.id) %>% summarise(m = max(Proportion), .groups = "drop")
sig <- merge(sig, ymax, by.x = "cell_type", by.y = "predicted.id")
sig$x <- match(sig$cell_type, lev); sig$y.bracket <- sig$m + 0.012; sig$y.position <- sig$m + 0.03
sig$sig_label <- stars(sig$fdr_exact_30); sig$sig_color <- "#D55E00"
th_comp <- theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
        axis.text.y = element_text(size = 11, face = "bold", color = "black"),
        axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
        plot.background = element_rect(fill = "white", color = NA), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5), legend.position = "top",
        legend.title = element_blank(), legend.text = element_text(size = 11),
        legend.background = element_rect(fill = "gray95", color = NA))
pC <- ggplot(long, aes(x = predicted.id, y = Proportion, fill = Group)) +
  geom_rect(data = sig, aes(xmin = x - 0.45, xmax = x + 0.45, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#FFF3CD", alpha = 0.35) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.7, alpha = 0.8, color = "gray30") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.5, alpha = 0.6, color = "black") +
  geom_segment(data = sig, aes(x = x - 0.25, xend = x + 0.25, y = y.bracket, yend = y.bracket), inherit.aes = FALSE, linewidth = 0.6, color = "black") +
  geom_text(data = sig, aes(x = x, y = y.position, label = sig_label, color = sig_color), inherit.aes = FALSE, size = 8, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() + scale_fill_manual(values = c("MIS-C flare" = "#E41A1C", "Remission" = "#377EB8")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) + labs(x = "Cell type", y = "Proportion") + th_comp
grp <- c("CD4 Naive" = "CD4 T", "CD4 TCM" = "CD4 T", "CD4 TEM" = "CD4 T", "CD4 CTL" = "CD4 T", "CD4 Proliferating" = "CD4 T",
         "Treg" = "Treg", "CD8 Naive" = "CD8 T", "CD8 TCM" = "CD8 T", "CD8 TEM" = "CD8 T", "CD8 Proliferating" = "CD8 T",
         "MAIT" = "MAIT / gdT / dnT", "gdT" = "MAIT / gdT / dnT", "dnT" = "MAIT / gdT / dnT",
         "NK" = "NK", "NK_CD56bright" = "NK", "NK Proliferating" = "NK", "ILC" = "NK",
         "B naive" = "B", "B intermediate" = "B", "B memory" = "B", "Plasmablast" = "Plasmablast",
         "CD14 Mono" = "Monocytes", "CD16 Mono" = "Monocytes", "cDC1" = "DC", "cDC2" = "DC", "pDC" = "DC", "ASDC" = "DC")
lk <- read.delim(file.path(P, "wes", "patient_key.tsv"), stringsAsFactors = FALSE); lk <- lk[lk$manuscript_code != "-", ]
bar <- long %>% mutate(lineage = ifelse(as.character(predicted.id) %in% names(grp), grp[as.character(predicted.id)], "Other"),
                       Patient = factor(lk$manuscript_code[match(patient, lk$file_code)], levels = paste0("P", 1:10))) %>%
  group_by(Group, Patient, lineage) %>% summarise(Proportion = sum(Proportion), .groups = "drop")
pCbar <- ggplot(bar, aes(x = Patient, y = Proportion, fill = lineage)) + geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_wrap(~ Group, nrow = 1) + scale_fill_brewer(palette = "Set3", name = "Lineage") +
  labs(x = "Patient", y = "Proportion") + th_comp + theme(legend.position = "right", strip.text = element_text(face = "bold", size = 13))
cd4 <- c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM"); cd8 <- c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM")
rat <- long %>% group_by(Group, patient) %>% summarise(CD4_CD8_ratio = sum(n[cell_type %in% cd4]) / sum(n[cell_type %in% cd8]), .groups = "drop")
w <- pivot_wider(rat, names_from = Group, values_from = CD4_CD8_ratio); p_rat <- exact_paired(w$`MIS-C flare` - w$Remission)
pD <- ggplot(rat, aes(x = Group, y = CD4_CD8_ratio, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", width = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "black") +
  annotate("text", x = 1.5, y = max(rat$CD4_CD8_ratio) + 0.3, label = sprintf("p = %.3f", p_rat), size = 5) +
  scale_fill_manual(values = c("MIS-C flare" = "#E41A1C", "Remission" = "#377EB8")) +
  labs(x = NULL, y = "CD4⁺/CD8⁺ T cell ratio") + theme_minimal(base_size = 13) +
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size = 14, face = "bold"), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5), legend.position = "none")
if (file.exists(OBJ)) {
  suppressPackageStartupMessages(library(Seurat))
  obj <- readRDS(OBJ)
  obj$predicted.id <- nice_ct(as.character(obj$predicted.id))
  cond <- as.character(obj$Condition)
  obj$Condition <- factor(ifelse(cond %in% c("MISC", "MISC flare", "Flare"), "MIS-C flare", "Remission"), levels = c("Remission", "MIS-C flare"))
  thU <- theme_minimal(base_size = 16) + theme(legend.position = "right", legend.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"), plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"))
  pA <- DimPlot(obj, reduction = "umap_harmony", group.by = "predicted.id", label = TRUE, label.size = 3.5, label.box = FALSE,
                pt.size = 0.5, repel = TRUE, raster = TRUE) + thU + labs(title = NULL, color = "Cell type", x = "UMAP 1", y = "UMAP 2")
  pB <- DimPlot(obj, reduction = "umap_harmony", group.by = "Condition", cols = c("Remission" = "#377EB8", "MIS-C flare" = "#E41A1C"),
                pt.size = 0.5, raster = TRUE) + thU + theme(legend.text = element_text(size = 13, face = "bold")) +
    labs(title = NULL, color = "Disease state", x = "UMAP 1", y = "UMAP 2")
  rm(obj); gc()
  top <- pA + pB + plot_layout(widths = c(1.35, 1))
  tagth <- theme(plot.tag = element_text(size = 20, face = "bold"))
  save2((top / (pC + pD + plot_layout(widths = c(3, 1)))) + plot_layout(heights = c(1, 1.05)) +
          plot_annotation(tag_levels = "A") & tagth, file.path(OUT, "Figure3_composition"), 20, 16)
  save2((top / (pCbar + pD + plot_layout(widths = c(3, 1)))) + plot_layout(heights = c(1, 1.05)) +
          plot_annotation(tag_levels = "A") & tagth, file.path(OPT, "Figure3_optB_stackedbars"), 20, 16)
  cat("Figure 3: options A (box plots, as submitted), B (stacked bar plot per patient)\n")
}

## ======================================================================= Figure 4
hm <- msigdbr(species = "Homo sapiens", category = "H"); hm_sets <- split(hm$gene_symbol, hm$gs_name)
PRETTY <- c(HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF-α/NF-κB Signalling", HALLMARK_TGF_BETA_SIGNALING = "TGF-β Signalling",
            HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Oxidative Phosphorylation", HALLMARK_MYC_TARGETS_V1 = "MYC Targets v1",
            HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon-γ Response", HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon-α Response",
            HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory Response", HALLMARK_GLYCOLYSIS = "Glycolysis",
            HALLMARK_G2M_CHECKPOINT = "G2M Checkpoint", HALLMARK_E2F_TARGETS = "E2F Targets", HALLMARK_DNA_REPAIR = "DNA Repair")
keep_sets <- function(x) (grepl("inflamm|immune|cytokine|interferon|NF-kB|IL-|TNF-", x, ignore.case = TRUE)) |
  x %in% c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_GLYCOLYSIS", "HALLMARK_DNA_REPAIR")
run_fgsea <- function(stat) {
  stat <- sort(stat[is.finite(stat)], decreasing = TRUE)
  set.seed(20261025)
  r <- fgsea(pathways = hm_sets, stats = stat, minSize = 15, maxSize = 500, eps = 0)
  as.data.frame(r[, c("pathway", "pval", "padj", "NES", "size")])
}
bulk <- rd(R, "bulk_rerun", "bulk_deseq2_all_genes.tsv")
gsea_b <- do.call(rbind, lapply(c("main A_as_published", "validation B_pairs_fixed"), function(a) {
  b <- bulk[bulk$analysis == a & !is.na(bulk$lfcSE) & bulk$lfcSE > 0, ]
  r <- run_fgsea(setNames(b$log2FoldChange / b$lfcSE, b$gene)); r$cluster <- ifelse(grepl("^main", a), "Main", "Validation"); r }))
de_files <- list.files(file.path(R, "pseudobulk", "de"), full.names = TRUE)
gsea_s <- do.call(rbind, lapply(de_files, function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  r <- run_fgsea(setNames(sign(d$logFC) * sqrt(pmax(d$F, 0)), d$gene)); r$cluster <- d$cell_type[1]; r }))
write.table(rbind(transform(gsea_b, source = "Bulk"), transform(gsea_s, source = "snRNA-seq pseudobulk")),
            file.path(OUT, "fgsea_hallmark_bulk_pseudobulk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p4ab <- function(B, S, metric_lab) {
  a <- B %>% filter(keep_sets(pathway), padj < 0.05) %>% mutate(pathway_label = PRETTY[pathway], sig = padj < 0.05,
                                                                  Cohort = factor(cluster, levels = c("Main", "Validation")))
  ord <- a %>% group_by(pathway_label) %>% summarise(m = mean(NES)) %>% arrange(m) %>% pull(pathway_label)
  a$pathway_label <- factor(a$pathway_label, levels = ord)
  pa <- ggplot(a, aes(x = pathway_label, y = NES, fill = sig)) + geom_col(width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey55") + coord_flip(clip = "off") +
    facet_wrap(~ Cohort, nrow = 1, scales = "fixed") + scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "grey85"), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR ≥ 0.05"), name = NULL) +
    labs(tag = "A", x = "Pathway", y = metric_lab) + theme_minimal(base_size = 13.5) +
    theme(plot.tag = element_text(size = 18, face = "bold"), axis.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 11, face = "bold"),
          strip.text = element_text(size = 13, face = "bold", colour = "black"), strip.background = element_rect(fill = "#eaeef5", colour = NA),
          panel.grid.minor = element_blank(), panel.spacing.x = unit(1, "lines"))
  cts <- c("B memory", "B naive", "CD14 Mono", "CD16 Mono", "CD4 Naive", "CD4 TCM", "CD8 Naive", "CD8 TCM", "CD8 TEM", "gdT", "NK", "Plasmablast", "Treg")
  b <- rbind(transform(B, source = "Bulk"), transform(S[S$cluster %in% cts, ], source = "snRNA-seq pseudobulk")) %>%
    filter(keep_sets(pathway), pathway %in% names(PRETTY), padj < 0.05) %>%
    mutate(pathway_label = factor(PRETTY[pathway], levels = rev(PRETTY)), log_padj = pmin(-log10(padj), 15),
           enrichment_color = factor(ifelse(NES > 0, "red", "blue"), levels = c("red", "blue")),
           cluster = factor(cluster, levels = c("Main", "Validation", cts)), source = factor(source, levels = c("Bulk", "snRNA-seq pseudobulk")))
  pb <- ggplot(b, aes(x = cluster, y = pathway_label)) +
    geom_point(aes(size = log_padj, fill = enrichment_color), shape = 21, color = "black") +
    scale_size_continuous(range = c(3, 12), name = expression(-log[10] * "(adj. p)")) +
    scale_fill_manual(values = c("red", "blue"), labels = c("Higher at flare", "Lower at flare"), name = "Direction", drop = FALSE) +
    facet_grid(~ source, scales = "free_x", space = "free_x") +
    labs(tag = "B", x = "Cohort or population", y = "Pathway") + scale_x_discrete(labels = function(x) nice_ct(x)) + theme_minimal(base_size = 13) +
    theme(plot.tag = element_text(size = 18, face = "bold"), strip.text.x = element_text(face = "bold", size = 13, color = "black"),
          strip.background.x = element_rect(fill = "#eaeef5", color = "#eaeef5"), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 11, face = "bold"), axis.title = element_text(size = 13, face = "bold"))
  pa + pb + plot_layout(widths = c(1.05, 2.6))
}
prog <- list("Activation-associated TFs" = c("BATF", "PRDM1", "STAT1", "IRF4", "STAT3"),
             "IFN-γ core response (ISGs)" = c("IFITM1", "GBP1", "SOCS1", "PARP9"),
             "Classical pathway" = c("TBX21", "STAT4"), "IFNG / IFNG-AS1" = c("IFNG", "IFNG-AS1"),
             "Cytotoxic effectors" = c("GZMA", "GZMB", "PRF1"), "Checkpoint" = c("CD38", "CTLA4", "LAG3", "HAVCR2"),
             "Exhaustion-associated TFs" = c("NR4A3", "NR4A1", "NR4A2"), "Progenitor/repressor" = c("LEF1", "BACH2", "TCF7"),
             "TGF-β pathway" = c("TGFBR2", "SKI", "TGFB1", "SMAD7", "SKIL"),
             "Upstream gp130 / γc-receptor" = c("SOCS3", "IL6ST", "IL2RG", "IL15RA", "IL12RB1", "IL12RB2", "EBI3", "IL27"))
gt <- do.call(rbind, lapply(names(prog), function(k) data.frame(claim = k, gene = prog[[k]], y_pos = rev(seq_along(prog[[k]])))))
gt$claim <- factor(gt$claim, levels = names(prog))
bk <- rbind(transform(bulk[bulk$analysis == "main A_as_published", ], col = "Main cohort"),
            transform(bulk[bulk$analysis == "validation B_pairs_fixed", ], col = "Validation cohort"))
bk <- merge(gt, data.frame(gene = bk$gene, col = bk$col, lfc = bk$log2FoldChange, p = bk$padj, mode = "gw"), by = "gene", all.x = TRUE)
kg <- rd(R, "pseudobulk_followup", "key_gene_panel_pseudobulk.tsv")
sc_ct <- c("CD8 TEM" = "CD8 TEM", "NK" = "NK", "gdT" = "γδ T", "CD4 Naive" = "CD4 Naive", "CD8 Naive" = "CD8 Naive", "CD14 Mono" = "CD14 Mono")
sc <- do.call(rbind, lapply(names(sc_ct), function(ct) {
  m <- merge(gt, kg[kg$cell_type == ct, ], by = "gene", all.x = TRUE)
  gw <- suppressWarnings(as.numeric(m$FDR_genomewide)); pn <- suppressWarnings(as.numeric(m$FDR_panel))
  data.frame(m[, c("gene", "claim", "y_pos")], col = sc_ct[[ct]], lfc = as.numeric(m$logFC),
             p = ifelse(!is.na(gw) & gw < 0.05, gw, pn),
             mode = ifelse(!is.na(gw) & gw < 0.05, "gw", ifelse(!is.na(pn) & pn < 0.05, "panel", "ns")))
}))
dot <- function(df, levels_x, title, show_genes, tag = NULL) {
  df$x <- match(df$col, levels_x)
  df$sig <- ifelse(is.na(df$lfc), "na", ifelse(df$mode != "ns" & !is.na(df$p) & df$p < 0.05, "sig", "ns"))
  df$size <- pmin(-log10(pmax(df$p, 1e-50)), 50); df$fc <- pmax(pmin(df$lfc, 3), -3)
  g <- ggplot(df, aes(x = x, y = y_pos)) +
    geom_point(data = df[df$sig == "sig" & df$mode != "panel", ], aes(size = size, fill = fc), shape = 21, colour = "grey25", stroke = 0.4) +
    geom_point(data = df[df$sig == "sig" & df$mode == "panel", ], aes(size = size, fill = fc), shape = 23, colour = "grey25", stroke = 0.4) +
    geom_point(data = df[df$sig == "ns", ], shape = 1, size = 1.4, colour = "grey55", stroke = 0.6) +
    geom_point(data = df[df$sig == "na", ], shape = 22, size = 1.8, fill = "grey85", colour = "grey60", stroke = 0.4) +
    scale_fill_gradientn(colours = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
                         limits = c(-3, 3), breaks = c(-3, 0, 3), labels = c("≤ -3", "0", "≥ 3"), name = expression("log"[2] * "FC")) +
    scale_size_continuous(limits = c(0, 50), breaks = c(2, 5, 10, 25, 50), range = c(1.5, 6), name = expression(-log[10] * "(adj. p)")) +
    scale_x_continuous(breaks = seq_along(levels_x), labels = levels_x, limits = c(if (show_genes) -0.6 else 0.4, length(levels_x) + 0.6)) +
    scale_y_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    facet_grid(claim ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = title, tag = tag, x = NULL, y = NULL) + theme_minimal(base_size = 11) +
    theme(panel.grid.major = element_line(colour = "grey94", linewidth = 0.25), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), panel.spacing.y = unit(0.6, "lines"), panel.border = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9.5, face = "bold", colour = "#1F4E78", margin = margin(r = 8)),
          strip.placement = "outside", axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5, colour = "grey15"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11), plot.tag = element_text(size = 18, face = "bold"))
  if (show_genes) g + geom_text(aes(x = 0.35, label = gene), hjust = 1, size = 3.3, fontface = "italic") + coord_cartesian(clip = "off")
  else g + theme(strip.text.y.left = element_blank())
}
p4c <- dot(bk, c("Main cohort", "Validation cohort"), "Whole blood", TRUE, tag = "C") +
  dot(sc, unname(sc_ct), "Single-nucleus pseudobulk", FALSE) + plot_layout(widths = c(2.6, 5.9), guides = "collect")
fig4 <- function(top) {
  bottom <- wrap_elements(p4c)
  (wrap_elements(top) / bottom) + plot_layout(heights = c(0.85, 2.2))
}
save2(fig4(p4ab(gsea_b, gsea_s, "Normalised enrichment score")), file.path(OUT, "Figure4_programmes"), 15, 20, dpi = 120)
cat("Figure 4: options A (fgsea, as submitted), B (camera)\n")
print(gsea_b %>% filter(keep_sets(pathway)) %>% arrange(cluster, padj) %>% select(cluster, pathway, NES, padj))

## ======================================================================= Figure 5
kg5 <- kg[kg$gene %in% c("IFNG", "IFNG-AS1") & kg$cell_type %in% c("B naive", "Plasmablast", "B memory", "CD14 Mono", "CD4 Naive", "CD4 TCM",
                                                                  "CD8 Naive", "CD8 TCM", "CD8 TEM", "gdT", "NK", "NK_CD56bright"), ]
gw <- suppressWarnings(as.numeric(kg5$FDR_genomewide)); pn <- suppressWarnings(as.numeric(kg5$FDR_panel))
kg5$p_label <- ifelse(!is.na(gw) & gw < 0.05, stars(gw), ifelse(!is.na(pn) & pn < 0.05, paste0(stars(pn), "†"), "ns"))
kg5$logFC <- as.numeric(kg5$logFC); kg5$cell_type <- nice_ct(kg5$cell_type)
ord5 <- kg5 %>% filter(gene == "IFNG") %>% arrange(logFC) %>% pull(cell_type)
kg5$cell_type <- factor(kg5$cell_type, levels = unique(c(ord5, kg5$cell_type)))
p5a <- ggplot(kg5, aes(x = logFC, y = cell_type, fill = logFC > 0)) + geom_col(width = 0.7) +
  geom_text(aes(label = p_label, hjust = ifelse(logFC > 0, -0.3, 1.3)), size = 3) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") + facet_wrap(~ gene, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "#E15759", "FALSE" = "#4C72B0"), labels = c("Lower at flare", "Higher at flare"), guide = guide_legend(title = NULL)) +
  scale_x_continuous(labels = function(x) sprintf("%.1f", x), expand = expansion(add = 0.6)) +
  labs(x = expression(log[2] * " fold change"), y = NULL, tag = "A") + theme_bw(base_size = 14) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), strip.background = element_rect(fill = "#eaeef5"),
        axis.text.y = element_text(face = "bold"), panel.border = element_blank(), axis.line = element_line(colour = "grey40"),
        plot.tag = element_text(size = 18, face = "bold"))
pbd <- rd(R, "figure_data", "pb_ifng.tsv")
cts5 <- c("CD4 Naive", "CD4 TCM", "CD8 TCM", "CD8 TEM", "NK", "gdT", "B naive", "B memory", "CD14 Mono")
pbd <- pbd[pbd$cell_type %in% cts5 & pbd$n_cells >= 10, ]
pbd$Condition <- factor(ifelse(pbd$condition == "MISC", "Flare", "Remission"), levels = c("Remission", "Flare"))
pbd <- pbd %>% group_by(cell_type, gene, patient) %>% filter(n() == 2) %>% ungroup() %>% mutate(cpm = 2^log2cpm)
plots5 <- lapply(cts5, function(ct) {
  d <- pbd[pbd$cell_type == ct, ]; d$gene <- factor(d$gene, levels = c("IFNG", "IFNG-AS1"))
  ggplot(d, aes(x = Condition, y = cpm)) + geom_violin(aes(fill = Condition), width = 0.6, alpha = 0.3, color = NA) +
    geom_line(aes(group = patient), color = "grey30", linewidth = 0.3, alpha = 0.6) + geom_point(aes(color = Condition), size = 1.2) +
    facet_wrap(~ gene, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(Remission = "#4C72B0", Flare = "#E15759"), guide = "none") +
    scale_color_manual(values = c(Remission = "#4C72B0", Flare = "#E15759"), guide = "none") +
    labs(title = nice_ct(ct), x = NULL, y = NULL) + theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10), strip.text = element_text(face = "bold", size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8), panel.grid = element_blank(),
          panel.border = element_blank(), axis.line = element_line(colour = "grey40"), plot.tag = element_text(size = 18, face = "bold"))
})
plots5[[1]] <- plots5[[1]] + labs(tag = "B")
p5b <- wrap_plots(plots5, nrow = 1)
bk5 <- rd(R, "figure_data", "bulk_ifng.tsv"); bk5$Expression <- 2^bk5$log2norm - 1
bk5$Condition <- factor(ifelse(bk5$condition == "flare", "Flare", "Remission"), levels = c("Remission", "Flare"))
bk5$Cohort <- factor(ifelse(bk5$cohort == "main", "Main", "Validation"), levels = c("Main", "Validation"))
bk5$Gene <- factor(bk5$gene, levels = c("IFNG", "IFNG-AS1"))
pl <- bk5 %>% select(pair, Gene, Cohort, Condition, Expression) %>% pivot_wider(names_from = Condition, values_from = Expression) %>%
  filter(!is.na(Remission), !is.na(Flare)) %>% mutate(x_rem = 1 - 0.14, x_flr = 2 + 0.14)
p5c <- ggplot() +
  geom_half_violin(data = bk5[bk5$Condition == "Remission", ], aes(x = Condition, y = Expression, fill = Condition), side = "l", trim = FALSE, width = 0.4, position = position_nudge(x = -0.1), alpha = 0.5) +
  geom_half_violin(data = bk5[bk5$Condition == "Flare", ], aes(x = Condition, y = Expression, fill = Condition), side = "r", trim = FALSE, width = 0.4, position = position_nudge(x = 0.1), alpha = 0.5) +
  geom_segment(data = pl, aes(x = x_rem, xend = x_flr, y = Remission, yend = Flare), color = "grey60", linewidth = 0.3, alpha = 0.7) +
  geom_point(data = bk5, aes(x = Condition, y = Expression, color = Condition), size = 1.5, position = position_nudge(x = ifelse(bk5$Condition == "Remission", -0.1, 0.1))) +
  facet_grid(Gene ~ Cohort, scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(Remission = "#4C72B0", Flare = "#E15759"), guide = "none") +
  scale_color_manual(values = c(Remission = "#4C72B0", Flare = "#E15759"), guide = "none") +
  scale_y_continuous(labels = scales::label_comma(accuracy = 1)) +
  labs(x = "Condition", y = "Normalised counts", tag = "C") +
  theme_bw(base_size = 12) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0), plot.tag = element_text(size = 18, face = "bold"),
        panel.border = element_blank(), axis.line = element_line(colour = "grey40"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        strip.background.x = element_rect(fill = "#eaeef5", colour = NA), strip.text.x = element_text(face = "bold", size = 11, colour = "black"))
f5 <- function(extra = NULL) {
  parts <- list(wrap_elements(p5a), wrap_elements(p5b), wrap_elements(p5c)); h <- c(1.5, 1.1, 1.3)
  if (!is.null(extra)) { parts <- c(parts, list(wrap_elements(extra))); h <- c(h, 0.9) }
  wrap_plots(parts, ncol = 1, heights = h)
}
save2(f5(), file.path(OUT, "Figure5_IFNG_IFNGAS1"), 12, 15, dpi = 130)
