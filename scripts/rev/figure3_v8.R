# Figure 3 of the article: (A, B) Harmony UMAP coloured by population and by disease state,
# (C) population proportions per sample with exact paired Wilcoxon tests (BH over 30 populations), (D) CD4/CD8 ratio.
# usage: Rscript figure3_v8.R <project dir> <harmony-UMAP seurat .rds> <out dir> <R lib>
for (l in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(Sys.setlocale("LC_CTYPE", l))) break
args <- commandArgs(trailingOnly = TRUE)
P <- args[1]; OBJ <- args[2]; OUT <- args[3]; .libPaths(c(args[4], .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(gghalves); library(dplyr); library(tidyr)
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
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, seed = 20261110), size = 1.5, alpha = 0.6, color = "black") +
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
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed = 20261110), size = 2, alpha = 0.8, color = "black") +
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
  suppressPackageStartupMessages(library(ggrepel))
  emb <- as.data.frame(Embeddings(obj, "umap_harmony")); colnames(emb) <- c("u1", "u2")
  emb$ct <- obj$predicted.id; emb$cond <- obj$Condition
  cen <- emb %>% group_by(ct) %>% summarise(u1 = median(u1), u2 = median(u2), .groups = "drop")
  RDPI <- c(2400, 2400)
  pA <- DimPlot(obj, reduction = "umap_harmony", group.by = "predicted.id", label = FALSE,
                pt.size = 5, raster = TRUE, raster.dpi = RDPI) + thU +
    geom_label_repel(data = cen, aes(x = u1, y = u2, label = ct), inherit.aes = FALSE, size = 4.2, fontface = "bold",
                     fill = alpha("white", 0.75), label.size = NA, label.padding = unit(0.12, "lines"),
                     box.padding = 0.25, min.segment.length = 0.2, segment.colour = "grey30", max.overlaps = Inf, seed = 1) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 4))) +
    labs(title = NULL, color = "Cell type", x = "UMAP 1", y = "UMAP 2")
  set.seed(20261110)
  pB <- DimPlot(obj, reduction = "umap_harmony", group.by = "Condition", cols = c("Remission" = "#377EB8", "MIS-C flare" = "#E41A1C"),
                pt.size = 5, raster = TRUE, raster.dpi = RDPI, shuffle = TRUE, seed = 20261110) + thU +
    theme(legend.text = element_text(size = 13, face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    labs(title = NULL, color = "Disease state", x = "UMAP 1", y = "UMAP 2")
  rm(obj); gc()
  top <- pA + pB + plot_layout(widths = c(1.35, 1))
  tagth <- theme(plot.tag = element_text(size = 20, face = "bold"))
  save2((top / (pC + pD + plot_layout(widths = c(3, 1)))) + plot_layout(heights = c(1, 1.05)) +
          plot_annotation(tag_levels = "A") & tagth, file.path(OUT, "Figure3_composition"), 20, 16)
  cat("Figure 3: options A (box plots, as submitted), B (stacked bar plot per patient)\n")
}

