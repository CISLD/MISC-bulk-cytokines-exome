# v11: Figure 5A labels without the panel FDR. Genes below the primary expression filter take the genome-wide FDR of the
#      alternative full-design filter (pseudobulk_designfilter_panel.R in the companion repository, marked ‡) or are not tested (nt).
# Figures 2 and 5 of the article.
#   Fig 2  eight flare-remission pairs, plasma 1:2, values below LLOQ set to LLOQ/sqrt(2); exact paired Wilcoxon, BH over 26 analytes
#   Fig 5  paired pseudobulk IFNG / IFNG-AS1; per-patient pseudobulk CPM; whole-blood normalised counts
# Figure 3 is drawn by figure3_v8.R and Figure 4 by figure4_camera_v11.R.
# usage: Rscript figures_article_v11.R <project dir> <out dir> <R lib>
for (l in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(Sys.setlocale("LC_CTYPE", l))) break
args <- commandArgs(trailingOnly = TRUE)
P <- args[1]; OUT <- args[2]; .libPaths(c(args[3], .libPaths()))
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

nice_ct <- function(x) gsub("^gdT$", "γδ T", gsub("_", " ", x))
kg <- rd(R, "pseudobulk_followup", "key_gene_panel_pseudobulk.tsv")

## ======================================================================= Figure 5
kg5 <- kg[kg$gene %in% c("IFNG", "IFNG-AS1") & kg$cell_type %in% c("B naive", "Plasmablast", "B memory", "CD14 Mono", "CD4 Naive", "CD4 TCM",
                                                                  "CD8 Naive", "CD8 TCM", "CD8 TEM", "gdT", "NK", "NK_CD56bright"), ]
# v11: no panel FDR. Below the primary expression filter: genome-wide FDR of the alternative full-design filter (‡) or not tested (nt)
ds <- rd(R, "pseudobulk_designfilter_panel", "key_genes_designfilter_panel.tsv"); ds <- ds[ds$passed_design_filter %in% c("TRUE", TRUE), ]
gw <- suppressWarnings(as.numeric(kg5$FDR_genomewide))
al <- suppressWarnings(as.numeric(ds$FDR_genomewide[match(paste(kg5$cell_type, kg5$gene), paste(ds$cell_type, ds$gene))]))
kg5$p_label <- ifelse(!is.na(gw), ifelse(gw < 0.05, stars(gw), "ns"), ifelse(is.na(al), "nt", ifelse(al < 0.05, paste0(stars(al), "‡"), "ns‡")))
print(kg5[, c("cell_type", "gene", "logFC", "FDR_genomewide", "p_label")])
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
