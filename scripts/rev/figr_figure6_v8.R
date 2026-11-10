#!/usr/bin/env Rscript
# v8 (2026-11-08): Figure 6 keeps panels B-C (heatmaps); the TBX21-target and three-level panels move to
# Supplementary Figure S4 (A, B). Writes <out_prefix>_BC.* and <out_prefix>_S4.* (A4 portrait page).
# v5: panels relabelled B-E (panel A of Figure 6 is the STAT-motif chromVAR panel, make_figures_v5.py)
# Revised Figure 6: motif-supported TF-DORC associations in flare and remission.
#   A, B  Spearman rho (smoothed TF RNA vs DORC score) for the most connected TFs in CD8 TEM and
#         gdT, flare and remission side by side, DORCs grouped by Table S8a panel.
#   C     TBX21 targets in CD8 TEM: rho in remission -> rho in flare, patient-level permutation p.
#   D     Motif-supported Table S8 edges with a flare-remission difference (FDR < 0.05) at three
#         units of analysis.
#
# usage: Rscript figr_figure6.R <CD8_TEM_differential_edges.tsv> <gdT_differential_edges.tsv>
#                               <Table_S8a_long.tsv> <three_level_SUMMARY.tsv> <out_prefix>
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args <- commandArgs(trailingOnly = TRUE)
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
for (l in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(Sys.setlocale("LC_CTYPE", l))) break
edges <- rbind(rd(args[1]), rd(args[2]))
panels <- read.delim(args[3], skip = 1, stringsAsFactors = FALSE)
summ <- rd(args[4])
out <- args[5]

MIN_RHO <- 0.2    # edge shown if |rho| >= MIN_RHO in either condition
N_TF <- 12        # TFs per heatmap

edges <- subset(edges, Enrichment.P < 0.05)
edges$ct <- ifelse(grepl("CD8", edges$cell_type), "CD8 TEM", "gdT")
panels$ct <- ifelse(grepl("CD8", panels$cell_type), "CD8 TEM", "gdT")
panels <- panels[!duplicated(panels[, c("ct", "gene")]), ]
edges$panel <- panels$gene_set[match(paste(edges$ct, edges$DORC), paste(panels$ct, panels$gene))]
edges$panel[edges$panel == "IFNG"] <- "IFN-γ"
edges$panel <- factor(edges$panel, c("IFN-γ", "TGF-β", "Checkpoint"))
edges$strong <- pmax(abs(edges$rho_flare), abs(edges$rho_remission)) >= MIN_RHO

heat <- function(cell, title) {
  e <- subset(edges, ct == cell & strong & !is.na(panel))
  top <- names(sort(table(e$Motif), decreasing = TRUE))[seq_len(min(N_TF, length(unique(e$Motif))))]
  e <- subset(e, Motif %in% top)
  long <- rbind(data.frame(e[, c("Motif", "DORC", "panel")], condition = "Flare", rho = e$rho_flare),
                data.frame(e[, c("Motif", "DORC", "panel")], condition = "Remission", rho = e$rho_remission))
  long$condition <- factor(long$condition, c("Flare", "Remission"))
  long$Motif <- factor(long$Motif, rev(top))
  ggplot(long, aes(DORC, Motif, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_grid(condition ~ panel, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-0.8, 0.8),
                         oob = scales::squish, name = "Spearman ρ") +
    labs(tag = title, x = "DORC (target gene)", y = "Transcription factor") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = element_blank(),
          strip.text = element_text(face = "bold"), plot.tag = element_text(face = "bold", size = 13))
}

pA <- heat("CD8 TEM", "B")
pB <- heat("gdT", "C")

tb <- subset(edges, ct == "CD8 TEM" & Motif == "TBX21" & strong)
tb <- tb[order(tb$rho_flare), ]
tb$DORC <- factor(tb$DORC, tb$DORC)
tb$lab <- sprintf("p = %.2f", tb$p_perm)
pC <- ggplot(tb, aes(y = DORC)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(x = rho_remission, xend = rho_flare, yend = DORC), colour = "grey55") +
  geom_point(aes(x = rho_remission, colour = "Remission"), size = 2) +
  geom_point(aes(x = rho_flare, colour = "Flare"), size = 2) +
  geom_text(aes(x = 1.02, label = lab), hjust = 0, size = 2.4, colour = "grey30") +
  scale_colour_manual(values = c(Flare = "#b2182b", Remission = "#2166ac"), name = NULL) +
  scale_x_continuous(limits = c(-0.8, 1.25), breaks = seq(-0.75, 0.75, 0.25)) +
  labs(tag = "A", x = "Spearman ρ (TBX21 RNA vs DORC score)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 13), panel.grid.minor = element_blank())

lvl <- c("Nucleus", "Nucleus, k-NN corrected", "Patient")
d <- rbind(data.frame(ct = summ$cell_type, unit = lvl[1], n = summ$motif_sig_cell, tested = summ$motif_supported_tested),
           data.frame(ct = summ$cell_type, unit = lvl[2], n = summ$motif_sig_knn, tested = summ$motif_supported_tested),
           data.frame(ct = summ$cell_type, unit = lvl[3], n = summ$motif_sig_patient, tested = summ$patient_tested))
d$ct <- ifelse(grepl("CD8", d$ct), "CD8 TEM", "γδ T")
d$unit <- factor(d$unit, lvl)
pD <- ggplot(d, aes(unit, n, fill = unit)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d / %d", n, tested)), vjust = -0.3, size = 2.4) +
  facet_wrap(~ ct) +
  scale_fill_manual(values = c("#bdbdbd", "#737373", "#252525"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(tag = "B", x = "Unit of analysis", y = "Edges, FDR < 0.05") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), plot.tag = element_text(face = "bold", size = 13),
        panel.grid.major.x = element_blank())

fig <- pA / pB
ggsave(paste0(out, "_BC.pdf"), fig, width = 9, height = 9.5, device = cairo_pdf)
ggsave(paste0(out, "_BC.png"), fig, width = 9, height = 9.5, dpi = 400)
s4 <- (pC / pD / plot_spacer()) + plot_layout(heights = c(1.3, 1, 0.9))
s4 <- wrap_elements(s4) + theme(plot.margin = margin(50, 60, 40, 60))
ggsave(paste0(out, "_S4.pdf"), s4, width = 8.27, height = 11.69, device = cairo_pdf)
ggsave(paste0(out, "_S4.png"), s4, width = 8.27, height = 11.69, dpi = 300)

cat("TBX21 targets:", nrow(tb), "\n")
