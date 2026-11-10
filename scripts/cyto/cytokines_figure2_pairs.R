#!/usr/bin/env Rscript
# Revised Figure 2, 8 complete pairs only: flare-remission plasma cytokines (2x dilution), one line per patient,
# dashed line at the lower limit of quantification, BH-adjusted exact paired Wilcoxon p-values.
#
# usage: Rscript cytokines_figure2.R <sample_values_DIL2X.tsv> <cytokines_primary_LODsqrt2.tsv> <out_prefix>
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
sv  <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE)
res <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE)
out <- args[3]

show <- c("IFN-y", "CXCL10 (IP-10)", "IL-18", "IL-6", "sCD25 (IL-2Ra)", "IL-10", "VEGF", "IL-12p70")
label <- c("IFN-y" = "IFN-gamma", "CXCL10 (IP-10)" = "CXCL10", "IL-18" = "IL-18", "IL-6" = "IL-6",
           "sCD25 (IL-2Ra)" = "sCD25 (IL-2Ralpha)", "IL-10" = "IL-10", "VEGF" = "VEGF",
           "IL-12p70" = "IL-12p70")

sv <- subset(sv, analyte %in% show)
sv$pair <- ifelse(sv$side == "flare", sv$sample, sv$pair_flare)
sv$paired <- sv$pair_flare != ""
sv$below <- sv$wells_below == sv$n_wells
sv$condition <- factor(ifelse(sv$side == "flare", "Flare", "Remission"), c("Flare", "Remission"))

res <- subset(res, analyte %in% show)
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))
res$facet <- sprintf("%s\nFC %.2f, adj. p = %s", label[res$analyte], 2^res$median_log2FC, fmt_p(res$p_adj_BH))
sv$facet  <- factor(res$facet[match(sv$analyte, res$analyte)], levels = res$facet[match(show, res$analyte)])
lod <- data.frame(facet = factor(res$facet, levels = levels(sv$facet)), lloq = res$lloq_2x)
lod <- lod[!is.na(lod$lloq), ]

pp <- subset(sv, paired)
g <- ggplot(pp, aes(condition, value)) +
  geom_hline(data = lod, aes(yintercept = lloq), linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_line(aes(group = pair), colour = "grey60", linewidth = 0.4) +
  geom_point(aes(fill = condition, shape = below), size = 2.2, stroke = 0.5) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), labels = c("Quantified", "Below LLOQ (plotted at LLOQ/sqrt2)"),
                     name = NULL) +
  scale_fill_manual(values = c(Flare = "#C0392B", Remission = "#2E86C1"), guide = "none") +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
  facet_wrap(~ facet, nrow = 2, scales = "free_y") +
  labs(x = NULL, y = "Concentration (pg/mL, log scale)",
       caption = paste("Lines: 8 paired patients.",
                       "Dashed line: lower limit of quantification. FC: median paired fold change;",
                       "adj. p: exact paired Wilcoxon, Benjamini-Hochberg over 26 analytes.", sep = "\n")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"),
        panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0, size = 7))
ggsave(paste0(out, ".pdf"), g, width = 7.2, height = 5.3, device = cairo_pdf)
ggsave(paste0(out, ".png"), g, width = 7.2, height = 5.3, dpi = 300)
cat("written", out, "\n")
