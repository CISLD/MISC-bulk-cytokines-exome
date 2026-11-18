# v11: panel C without the panel FDR (the targeted panel only selects the genes shown). Genes below the primary
#      expression filter take the genome-wide FDR of the alternative full-design filter (triangle) or are not tested.
# v9: whole-blood genes filtered on the flare/remission groups (filterByExpr(group = condition)), as in the
#     primary pseudobulk analysis; key genes from the corrected pseudobulk_followup.R (primary fit).
# Figure 4 of the article, gene sets tested with camera throughout.
#   A  whole blood, main and validation cohorts: camera on paired edgeR models (~ patient + condition), Hallmark sets;
#      bar = signed -log10 FDR. In the validation cohort the two time points of one child carry different patient codes
#      in the count file and are re-joined as one pair (as in scripts/bulk/bulk_rerun.R).
#   B  the same whole-blood camera results with pseudobulk camera (results/pseudobulk/camera_hallmark_all.tsv from
#      pseudobulk_edgeR.R in the companion repository, paired edgeR QL models)
#   C  key genes: whole-blood DESeq2 (results/bulk_rerun) and paired pseudobulk (results/pseudobulk_followup), genes below the
#      primary expression filter from results/pseudobulk_designfilter_panel (triangles)
# Also writes camera_hallmark_bulk_pseudobulk.tsv (source data of Figure 4A-B).
# usage: Rscript figure4_camera_v11.R <project dir> <out dir> <R lib> <main cohort .RData> <validation cohort .RData>
#   each .RData holds coldata (Sample_ID, Patient_ID, Condition) and countdata (salmon gene counts)
for (l in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(Sys.setlocale("LC_CTYPE", l))) break
args <- commandArgs(trailingOnly = TRUE)
P <- args[1]; OUT <- args[2]; .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(edgeR); library(limma); library(msigdbr)
})
OPT <- file.path(OUT, "options")
R <- file.path(P, "results")
rd <- function(...) read.delim(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)
save2 <- function(p, file, w, h, dpi = 150) {
  if (grepl("/options/", file)) return(invisible(NULL))
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 300)
}
nice_ct <- function(x) gsub("^gdT$", "γδ T", gsub("_", " ", x))
hm <- msigdbr(species = "Homo sapiens", category = "H"); hm_sets <- split(hm$gene_symbol, hm$gs_name)
PRETTY <- c(HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF-α/NF-κB Signalling", HALLMARK_TGF_BETA_SIGNALING = "TGF-β Signalling",
            HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Oxidative Phosphorylation", HALLMARK_MYC_TARGETS_V1 = "MYC Targets v1",
            HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon-γ Response", HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon-α Response",
            HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory Response", HALLMARK_GLYCOLYSIS = "Glycolysis",
            HALLMARK_G2M_CHECKPOINT = "G2M Checkpoint", HALLMARK_E2F_TARGETS = "E2F Targets", HALLMARK_DNA_REPAIR = "DNA Repair")
keep_sets <- function(x) (grepl("inflamm|immune|cytokine|interferon|NF-kB|IL-|TNF-", x, ignore.case = TRUE)) |
  x %in% c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_GLYCOLYSIS", "HALLMARK_DNA_REPAIR")
# camera, whole blood: paired edgeR models on the salmon count matrices (same inputs as bulk_rerun.R)
envs <- c(Main = args[4], Validation = args[5])
cam_b <- do.call(rbind, lapply(names(envs), function(co) {
  e <- new.env(); load(envs[[co]], envir = e)
  cd <- e$coldata; cnt <- as.matrix(e$countdata); cd <- cd[match(colnames(cnt), cd$Sample_ID), ]
  rownames(cnt) <- make.unique(sub("^gene-", "", rownames(cnt)))
  pid <- as.character(cd$Patient_ID); tb <- table(pid); single <- names(tb)[tb == 1]
  if (length(single) == 2) pid[pid == single[2]] <- single[1]
  cond <- factor(ifelse(tolower(cd$Condition) %in% c("misc", "flare"), "Flare", "Remission"), levels = c("Remission", "Flare"))
  ks <- pid %in% names(which(table(pid) == 2))
  cat(co, ": pairs =", sum(ks) / 2, "\n")
  y <- DGEList(round(cnt[, ks])); des <- model.matrix(~ factor(pid[ks]) + cond[ks])
  y <- calcNormFactors(y[filterByExpr(y, group = cond[ks]), , keep.lib.sizes = FALSE]); y <- estimateDisp(y, des)
  cam <- camera(y, ids2indices(hm_sets, rownames(y)), des, contrast = ncol(des))
  data.frame(pathway = rownames(cam), NGenes = cam$NGenes, Direction = cam$Direction, pval = cam$PValue, padj = cam$FDR,
             NES = ifelse(cam$Direction == "Up", 1, -1) * -log10(cam$FDR), size = cam$NGenes, cluster = co)
}))
bulk <- rd(R, "bulk_rerun", "bulk_deseq2_all_genes.tsv")  # gene-level whole-blood results for panel C
cs <- rd(R, "pseudobulk", "camera_hallmark_all.tsv")
cam_s <- data.frame(pathway = cs$set, NGenes = cs$NGenes, Direction = cs$Direction, pval = cs$PValue, padj = cs$FDR,
                    NES = ifelse(cs$Direction == "Up", 1, -1) * -log10(cs$FDR), size = cs$NGenes, cluster = cs$cell_type)
out_tab <- rbind(transform(cam_b, source = "Whole blood"), transform(cam_s, source = "snRNA-seq pseudobulk"))
names(out_tab)[names(out_tab) == "NES"] <- "signed_log10_FDR"
write.table(out_tab[, c("source", "cluster", "pathway", "NGenes", "Direction", "pval", "padj", "signed_log10_FDR")],
            file.path(OUT, "camera_hallmark_bulk_pseudobulk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
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
    scale_fill_manual(values = c(red = "red", blue = "blue"), labels = c(red = "Higher at flare", blue = "Lower at flare"), name = "Direction") +
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
# v11: no panel FDR. Genes below the primary expression filter take the genome-wide FDR of the alternative full-design filter
# analysis (mode "alt", triangle) if they pass that filter, otherwise they are drawn as not tested.
ds <- rd(R, "pseudobulk_designfilter_panel", "key_genes_designfilter_panel.tsv")
ds <- ds[ds$passed_design_filter %in% c("TRUE", TRUE), ]
kg$alt <- suppressWarnings(as.numeric(ds$FDR_genomewide[match(paste(kg$cell_type, kg$gene), paste(ds$cell_type, ds$gene))]))
sc_ct <- c("CD8 TEM" = "CD8 TEM", "NK" = "NK", "gdT" = "γδ T", "CD4 Naive" = "CD4 Naive", "CD8 Naive" = "CD8 Naive", "CD14 Mono" = "CD14 Mono")
sc <- do.call(rbind, lapply(names(sc_ct), function(ct) {
  m <- merge(gt, kg[kg$cell_type == ct, ], by = "gene", all.x = TRUE)
  gw <- suppressWarnings(as.numeric(m$FDR_genomewide)); al <- m$alt
  data.frame(m[, c("gene", "claim", "y_pos")], col = sc_ct[[ct]],
             lfc = ifelse(is.na(gw) & is.na(al), NA, as.numeric(m$logFC)),   # below both filters: not tested
             p = ifelse(!is.na(gw), gw, al),
             mode = ifelse(!is.na(gw), ifelse(gw < 0.05, "gw", "ns"), ifelse(!is.na(al) & al < 0.05, "alt", "ns")))
}))
dot <- function(df, levels_x, title, show_genes, tag = NULL) {
  df$x <- match(df$col, levels_x)
  df$sig <- ifelse(is.na(df$lfc), "na", ifelse(df$mode != "ns" & !is.na(df$p) & df$p < 0.05, "sig", "ns"))
  df$size <- pmin(-log10(pmax(df$p, 1e-50)), 50); df$fc <- pmax(pmin(df$lfc, 3), -3)
  g <- ggplot(df, aes(x = x, y = y_pos)) +
    geom_point(data = df[df$sig == "sig" & df$mode != "alt", ], aes(size = size, fill = fc), shape = 21, colour = "grey25", stroke = 0.4) +
    geom_point(data = df[df$sig == "sig" & df$mode == "alt", ], aes(size = size, fill = fc), shape = 24, colour = "grey25", stroke = 0.4) +
    geom_point(data = df[df$sig == "ns", ], shape = 1, size = 1.4, colour = "grey55", stroke = 0.6) +
    geom_point(data = df[df$sig == "na", ], shape = 22, size = 1.8, fill = "grey85", colour = "grey60", stroke = 0.4) +
    scale_fill_gradientn(colours = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
                         limits = c(-3, 3), breaks = c(-3, 0, 3), labels = c("≤ -3", "0", "≥ 3"), name = expression("log"[2] * "FC")) +
    scale_size_continuous(limits = c(0, 50), breaks = c(2, 5, 10, 25, 50), range = c(1.5, 6), name = expression(-log[10] * "(adj. p)"),
                          guide = guide_legend(override.aes = list(shape = 21, fill = "grey70", colour = "grey25"))) +
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
save2(fig4(p4ab(cam_b, cam_s, expression("Signed " * -log[10] * " FDR (camera)"))), file.path(OUT, "Figure4_programmes"), 15, 20, dpi = 120)
print(cam_b %>% filter(keep_sets(pathway), padj < 0.05) %>% arrange(cluster, padj) %>% select(cluster, pathway, Direction, padj))
