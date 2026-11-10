# Bulk whole-blood RNA-seq (PAXgene, Salmon), flare vs remission: re-run of the
# published DESeq2 analysis with two corrections.
#
#  A  published model, reproduced from the saved pipeline environment
#     (design ~ Patient_ID + Condition, DESeq2 Wald, BH)
#  B  validation cohort with the split pair re-joined (one patient's two samples carry
#     differently spelled Patient_IDs, so the published design treats them as two
#     unpaired individuals)
#  C  composition-adjusted model: ~ Patient_ID + neutrophil + lymphocyte + Condition,
#     using LM22 NNLS fractions estimated from the same samples
#
# Patient identifiers are never written out: pairs are recoded as pair_01, pair_02, ...
#
# Usage: bulk_rerun.R MAIN_ENV VALIDATION_ENV COMPOSITION_TSV S7_TSV OUT_DIR

suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
MAIN <- args[1]; VAL <- args[2]; COMP <- args[3]; S7 <- args[4]; OUT <- args[5]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

panel <- c("BATF", "STAT1", "STAT3", "PRDM1", "IRF4", "STAT4", "TBX21", "IFNG", "IFNG-AS1", "IL6ST",
           "IL12RB1", "IL12RB2", "GBP1", "SOCS1", "SOCS3", "CTLA4", "HAVCR2", "GZMB", "TCF7", "CXCL10")

load_env <- function(f) {
  e <- new.env(); load(f, envir = e)
  cd <- e$coldata; cnt <- as.matrix(e$countdata)
  stopifnot(identical(colnames(cnt), as.character(cd$Sample_ID)) || all(colnames(cnt) %in% cd$Sample_ID))
  cd <- cd[match(colnames(cnt), cd$Sample_ID), ]
  sym <- sub("^gene-", "", rownames(cnt))
  storage.mode(cnt) <- "integer"
  rownames(cnt) <- make.unique(sym)
  list(cd = cd, cnt = cnt, formula = e$formula, min_count = e$min_count)
}
recode_pairs <- function(pid) {
  u <- unique(pid); setNames(sprintf("pair_%02d", seq_along(u)), u)[pid]
}
run_deseq <- function(cnt, cd, design, label) {
  dds <- DESeqDataSetFromMatrix(cnt, cd, design)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds, quiet = TRUE)
  rn <- resultsNames(dds); coef <- rn[grepl("^Condition", rn)]
  r <- as.data.frame(results(dds, name = coef))
  r$gene <- rownames(r); r$analysis <- label
  r
}
cond_factor <- function(x) factor(ifelse(tolower(x) %in% c("misc", "flare"), "flare", "remission"),
                                  levels = c("remission", "flare"))

comp <- read.delim(COMP, check.names = FALSE, stringsAsFactors = FALSE)
res <- list(); info <- list()
for (cohort in c("main", "validation")) {
  E <- load_env(if (cohort == "main") MAIN else VAL)
  cd <- E$cd
  cd$Condition <- cond_factor(cd$Condition)
  tb <- table(cd$Patient_ID)
  info[[cohort]] <- data.frame(cohort = cohort, samples = nrow(cd), patient_ids = length(tb),
                               complete_pairs = sum(tb == 2), singletons = sum(tb == 1),
                               saved_formula = paste(deparse(E$formula), collapse = ""))
  # A: as published
  cdA <- data.frame(row.names = cd$Sample_ID, Patient_ID = factor(recode_pairs(cd$Patient_ID)), Condition = cd$Condition)
  res[[paste(cohort, "A")]] <- run_deseq(E$cnt, cdA, ~ Patient_ID + Condition, paste(cohort, "A_as_published"))
  # B: re-join split pair (exactly one flare-only and one remission-only Patient_ID)
  pid <- as.character(cd$Patient_ID)
  single <- names(tb)[tb == 1]
  if (length(single) == 2) {
    conds <- sapply(single, function(s) as.character(cd$Condition[pid == s]))
    stopifnot(setequal(conds, c("flare", "remission")))
    pid[pid == single[2]] <- single[1]
  }
  cdB <- data.frame(row.names = cd$Sample_ID, Patient_ID = factor(recode_pairs(pid)), Condition = cd$Condition)
  if (length(single) == 2) {
    res[[paste(cohort, "B")]] <- run_deseq(E$cnt, cdB, ~ Patient_ID + Condition, paste(cohort, "B_pairs_fixed"))
  }
  # C: composition-adjusted (on the fixed pairing)
  # Sample_ID starts with the sample code used in the composition table (e.g. P<n>-<yy>)
  m <- match(sub("^([A-Za-z][0-9]+-[0-9]+).*$", "\\1", cd$Sample_ID), comp$sample_code)
  if (all(!is.na(m))) {
    cdC <- cdB
    cdC$neut <- as.numeric(scale(comp$neut[m])); cdC$lymph <- as.numeric(scale(comp$lymph[m]))
    res[[paste(cohort, "C")]] <- run_deseq(E$cnt, cdC, ~ Patient_ID + neut + lymph + Condition,
                                           paste(cohort, "C_composition_adjusted"))
    info[[cohort]]$composition_matched <- sum(!is.na(m))
    d <- tapply(comp$neut[m], cd$Condition, mean); l <- tapply(comp$lymph[m], cd$Condition, mean)
    info[[cohort]]$neut_flare_minus_rem <- round(d["flare"] - d["remission"], 3)
    info[[cohort]]$lymph_flare_minus_rem <- round(l["flare"] - l["remission"], 3)
  } else {
    info[[cohort]]$composition_matched <- sum(!is.na(m))
  }
}
info <- do.call(rbind, info)
all <- do.call(rbind, res)
write.table(info, file.path(OUT, "cohort_design_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all[, c("analysis", "gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")],
            file.path(OUT, "bulk_deseq2_all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# reproduction check against published Table S7 (main and validation log2FC / adj.p)
s7 <- read.delim(S7, skip = 3, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
names(s7)[1:6] <- c("Category", "Gene", "lfc_main", "padj_main", "lfc_val", "padj_val")
s7 <- s7[!is.na(suppressWarnings(as.numeric(s7$lfc_main))) & s7$Gene != "", ]
s7$lfc_main <- as.numeric(s7$lfc_main); s7$lfc_val <- as.numeric(s7$lfc_val)
s7 <- s7[!duplicated(s7$Gene), ]
A_main <- all[all$analysis == "main A_as_published", ]; A_val <- all[all$analysis == "validation A_as_published", ]
rep <- data.frame(gene = s7$Gene,
                  published_main = s7$lfc_main, rerun_main = A_main$log2FoldChange[match(s7$Gene, A_main$gene)],
                  published_val = s7$lfc_val, rerun_val = A_val$log2FoldChange[match(s7$Gene, A_val$gene)])
write.table(rep, file.path(OUT, "reproduction_vs_TableS7.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

key <- all[all$gene %in% panel, c("analysis", "gene", "log2FoldChange", "padj")]
w <- reshape(key, idvar = "gene", timevar = "analysis", direction = "wide")
write.table(w, file.path(OUT, "key_genes_by_analysis.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

options(width = 250)
cat("== Design check\n"); print(info, row.names = FALSE)
cat("\n== Reproduction of Table S7 log2FC (as-published model)\n")
cat("main: r =", round(cor(rep$published_main, rep$rerun_main, use = "complete.obs"), 4),
    "| max |diff| =", round(max(abs(rep$published_main - rep$rerun_main), na.rm = TRUE), 3), "| n =", sum(complete.cases(rep[, 2:3])), "\n")
cat("validation: r =", round(cor(rep$published_val, rep$rerun_val, use = "complete.obs"), 4),
    "| max |diff| =", round(max(abs(rep$published_val - rep$rerun_val), na.rm = TRUE), 3), "| n =", sum(complete.cases(rep[, 4:5])), "\n")
cat("\n== Genome-wide DE counts (padj < 0.05)\n")
print(tapply(all$padj < 0.05, all$analysis, sum, na.rm = TRUE))
cat("\n== Key genes: log2FC (padj) by analysis\n")
for (g in panel) {
  x <- key[key$gene == g, ]
  if (!nrow(x)) next
  cat(sprintf("%-9s", g), paste(sprintf("%s %+.2f (%.1e)", sub(" .*", "", sub("_.*", "", sub("^(main|validation) ", "\\1:", x$analysis))),
                                        x$log2FoldChange, x$padj), collapse = " | "), "\n")
}
