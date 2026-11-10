# Per-sample whole-blood expression of IFNG and IFNG-AS1 (Figure 5C): DESeq2 size-factor-normalised
# log2(count + 1), both cohorts, validation split pair re-joined as in bulk_rerun.R. Pairs are recoded.
# usage: Rscript figure_data_bulk.R MAIN_ENV VALIDATION_ENV OUT_TSV
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
out <- list()
for (k in 1:2) {
  e <- new.env(); load(args[k], envir = e)
  cd <- e$coldata; cnt <- as.matrix(e$countdata)
  cd <- cd[match(colnames(cnt), cd$Sample_ID), ]
  rownames(cnt) <- make.unique(sub("^gene-", "", rownames(cnt)))
  storage.mode(cnt) <- "integer"
  pid <- as.character(cd$Patient_ID)
  tb <- table(pid); single <- names(tb)[tb == 1]
  if (length(single) == 2) pid[pid == single[2]] <- single[1]
  u <- unique(pid); pair <- setNames(sprintf("pair_%02d", seq_along(u)), u)[pid]
  cond <- ifelse(tolower(cd$Condition) %in% c("misc", "flare"), "flare", "remission")
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt, data.frame(row.names = colnames(cnt), c = factor(cond)), ~ c))
  nc <- counts(dds, normalized = TRUE)
  for (g in c("IFNG", "IFNG-AS1"))
    out[[length(out) + 1]] <- data.frame(cohort = c("main", "validation")[k], pair = pair, condition = cond,
                                         gene = g, log2norm = log2(nc[g, ] + 1))
}
write.table(do.call(rbind, out), args[3], sep = "\t", quote = FALSE, row.names = FALSE)
