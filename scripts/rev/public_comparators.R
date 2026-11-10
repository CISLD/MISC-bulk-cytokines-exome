# Specificity of the flare transcriptional pattern: public paediatric comparator cohorts.
#
# Question (reviewers R1.3, R1.5, R2.9, R3.3): is the BATF/STAT1/PRDM1(/IRF4) induction,
# STAT4 reduction and IFNG/IFNG-AS1 pattern specific to MIS-C, or shared with other acute
# paediatric inflammation?
#
# Datasets (GEO, downloaded processed matrices):
#   GSE225220  whole-blood RNA-seq counts: MIS-C acute, paediatric COVID-19 acute,
#              MIS-C >= 3 months, non-inflammatory controls
#   GSE73461   whole-blood Illumina HT-12: Kawasaki disease, definite viral, definite
#              bacterial, other inflammatory (JIA + HSP), healthy controls
#   GSE66099   whole-blood Affymetrix U133 Plus 2: septic shock, sepsis, SIRS, controls
#   GSE167028  PBMC RNA-seq (VST): MIS-C, MIS-C with myocarditis, Kawasaki disease, controls
#
# Analysis (identical for all datasets):
#   - one value per gene (highest-mean probe for arrays)
#   - limma, each disease group vs the dataset's own controls
#   - model A: ~ group; model B: ~ group + composition scores (mean z of marker genes for
#     neutrophils, T cells, NK cells, monocytes; whole blood only)
#   - "effector TF module" = mean z (vs controls) of BATF, STAT1, PRDM1, IRF4

suppressPackageStartupMessages({ library(limma); library(edgeR) })
args <- commandArgs(trailingOnly = TRUE)
DATA <- args[1]; META <- args[2]; OUT <- args[3]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

genes <- c("BATF", "STAT1", "PRDM1", "IRF4", "STAT3", "IL6ST", "STAT4", "TBX21", "IFNG", "IFNG-AS1",
           "GBP1", "SOCS1", "SOCS3", "IL12RB1", "IL12RB2", "CTLA4", "HAVCR2", "TCF7", "CXCL10")
module <- c("BATF", "STAT1", "PRDM1", "IRF4")
markers <- list(neutrophil = c("CSF3R", "FCGR3B", "CXCR2", "FPR1"),
                T_cell = c("CD3E", "CD3D", "CD2", "LCK"),
                NK = c("NKG7", "GNLY", "KLRD1"),
                monocyte = c("CD14", "LYZ", "CSF1R"))

collapse_probes <- function(mat, sym) {
  ok <- !is.na(sym) & sym != "" & !grepl("///", sym)
  mat <- mat[ok, , drop = FALSE]; sym <- sym[ok]
  m <- rowMeans(mat)
  o <- order(sym, -m); keep <- o[!duplicated(sym[o])]
  out <- mat[keep, , drop = FALSE]; rownames(out) <- sym[keep]; out
}
read_annot <- function(f) {
  x <- readLines(gzfile(f)); start <- grep("^ID\t", x)[1]
  a <- read.delim(text = x[start:length(x)], quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  a <- a[!grepl("^!", a$ID), ]
  setNames(a[["Gene symbol"]], a$ID)
}
series_groups <- function(f, field_regex) {
  x <- readLines(gzfile(f))
  ids <- strsplit(x[grep("^!Sample_geo_accession", x)], "\t")[[1]][-1]
  ch <- x[grep("^!Sample_characteristics_ch1|^!Sample_source_name_ch1|^!Sample_title", x)]
  for (l in ch) {
    v <- gsub('"', "", strsplit(l, "\t")[[1]][-1])
    if (any(grepl(field_regex, v))) return(list(ids = gsub('"', "", ids), values = sub(paste0("^", field_regex), "", v)))
  }
  stop("field not found")
}
analyse <- function(name, E, group, control, whole_blood) {
  group <- factor(group); group <- relevel(group, ref = control)
  sc <- sapply(markers, function(g) { g <- intersect(g, rownames(E)); colMeans(t(scale(t(E[g, , drop = FALSE])))) })
  res <- list(); modres <- list()
  for (adj in c(FALSE, if (whole_blood) TRUE)) {
    design <- if (adj) model.matrix(~ group + sc) else model.matrix(~ group)
    fit <- eBayes(lmFit(E, design))
    for (lv in levels(group)[-1]) {
      cf <- paste0("group", lv)
      tt <- topTable(fit, coef = cf, number = Inf, sort.by = "none")
      g <- intersect(genes, rownames(tt))
      res[[paste(lv, adj)]] <- data.frame(dataset = name, contrast = paste(lv, "vs", control),
                                          composition_adjusted = adj, gene = g, logFC = tt[g, "logFC"],
                                          P = tt[g, "P.Value"], FDR_panel = p.adjust(tt[g, "P.Value"], "BH"))
    }
    # module score
    ctl <- group == control
    z <- (E[intersect(module, rownames(E)), , drop = FALSE] - rowMeans(E[intersect(module, rownames(E)), ctl, drop = FALSE])) /
         apply(E[intersect(module, rownames(E)), ctl, drop = FALSE], 1, sd)
    ms <- colMeans(z)
    mfit <- eBayes(lmFit(matrix(ms, nrow = 1), design))
    for (lv in levels(group)[-1]) {
      cf <- paste0("group", lv)
      modres[[paste(lv, adj)]] <- data.frame(dataset = name, contrast = paste(lv, "vs", control),
                                             composition_adjusted = adj, module_mean_z_cases = mean(ms[group == lv]),
                                             delta = mfit$coefficients[1, cf], P = mfit$p.value[1, cf],
                                             n_cases = sum(group == lv), n_controls = sum(ctl))
    }
  }
  list(genes = do.call(rbind, res), module = do.call(rbind, modres))
}

all_g <- list(); all_m <- list()

## GSE225220: whole-blood RNA-seq counts
cnt <- read.delim(file.path(DATA, "GSE225220_wbrna_counts.tsv.gz"), check.names = FALSE, row.names = 1)
cnt <- t(as.matrix(cnt))
sym <- sub("^[^_]+_", "", rownames(cnt))
meta <- read.delim(file.path(DATA, "GSE225220_wbrna_supplementary.tsv.gz"), stringsAsFactors = FALSE)
meta <- meta[match(colnames(cnt), meta$wbrna_sample_id), ]
grp <- make.names(paste(meta$Diagnosis, meta$timepoint))
y <- DGEList(cnt); y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
E <- cpm(y, log = TRUE, prior.count = 2)
E <- collapse_probes(E, sub("^[^_]+_", "", rownames(E)))
cat("GSE225220 groups:\n"); print(table(grp))
keep <- grp %in% c("MIS.C.acute", "COVID.19.acute", "MIS.C.3m.", "Control_Non.inflammatory.Not.hospitalized")
r <- analyse("GSE225220 whole blood", E[, keep], grp[keep], "Control_Non.inflammatory.Not.hospitalized", TRUE)
all_g[[1]] <- r$genes; all_m[[1]] <- r$module

## GSE73461: whole-blood Illumina HT-12 (background-subtracted, not log)
x <- read.delim(file.path(DATA, "GSE73461_GEOupload_Discovery_Dataset_Normalised_Sept_15_n_459.txt.gz"),
                check.names = FALSE, row.names = 1)
x <- x[, !grepl("Detection_Pval", colnames(x))]
x <- as.matrix(x[, colSums(!is.na(x)) > 0])   # file ends with an empty column
x[x < 1] <- 1; x <- log2(x + 15)
# rows are Illumina array-address IDs; map through the full platform table
pl <- readLines(file.path(DATA, "GPL10558_table.txt"))
pl <- pl[grep("^ID\t", pl)[1]:length(pl)]; pl <- pl[!grepl("^!", pl)]
pl <- read.delim(text = pl, quote = "", stringsAsFactors = FALSE, check.names = FALSE)
ann <- setNames(pl$Symbol, as.character(pl$Array_Address_Id))
E <- collapse_probes(x, ann[rownames(x)])
grp <- sub("_[0-9]+$", "", colnames(E))
cat("GSE73461 groups:\n"); print(table(grp))
# "Inflammatory" = JIA and HSP in the series metadata (66 + 18)
keep <- grp %in% c("Control", "Kawasaki", "DefiniteViral", "DefiniteBacterial", "Inflammatory")
r <- analyse("GSE73461 whole blood", E[, keep], grp[keep], "Control", TRUE)
all_g[[2]] <- r$genes; all_m[[2]] <- r$module

## GSE66099: whole-blood U133 Plus 2 (series matrix values)
f <- file.path(META, "GSE66099_series_matrix.txt.gz")
L <- readLines(gzfile(f)); s <- grep("^!series_matrix_table_begin", L); e <- grep("^!series_matrix_table_end", L)
x <- read.delim(text = L[(s + 1):(e - 1)], row.names = 1, check.names = FALSE)
x <- as.matrix(x); if (max(x, na.rm = TRUE) > 100) x <- log2(pmax(x, 1))
ann <- read_annot(file.path(DATA, "GPL570.annot.gz"))
E <- collapse_probes(x, ann[rownames(x)])
gi <- series_groups(f, "disease: ")
grp <- gi$values[match(colnames(E), gi$ids)]
cat("GSE66099 groups:\n"); print(table(grp))
r <- analyse("GSE66099 whole blood", E, grp, "Control", TRUE)
all_g[[3]] <- r$genes; all_m[[3]] <- r$module

## GSE167028: PBMC RNA-seq (VST)
x <- read.csv2(file.path(DATA, "GSE167028_vst_normalised_matrix.csv.gz"), row.names = 1, check.names = FALSE)
x <- as.matrix(x); storage.mode(x) <- "numeric"
grp <- sub("^P[0-9]+_", "", colnames(x)); grp <- sub("\\(.*$", "", grp)
grp[grepl("^C[0-9]+$", grp)] <- "CTL"          # controls are named C20, C21, ...
cat("GSE167028 groups:\n"); print(table(grp))
keep <- grp %in% c("MIS-C", "MIS-C_MYO", "KD", "CTL")
grp2 <- make.names(grp)
r <- analyse("GSE167028 PBMC", x[, keep], grp2[keep], "CTL", FALSE)
all_g[[4]] <- r$genes; all_m[[4]] <- r$module

G <- do.call(rbind, all_g); M <- do.call(rbind, all_m)
write.table(G, file.path(OUT, "comparator_gene_panel.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(M, file.path(OUT, "comparator_module_score.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 220)
cat("\n== Effector TF module (BATF/STAT1/PRDM1/IRF4), each group vs controls\n")
print(M, row.names = FALSE, digits = 3)
cat("\n== Key genes, composition-adjusted where available\n")
key <- G[G$gene %in% c("BATF", "STAT1", "PRDM1", "IRF4", "STAT4", "IFNG", "IFNG-AS1", "IL6ST") &
         (G$composition_adjusted | G$dataset == "GSE167028 PBMC"), ]
w <- reshape(key[, c("contrast", "dataset", "gene", "logFC")], idvar = c("dataset", "contrast"),
             timevar = "gene", direction = "wide")
print(w, row.names = FALSE, digits = 2)
