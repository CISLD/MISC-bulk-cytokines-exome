# Whole-blood cell composition for the composition-adjusted DESeq2 model (bulk_rerun.R, model C).
#
# Method: the LM22 signature matrix (Newman et al., Nat Methods 2015) fitted to each sample's gene TPM by
# non-negative least squares (Lawson-Hanson active-set algorithm, written out below), fractions rescaled
# to sum to 1. This is not CIBERSORT (no nu-SVR, no permutation p values). LM22.txt is not part of this
# repository: download it from the CIBERSORTx website (https://cibersortx.stanford.edu) under its terms.
#
# Neutrophil fraction = "Neutrophils". Lymphocyte fraction = CD4 T + CD8 T (incl. gamma delta) + NK + B
# (incl. plasma cells), grouped as below. bulk_rerun.R uses the columns sample_code, neut and lymph.
#
# Usage: Rscript lm22_nnls_composition.R <salmon.merged.gene_tpm.tsv> <LM22.txt> <out.tsv>
#   The TPM table is the nf-core/rnaseq Salmon output (columns gene_id, gene_name, then one per sample).
#   sample_code = the leading "<letter><digits>-<digits>" part of each sample column name.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
tpm_file <- args[1]; lm22_file <- args[2]; out_file <- args[3]
set.seed(20260917L)   # only used by the synthetic-mixture check

tpm <- read.delim(tpm_file, check.names = FALSE)
tc <- colnames(tpm)[-(1:2)]
codes <- regmatches(tc, regexpr("^[A-Za-z][0-9]+-[0-9]+", tc))
stopifnot(length(codes) == length(tc), !anyDuplicated(codes))
tm <- as.matrix(tpm[, -(1:2)]); colnames(tm) <- codes
tm <- rowsum(tm, tpm$gene_name)[, codes]

lm22 <- as.matrix(read.delim(lm22_file, check.names = FALSE, row.names = 1))
genes <- intersect(rownames(lm22), rownames(tm))
cat("LM22 genes found:", length(genes), "of", nrow(lm22), "\n")
stopifnot(length(genes) >= 480)
S <- lm22[genes, ]

nnls_lh <- function(A, b, tol = 1e-10, maxit = 500) {
  n <- ncol(A); x <- numeric(n); P <- logical(n)
  w <- drop(crossprod(A, b - A %*% x)); it <- 0
  while (any(!P & w > tol * max(1, abs(w))) && it < maxit) {
    it <- it + 1
    cand <- which(!P)
    P[cand[which.max(w[cand])]] <- TRUE
    repeat {
      z <- numeric(n)
      z[P] <- qr.coef(qr(A[, P, drop = FALSE]), b)
      z[is.na(z)] <- 0
      if (all(z[P] > 0)) { x <- z; break }
      neg <- P & z <= 0
      a <- min(x[neg] / (x[neg] - z[neg]))
      x <- x + a * (z - x)
      P <- P & x > tol
      x[!P] <- 0
    }
    w <- drop(crossprod(A, b - A %*% x))
  }
  x
}
frac_of <- function(y) { x <- nnls_lh(S, y); x / sum(x) }

# implementation check: noise-free synthetic mixtures must be recovered
syn_err <- sapply(1:50, function(k) {
  f <- rgamma(ncol(S), 1); f <- f / sum(f)
  max(abs(frac_of(drop(S %*% f)) - f))
})
cat(sprintf("synthetic mixtures: max abs error %.2e\n", max(syn_err)))
stopifnot(max(syn_err) < 0.01)

fr <- t(sapply(codes, function(s) frac_of(tm[genes, s])))
colnames(fr) <- colnames(S)
fit_r <- sapply(codes, function(s) cor(tm[genes, s], drop(S %*% nnls_lh(S, tm[genes, s]))))
cat(sprintf("fit correlation: min %.3f\n", min(fit_r)))
stopifnot(all(fit_r >= 0.5), all(is.finite(fr)), all(abs(rowSums(fr) - 1) < 1e-8))

grp <- list(
  neut = "Neutrophils",
  mono = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
           "Dendritic cells resting", "Dendritic cells activated"),
  cd4 = c("T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
          "T cells follicular helper", "T cells regulatory (Tregs)"),
  cd8 = c("T cells CD8", "T cells gamma delta"),
  nk = c("NK cells resting", "NK cells activated"),
  bcell = c("B cells naive", "B cells memory", "Plasma cells"),
  othergran = c("Mast cells resting", "Mast cells activated", "Eosinophils"))
stopifnot(setequal(unlist(grp), colnames(fr)))
gf <- sapply(grp, function(g) rowSums(fr[, g, drop = FALSE]))
gf <- cbind(gf, lymph = gf[, "cd4"] + gf[, "cd8"] + gf[, "nk"] + gf[, "bcell"])

comp <- data.frame(sample_code = codes, fit_r = fit_r, gf, fr, check.names = FALSE)
write.table(comp, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("written:", out_file, "(", nrow(comp), "samples )\n")
