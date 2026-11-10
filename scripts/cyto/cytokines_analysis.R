#!/usr/bin/env Rscript
# Plasma cytokines (LEGENDplex COVID-19 Cytokine Storm Panels 1 and 2), flare vs remission.
#
# Input: wells_long.tsv from cytokines_parse.py (one row per well).
# Primary analysis: plasma at 2x dilution (the only matrix measured for every sample), complete
# flare-remission pairs, exact paired Wilcoxon signed-rank test, Benjamini-Hochberg over analytes.
#
# usage: Rscript cytokines_analysis.R <wells_long.tsv> <out_dir> [n_boot] [drop_pairs, e.g. S6]

args   <- commandArgs(trailingOnly = TRUE)
wells  <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE)
outdir <- args[2]
n_boot <- if (length(args) >= 3) as.integer(args[3]) else 10000L
drop   <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else character(0)  # flare IDs of pairs to drop
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261021)

w <- subset(wells, matrix == "DIL2X" & !(pair_flare %in% drop))
# Panel 1, unpaired flare samples: well 2 reads ~5x well 1 for every analyte (a dilution-factor
# error in one well); keep well 1, as the original summary sheet did. Paired samples are unaffected.
w <- subset(w, !(panel == "panel1_14plex" & side == "flare" & pair_flare == "" & well == 2))
w$value <- suppressWarnings(as.numeric(w$value))
w$limit <- suppressWarnings(as.numeric(w$limit))

# per-analyte lower limit: the largest '<x' printed at this dilution (conservative when plates differ)
lloq <- tapply(w$limit[w$flag == "below"], w$analyte[w$flag == "below"], max)
bp   <- w$flag == "below" & w$pair_flare != ""
lloq_paired <- tapply(w$limit[bp], w$analyte[bp], max)   # reported: limit on the paired samples' plate
uloq <- tapply(w$limit[w$flag == "above"], w$analyte[w$flag == "above"], min)

# ---- sample-level values under a given substitution for wells below the LLOQ -------------------
sample_values <- function(sub_fun) {
  x <- w
  lim <- ifelse(is.na(x$limit), lloq[x$analyte], x$limit)
  x$v <- ifelse(x$flag == "ok", x$value,
         ifelse(x$flag == "below", sub_fun(lim),
         ifelse(x$flag == "above", x$limit, NA)))
  agg <- aggregate(cbind(v, below = flag == "below", above = flag == "above") ~
                     analyte + side + sample + pair_flare, data = x,
                   FUN = function(z) c(mean = mean(z), n = length(z)))
  data.frame(analyte = agg$analyte, side = agg$side, sample = agg$sample,
             pair_flare = agg$pair_flare, value = agg$v[, "mean"],
             n_wells = agg$v[, "n"],
             wells_below = agg$below[, "mean"] * agg$below[, "n"],
             wells_above = agg$above[, "mean"] * agg$above[, "n"],
             stringsAsFactors = FALSE)
}

# ---- exact paired signed-rank test on log2 differences (ties: mid-ranks; zeros dropped) --------
exact_signed_rank <- function(d) {
  d <- d[abs(d) > 1e-12]
  n <- length(d)
  if (n == 0) return(list(p = 1, n = 0, V = NA))
  r <- rank(abs(d))
  V <- sum(r[d > 0])
  signs <- as.matrix(expand.grid(rep(list(c(0, 1)), n)))
  Vperm <- signs %*% r
  mu <- sum(r) / 2
  p <- mean(abs(Vperm - mu) >= abs(V - mu) - 1e-9)
  list(p = min(1, p), n = n, V = V)
}

sign_test <- function(k, n) if (n == 0) NA else binom.test(k, n, 0.5)$p.value

boot_ci <- function(d, B = n_boot) {
  if (length(d) < 3) return(c(NA, NA))
  m <- replicate(B, median(sample(d, replace = TRUE)))
  unname(quantile(m, c(0.025, 0.975)))
}

analytes <- sort(unique(w$analyte))
subs <- list(`LOD/sqrt2` = function(l) l / sqrt(2), `LOD/2` = function(l) l / 2, LOD = function(l) l)

all_res <- list()
for (s in names(subs)) {
  sv <- sample_values(subs[[s]])
  if (s == "LOD/sqrt2") write.table(sv, file.path(outdir, "sample_values_DIL2X.tsv"),
                                    sep = "\t", quote = FALSE, row.names = FALSE)
  rows <- lapply(analytes, function(a) {
    x  <- subset(sv, analyte == a)
    fl <- subset(x, side == "flare"); rm <- subset(x, side == "remission")
    pr <- merge(fl, rm, by = "pair_flare", suffixes = c(".f", ".r"))
    pr <- pr[pr$pair_flare != "", ]
    d  <- log2(pr$value.f / pr$value.r)
    sat <- sum(x$wells_above > 0)
    t  <- exact_signed_rank(d)
    # sign test on pairs whose order is decidable despite censoring
    fb <- pr$wells_below.f == pr$n_wells.f; rb <- pr$wells_below.r == pr$n_wells.r
    up <- (!fb & rb) | (!fb & !rb & pr$value.f > pr$value.r)
    dn <- (fb & !rb) | (!fb & !rb & pr$value.f < pr$value.r)
    ci <- boot_ci(d)
    data.frame(
      analyte = a, substitution = s,
      n_flare = nrow(fl), n_remission = nrow(rm), n_pairs = nrow(pr),
      pct_below_flare = round(100 * mean(fl$wells_below == fl$n_wells), 1),
      pct_below_remission = round(100 * mean(rm$wells_below == rm$n_wells), 1),
      samples_above_uloq = sat,
      lloq_2x = unname(lloq_paired[a]), uloq_2x = unname(uloq[a]),
      median_flare_all = median(fl$value), q1_flare_all = quantile(fl$value, .25, names = FALSE),
      q3_flare_all = quantile(fl$value, .75, names = FALSE),
      median_flare_paired = median(pr$value.f),
      q1_flare_paired = quantile(pr$value.f, .25, names = FALSE),
      q3_flare_paired = quantile(pr$value.f, .75, names = FALSE),
      median_remission = median(rm$value),
      q1_remission = quantile(rm$value, .25, names = FALSE),
      q3_remission = quantile(rm$value, .75, names = FALSE),
      pairs_up = sum(d > 0), pairs_down = sum(d < 0),
      median_log2FC = median(d), log2FC_ci_low = ci[1], log2FC_ci_high = ci[2],
      p_wilcoxon_exact = t$p,
      sign_up = sum(up), sign_down = sum(dn), p_sign = sign_test(sum(up), sum(up) + sum(dn)),
      p_mannwhitney_all_flare = suppressWarnings(wilcox.test(fl$value, rm$value, exact = FALSE)$p.value),
      stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  # an analyte is not quantifiable when > 20% of samples have a well above the upper limit
  res$testable <- res$samples_above_uloq <= 0.2 * (res$n_flare + res$n_remission)
  res$heavily_censored <- pmax(res$pct_below_flare, res$pct_below_remission) >= 50
  res$p_adj_BH <- NA
  res$p_adj_BH[res$testable] <- p.adjust(res$p_wilcoxon_exact[res$testable], "BH")
  res$p_sign_adj_BH <- NA
  res$p_sign_adj_BH[res$testable] <- p.adjust(res$p_sign[res$testable], "BH")
  all_res[[s]] <- res
}

primary <- all_res[["LOD/sqrt2"]]
primary <- primary[order(primary$p_wilcoxon_exact), ]
write.table(primary, file.path(outdir, "cytokines_primary_LODsqrt2.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
sens <- do.call(rbind, all_res)
write.table(sens[, c("analyte", "substitution", "median_log2FC", "p_wilcoxon_exact", "p_adj_BH",
                     "pairs_up", "pairs_down")],
            file.path(outdir, "cytokines_sensitivity_substitution.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# ---- duplicate-well agreement (both wells detected) --------------------------------------------
dw <- subset(w, flag == "ok")
cvt <- aggregate(value ~ analyte + side + sample, data = dw,
                 FUN = function(z) if (length(z) == 2) 100 * sd(z) / mean(z) else NA)
names(cvt)[4] <- "cv_pct"
cvt <- cvt[!is.na(cvt$cv_pct), ]
cv_sample <- aggregate(cv_pct ~ side + sample, data = cvt,
                       FUN = function(z) c(median = median(z), n = length(z), over30 = sum(z > 30)))
cv_sample <- data.frame(side = cv_sample$side, sample = cv_sample$sample,
                        median_cv_pct = round(cv_sample$cv_pct[, "median"], 1),
                        analytes_with_both_wells = cv_sample$cv_pct[, "n"],
                        analytes_cv_over_30pct = cv_sample$cv_pct[, "over30"])
write.table(cv_sample, file.path(outdir, "duplicate_cv_by_sample.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# ---- minimum attainable p and power note -------------------------------------------------------
n_pairs <- max(primary$n_pairs)
cat(sprintf("pairs: %d; smallest exact two-sided p: %.4f\n", n_pairs, 2 / 2^n_pairs))
cat(sprintf("testable analytes: %d of %d\n", sum(primary$testable), nrow(primary)))
print(primary[, c("analyte", "n_pairs", "pct_below_flare", "pct_below_remission", "samples_above_uloq",
                  "median_flare_paired", "median_remission", "pairs_up", "pairs_down",
                  "median_log2FC", "log2FC_ci_low", "log2FC_ci_high", "p_wilcoxon_exact",
                  "p_adj_BH", "p_sign")], row.names = FALSE, digits = 3)
print(cv_sample, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
