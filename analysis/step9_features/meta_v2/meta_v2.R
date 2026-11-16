#!/usr/bin/env Rscript
# meta_v2.R — extension of ../f1_meta.R (same inputs, same DL estimator).
# Adds per gene: tau2, 95% prediction interval (t-distribution, df = k-1),
# leave-one-cohort-out (LOO) meta logFC min/max and a sign-flip flag.
# The original columns are recomputed with the verbatim f1 dl_meta() code and
# ASSERTED equal to the frozen ../f1_meta.csv before anything is written.
suppressPackageStartupMessages(library(data.table))
setwd(".")

primary <- list(   # verbatim from f1_meta.R
  GSE136825     = "GSE136825_polyp_vs_ctrl",
  GSE338426     = "GSE338426_CRSwNP_vs_ctrl",
  GSE179265     = "GSE179265_CRSwNP_vs_ctrl",
  GSE36830      = "GSE36830_CRSwNP_vs_UTctrl",
  GSE194282     = "GSE194282_NP_vs_ctrl",
  GSE69093      = "GSE69093_case_vs_ctrl_celltypeAdj",
  GSE313186     = "GSE313186_exposed_vs_ctrl",
  GSE70898_70900= "GSE70898_70900_polypT_vs_nasalT",
  GSE198950     = "GSE198950_CRSwNP_vs_NonCRS",
  GSE72713      = "GSE72713_CRSwNP_vs_ctrl"
)
all_list <- lapply(names(primary), function(gse) {
  f <- file.path("analysis/step2_bulk/de", paste0(primary[[gse]], ".csv"))
  d <- fread(f, select = c("symbol", "logFC", "t"))
  d <- d[!is.na(t) & t != 0 & !is.na(logFC) & symbol != ""]
  d[, cohort := gse]
  d[, se := abs(logFC / t)]
  d[se < 1e-8, se := 1e-8]
  unique(d[, .(gene = symbol, cohort, logFC, se)], by = c("gene", "cohort"))
})
comb <- rbindlist(all_list)
cat("combined rows:", nrow(comb), " unique genes:", length(unique(comb$gene)), "\n")

## DL fit returning the full stat set (f1 formulas verbatim + tau2/PI/LOO)
dl_full <- function(y, v) {
  k <- length(y)
  w <- 1 / v
  mu_fe <- sum(w * y) / sum(w)
  Q <- sum(w * (y - mu_fe)^2)
  df <- k - 1
  tau2 <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
  w2 <- 1 / (v + tau2)
  mu <- sum(w2 * y) / sum(w2)
  se_mu <- sqrt(1 / sum(w2))
  z <- mu / se_mu
  p <- 2 * pnorm(-abs(z))
  ci_lo <- mu - 1.96 * se_mu; ci_hi <- mu + 1.96 * se_mu
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  # 95% prediction interval, t-distribution with df = k-1
  tcrit <- qt(0.975, df = k - 1)
  pi_lo <- mu - tcrit * sqrt(tau2 + se_mu^2)
  pi_hi <- mu + tcrit * sqrt(tau2 + se_mu^2)
  # leave-one-cohort-out
  loo <- sapply(seq_len(k), function(j) {
    yj <- y[-j]; vj <- v[-j]
    wj <- 1 / vj
    mfe <- sum(wj * yj) / sum(wj)
    Qj <- sum(wj * (yj - mfe)^2)
    dfj <- length(yj) - 1
    t2j <- max(0, (Qj - dfj) / (sum(wj) - sum(wj^2) / sum(wj)))
    w2j <- 1 / (vj + t2j)
    sum(w2j * yj) / sum(w2j)
  })
  flip <- any((mu > 0 & loo < 0) | (mu < 0 & loo > 0))
  c(k = k, meta_logFC = mu, ci_lo = ci_lo, ci_hi = ci_hi, meta_p = p, I2 = I2,
    tau2 = tau2, pi_lo = pi_lo, pi_hi = pi_hi,
    loo_min = min(loo), loo_max = max(loo), loo_sign_flip = as.numeric(flip))
}

res <- comb[, {
  r <- dl_full(logFC, se^2)
  list(k = as.integer(r["k"]), meta_logFC = r["meta_logFC"], ci_lo = r["ci_lo"],
       ci_hi = r["ci_hi"], meta_p = r["meta_p"], I2 = r["I2"], tau2 = r["tau2"],
       pi_lo = r["pi_lo"], pi_hi = r["pi_hi"], loo_min = r["loo_min"],
       loo_max = r["loo_max"], loo_sign_flip = as.logical(r["loo_sign_flip"]))
}, by = gene]
res <- res[k >= 3][order(meta_p)]
cat("genes with k>=3:", nrow(res), "\n")

## ---- assertion against frozen f1_meta.csv ---------------------------------
frz <- fread("analysis/step9_features/f1_meta.csv")
stopifnot(nrow(frz) == nrow(res), setequal(frz$gene, res$gene))
m <- merge(res, frz, by = "gene", suffixes = c("_v2", "_frz"))
stopifnot(all(m$k_v2 == m$k_frz))
for (cl in c("meta_logFC", "ci_lo", "ci_hi", "I2"))
  stopifnot(max(abs(m[[paste0(cl, "_v2")]] - m[[paste0(cl, "_frz")]]), na.rm = TRUE) < 1e-6)
rel_p <- abs(m$meta_p_v2 - m$meta_p_frz) / pmax(m$meta_p_frz, 1e-300)
stopifnot(max(rel_p, na.rm = TRUE) < 1e-4)
cat("ASSERTION PASSED: meta_logFC/CI/p/I2/k reproduce frozen f1_meta.csv for all",
    nrow(m), "genes\n")

dir.create("analysis/step9_features/meta_v2", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "analysis/step9_features/meta_v2/meta_extended.csv")
cat("wrote meta_extended.csv\n")

## key genes for the summary: CLC, IL33 + top consensus genes
cons <- fread("analysis/step6_validation/consensus_genes.csv")[version == "disease_cohorts_10"]
cons[, rank_key := n_sig_FDR005 * 100 + frac_consistent * 10 + pmin(abs(median_logFC), 9)]
top_up <- cons[direction == "up"][order(-rank_key)][1:10]$gene
top_dn <- cons[direction == "down"][order(-rank_key)][1:10]$gene
keys <- unique(c("CLC", "IL33", top_up, top_dn))
kt <- res[gene %in% keys][order(factor(gene, levels = keys))]
fwrite(kt, "analysis/step9_features/meta_v2/key_genes.csv")
print(kt)
cat("flips overall:", sum(res$loo_sign_flip), "of", nrow(res), "\n")
cat("CLC/IL33 rows:\n"); print(res[gene %in% c("CLC", "IL33")])
