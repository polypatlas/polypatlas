#!/usr/bin/env Rscript
# strict_meta.R — strict-tissue sensitivity meta-analysis (supplement B).
# Reuses dl_full() verbatim from analysis/step9_features/meta_v2/meta_v2.R
# (DL random effects + tau2 + 95% PI (t, df=k-1) + LOO; se = |logFC/t| from the
# same frozen de/ tables), restricted to the 7 whole-tissue cohorts
# (excludes GSE69093, GSE313186, GSE70898_70900). k>=3.
suppressPackageStartupMessages(library(data.table))
setwd(".")
OUT <- "analysis/step2_bulk/sensitivity_strict_tissue"

primary <- list(
  GSE136825     = "GSE136825_polyp_vs_ctrl",
  GSE338426     = "GSE338426_CRSwNP_vs_ctrl",
  GSE179265     = "GSE179265_CRSwNP_vs_ctrl",
  GSE36830      = "GSE36830_CRSwNP_vs_UTctrl",
  GSE194282     = "GSE194282_NP_vs_ctrl",
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

dl_full <- function(y, v) {   # verbatim from meta_v2.R
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
  tcrit <- qt(0.975, df = k - 1)
  pi_lo <- mu - tcrit * sqrt(tau2 + se_mu^2)
  pi_hi <- mu + tcrit * sqrt(tau2 + se_mu^2)
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
res[, FDR := p.adjust(meta_p, method = "BH")]
cat("genes with k>=3:", nrow(res), "\n")
fwrite(res, file.path(OUT, "strict_meta_extended.csv"))

sig <- res[FDR < 0.05]
pi_excl0 <- sig[(pi_lo > 0 & pi_hi > 0) | (pi_lo < 0 & pi_hi < 0)]
cat("BH-FDR<0.05 genes:", nrow(sig), "\n")
cat("of which 95% PI excludes 0:", nrow(pi_excl0),
    sprintf("(%.4f)\n", nrow(pi_excl0) / nrow(sig)))
cat("LOO flips overall:", sum(res$loo_sign_flip), "of", nrow(res), "\n")
cat("LOO flips among FDR<0.05:", sum(sig$loo_sign_flip), "of", nrow(sig), "\n")

## key genes: CLC, IL33 + top-10 up/down of the strict consensus (same rank key as meta_v2)
scons <- fread(file.path(OUT, "strict_consensus.csv"))
scons[, rank_key := n_sig_FDR005 * 100 + frac_consistent * 10 + pmin(abs(median_logFC), 9)]
top_up <- scons[direction == "up"][order(-rank_key)][1:10]$gene
top_dn <- scons[direction == "down"][order(-rank_key)][1:10]$gene
keys <- unique(c("CLC", "IL33", top_up, top_dn))
kt <- res[gene %in% keys][order(factor(gene, levels = keys))]
fwrite(kt, file.path(OUT, "strict_key_genes.csv"))
print(kt)
cat("CLC/IL33 rows:\n"); print(res[gene %in% c("CLC", "IL33")])

## global stats sidecar for the summary writer
gs <- data.table(
  n_genes_k_ge3 = nrow(res), n_fdr_sig = nrow(sig),
  n_pi_excl0_among_sig = nrow(pi_excl0),
  frac_pi_excl0_among_sig = nrow(pi_excl0) / nrow(sig),
  loo_flips_all = sum(res$loo_sign_flip), loo_flips_sig = sum(sig$loo_sign_flip))
fwrite(gs, file.path(OUT, "strict_meta_global_stats.csv"))
cat("wrote strict_meta_extended.csv + strict_key_genes.csv + strict_meta_global_stats.csv\n")
