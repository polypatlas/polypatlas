#!/usr/bin/env Rscript
# donor_meta.R — GSE36830 donor sensitivity meta-analyses.
# dl_full() verbatim from step9_features/meta_v2/meta_v2.R (as reused in
# sensitivity_strict_tissue/strict_meta.R): DL random effects + tau2 +
# 95% PI (t, df=k-1) + LOO; se = |logFC/t|; k>=3; BH-FDR genome-wide.
# Scenarios:
#   A_main9_no36830    : main 10 cohorts minus GSE36830 (k<=9)
#   A_strict6_no36830  : strict 7 cohorts minus GSE36830 (k<=6)
#   B_main10_polypOnly : main 10, GSE36830 = polyp-only DE table (k<=10)
suppressPackageStartupMessages(library(data.table))
setwd(".")
OUT <- "analysis/step2_bulk/sensitivity_gse36830_donor"
DEDIR <- "analysis/step2_bulk/de"

primary <- list(
  GSE136825     = "GSE136825_polyp_vs_ctrl",
  GSE338426     = "GSE338426_CRSwNP_vs_ctrl",
  GSE179265     = "GSE179265_CRSwNP_vs_ctrl",
  GSE36830      = "GSE36830_CRSwNP_vs_UTctrl",
  GSE194282     = "GSE194282_NP_vs_ctrl",
  GSE69093      = "GSE69093_case_vs_ctrl_celltypeAdj",
  GSE313186     = "GSE313186_exposed_vs_ctrl",
  GSE70898_70900 = "GSE70898_70900_polypT_vs_nasalT",
  GSE198950     = "GSE198950_CRSwNP_vs_NonCRS",
  GSE72713      = "GSE72713_CRSwNP_vs_ctrl"
)
strict7 <- setdiff(names(primary), c("GSE69093", "GSE313186", "GSE70898_70900"))

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

run_meta <- function(scenario, cohorts, polyp_only = FALSE) {
  all_list <- lapply(cohorts, function(gse) {
    f <- if (polyp_only && gse == "GSE36830")
      file.path(OUT, "GSE36830_CRSwNP_polypOnly_vs_UTctrl.csv")
    else file.path(DEDIR, paste0(primary[[gse]], ".csv"))
    d <- fread(f, select = c("symbol", "logFC", "t"))
    d <- d[!is.na(t) & t != 0 & !is.na(logFC) & symbol != ""]
    d[, cohort := gse]
    d[, se := abs(logFC / t)]
    d[se < 1e-8, se := 1e-8]
    unique(d[, .(gene = symbol, cohort, logFC, se)], by = c("gene", "cohort"))
  })
  comb <- rbindlist(all_list)
  cat(sprintf("[%s] combined rows: %d unique genes: %d\n", scenario, nrow(comb), length(unique(comb$gene))))
  res <- comb[, {
    r <- dl_full(logFC, se^2)
    list(k = as.integer(r["k"]), meta_logFC = r["meta_logFC"], ci_lo = r["ci_lo"],
         ci_hi = r["ci_hi"], meta_p = r["meta_p"], I2 = r["I2"], tau2 = r["tau2"],
         pi_lo = r["pi_lo"], pi_hi = r["pi_hi"], loo_min = r["loo_min"],
         loo_max = r["loo_max"], loo_sign_flip = as.logical(r["loo_sign_flip"]))
  }, by = gene]
  res <- res[k >= 3][order(meta_p)]
  res[, FDR := p.adjust(meta_p, method = "BH")]
  res[, scenario := scenario]
  fwrite(res, file.path(OUT, paste0("meta_extended_", scenario, ".csv")))
  sig <- res[FDR < 0.05]
  pi_excl0 <- sig[(pi_lo > 0 & pi_hi > 0) | (pi_lo < 0 & pi_hi < 0)]
  gs <- data.table(scenario = scenario, n_cohorts = length(cohorts),
                   n_genes_k_ge3 = nrow(res), n_fdr_sig = nrow(sig),
                   n_pi_excl0_among_sig = nrow(pi_excl0),
                   frac_pi_excl0_among_sig = if (nrow(sig)) nrow(pi_excl0) / nrow(sig) else NA_real_,
                   loo_flips_all = sum(res$loo_sign_flip), loo_flips_sig = sum(sig$loo_sign_flip))
  cat(sprintf("[%s] FDR<0.05: %d; PI excl 0: %d (%.4f)\n",
              scenario, nrow(sig), nrow(pi_excl0), gs$frac_pi_excl0_among_sig))
  list(res = res, gs = gs)
}

scen <- list(
  A_main9_no36830    = list(cohorts = setdiff(names(primary), "GSE36830"), polyp_only = FALSE),
  A_strict6_no36830  = list(cohorts = setdiff(strict7, "GSE36830"), polyp_only = FALSE),
  B_main10_polypOnly = list(cohorts = names(primary), polyp_only = TRUE)
)
outs <- lapply(names(scen), function(s) run_meta(s, scen[[s]]$cohorts, scen[[s]]$polyp_only))
names(outs) <- names(scen)
fwrite(rbindlist(lapply(outs, `[[`, "gs")), file.path(OUT, "meta_global_stats.csv"))

## key genes CLC / IL33 across scenarios + frozen references
keys <- rbindlist(lapply(names(outs), function(s) outs[[s]]$res[gene %in% c("CLC", "IL33")]))
froz_main_all <- fread("analysis/step9_features/meta_v2/meta_extended.csv")
froz_main_all[, FDR := p.adjust(meta_p, method = "BH")]  # frozen file carries no FDR column; BH recomputed from its frozen meta_p over all genes
froz_main <- froz_main_all[gene %in% c("CLC", "IL33")]
froz_main[, scenario := "frozen_main10"]
froz_strict <- fread("analysis/step2_bulk/sensitivity_strict_tissue/strict_meta_extended.csv")[gene %in% c("CLC", "IL33")]
froz_strict[, scenario := "frozen_strict7"]
cols <- c("scenario", "gene", "k", "meta_logFC", "ci_lo", "ci_hi", "meta_p", "FDR",
          "I2", "tau2", "pi_lo", "pi_hi", "loo_min", "loo_max", "loo_sign_flip")
allk <- rbindlist(list(keys[, ..cols], froz_main[, ..cols], froz_strict[, ..cols]), use.names = TRUE)
setorder(allk, gene, scenario)
fwrite(allk, file.path(OUT, "meta_key_genes.csv"))
print(allk)
cat("wrote meta_extended_*.csv + meta_global_stats.csv + meta_key_genes.csv\n")
