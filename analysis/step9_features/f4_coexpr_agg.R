#!/usr/bin/env Rscript
# F4 step B: aggregate per-cohort top-30 pairs -> Fisher-z mean (>=3 cohorts/pair)
# -> per gene TOP-10 positive partners. Writes f4_coexpr.csv + f4_coexpr.js.
suppressMessages(library(data.table))

cpdir <- "analysis/step9_features/checkpoints"
cps <- list.files(cpdir, pattern="^f4_pairs_.*[.]rds$", full.names=TRUE)
cat("checkpoints found:", length(cps), "\n")
stopifnot(length(cps) >= 3)

P <- rbindlist(lapply(cps, readRDS))
cat("total pairs:", nrow(P), "\n")
P[, z := atanh(pmin(r, 0.9999999))]
agg <- P[, .(sumz=sum(z), n_cohorts=.N), by=.(gene, partner)]
agg <- agg[n_cohorts >= 3]
agg[, mean_r := tanh(sumz/n_cohorts)]
setorder(agg, gene, -mean_r)
top10 <- agg[, head(.SD, 10), by=gene]
out <- top10[, .(gene, partner, mean_r=round(mean_r,4), n_cohorts)]
fwrite(out, "analysis/step9_features/f4_coexpr.csv")
cat("f4_coexpr.csv rows:", nrow(out), " genes:", length(unique(out$gene)), "\n")

con <- file("analysis/step9_features/f4_coexpr.js", "w")
writeLines("var F4COEXPR = {", con)
gs <- unique(out$gene); n <- length(gs)
for (i in seq_along(gs)) {
  d <- out[gene==gs[i]]
  inner <- paste(sprintf('["%s",%s,%d]', d$partner, d$mean_r, d$n_cohorts), collapse=",")
  writeLines(sprintf('"%s":[%s]%s', gs[i], inner, if (i<n) "," else ""), con)
}
writeLines("};", con); close(con)
cat("wrote f4_coexpr.js\n")
