#!/usr/bin/env Rscript
# F1: cross-cohort random-effects meta-analysis (DerSimonian-Laird) per gene
# over the PRIMARY disease-vs-control contrast of each bulk cohort.
library(data.table)

# --- primary contrast mapping (disease vs control; from runlog.csv inspection) ---
primary <- list(
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
# Excluded (no disease-vs-control primary contrast):
#   GSE189690 (NERD_vs_nonNERD: subtype contrast, no healthy control)
#   GSE267331 (dupilumab/placebo post_vs_pre: treatment contrasts, SINUS-52)
#   GSE174681 (sorted B-cell subsets: atlas design)
#   GSE150846 (sorted ILC2: tissue-of-origin/stimulation design)
#   GSE213809 (sorted ILC2 stimulation design)
#   GSE330761 (CRSwNP_MCTC vs CRSsNP_MCTC: no healthy control)

all_list <- lapply(names(primary), function(gse) {
  f <- file.path("analysis/step2_bulk/de", paste0(primary[[gse]], ".csv"))
  d <- fread(f, select=c("symbol","logFC","t"))
  d <- d[!is.na(t) & t != 0 & !is.na(logFC) & symbol != ""]
  d[, cohort := gse]
  d[, se := abs(logFC/t)]
  d[se < 1e-8, se := 1e-8]   # floor: logFC=0 -> se=0 would give infinite weight
  unique(d[, .(gene=symbol, cohort, logFC, se)], by=c("gene","cohort"))
})
comb <- rbindlist(all_list)
cat("combined rows:", nrow(comb), " unique genes:", length(unique(comb$gene)), "\n")

fwrite(comb[, .(gene, cohort, logFC, se)], "analysis/step9_features/f1_meta_cohorts.csv")
cat("wrote f1_meta_cohorts.csv\n")

# DerSimonian-Laird random-effects meta per gene (k >= 3)
dl_meta <- function(y, v) {
  k <- length(y)
  w <- 1/v
  mu_fe <- sum(w*y)/sum(w)
  Q <- sum(w*(y-mu_fe)^2)
  df <- k-1
  tau2 <- max(0, (Q-df)/(sum(w)-sum(w^2)/sum(w)))
  w2 <- 1/(v+tau2)
  mu <- sum(w2*y)/sum(w2)
  se_mu <- sqrt(1/sum(w2))
  z <- mu/se_mu
  p <- 2*pnorm(-abs(z))
  ci_lo <- mu - 1.96*se_mu; ci_hi <- mu + 1.96*se_mu
  I2 <- if (Q>0) max(0,(Q-df)/Q)*100 else 0
  c(k=k, meta_logFC=mu, ci_lo=ci_lo, ci_hi=ci_hi, meta_p=p, I2=I2)
}

res <- comb[, {
  r <- dl_meta(logFC, se^2)
  list(k=as.integer(r["k"]), meta_logFC=r["meta_logFC"], ci_lo=r["ci_lo"],
       ci_hi=r["ci_hi"], meta_p=r["meta_p"], I2=r["I2"])
}, by=gene]
res <- res[k >= 3][order(meta_p)]
cat("genes with k>=3:", nrow(res), "\n")
fwrite(res, "analysis/step9_features/f1_meta.csv")
cat("wrote f1_meta.csv\n")

# CLC sanity check
clc <- res[gene=="CLC"]; il33 <- res[gene=="IL33"]
cat("CLC: ", paste(round(unlist(clc[, -1]),4), collapse=" "), "\n")
cat("IL33:", paste(round(unlist(il33[, -1]),4), collapse=" "), "\n")

# compact JS: all k>=3 genes, rounded
js <- res[, .(gene,
  lfc=round(meta_logFC,3), lo=round(ci_lo,3), hi=round(ci_hi,3),
  p=signif(meta_p,3), i2=round(I2,1), k=as.integer(k))]
con <- file("analysis/step9_features/f1_meta.js", "w")
writeLines("var F1META = {", con)
n <- nrow(js)
for (i in seq_len(n)) {
  r <- js[i]
  ln <- sprintf('"%s":{lfc:%s,lo:%s,hi:%s,p:%s,i2:%s,k:%d}%s',
    r$gene, r$lfc, r$lo, r$hi, format(r$p, scientific=FALSE, trim=TRUE), r$i2, r$k,
    if (i<n) "," else "")
  writeLines(ln, con)
}
writeLines("};", con); close(con)
cat("wrote f1_meta.js\n")
