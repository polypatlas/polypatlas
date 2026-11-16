#!/usr/bin/env Rscript
# 01b_fetch_metadata.R — direct HTTPS matrix download + local parse (replaces 01)
suppressPackageStartupMessages({ library(GEOquery) })
OUT <- "analysis/step2_bulk"
dl <- file.path(OUT, "meta", "matrix"); dir.create(dl, showWarnings = FALSE, recursive = TRUE)
gdir <- function(gse) paste0("GSE", substr(gse, 4, 6), "nnn")

gses <- c("GSE338426","GSE198950","GSE136825","GSE36830","GSE194282","GSE158277",
          "GSE72713","GSE179265","GSE189690","GSE267331","GSE69093","GSE313186",
          "GSE70898","GSE70900","GSE174681","GSE150846","GSE330761","GSE213809",
          "GSE270953")

for (g in gses) {
  cat("=====", g, "=====\n")
  rds <- file.path(dl, paste0(g, ".rds"))
  if (file.exists(rds)) { cat("cached\n"); es <- readRDS(rds) }
  else {
    dirurl <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", gdir(g), g)
    lst <- try(system(sprintf("curl -s --max-time 60 '%s'", dirurl), intern = TRUE), silent = TRUE)
    if (inherits(lst, "try-error") || !length(lst)) { cat("DIR LISTING FAILED\n"); next }
    fns <- regmatches(lst, gregexpr(paste0(g, "_series_matrix[^\"<>]*\\.txt\\.gz"), lst))
    fns <- unique(unlist(fns))
    if (!length(fns)) { cat("no matrix files found\n"); next }
    cat("files:", paste(fns, collapse = ", "), "\n")
    es <- list()
    for (fn in fns) {
      dest <- file.path(dl, fn)
      if (!file.exists(dest))
        system(sprintf("curl -s --max-time 600 -o '%s' '%s%s'", dest, dirurl, fn))
      e <- try(getGEO(filename = dest, AnnotGPL = FALSE), silent = TRUE)
      if (inherits(e, "try-error")) { cat("parse FAILED for", fn, "\n"); next }
      es <- c(es, if (is.list(e) && !inherits(e, "ExpressionSet")) e else list(e))
    }
    if (!length(es)) { cat("ALL PARSE FAILED\n"); next }
    saveRDS(es, rds)
  }
  for (j in seq_along(es)) {
    e <- es[[j]]
    cat(sprintf("  [%s] %s: %d samples x %d features\n", g, annotation(e), ncol(e), nrow(e)))
    pd <- pData(e)
    keep <- grep("characteristics|source_name|title$|description|organism", names(pd), value = TRUE)
    write.csv(pd[, intersect(keep, names(pd)), drop = FALSE],
              file.path(OUT, "meta", sprintf("%s_%s_pheno.csv", g, gsub("\\W", "_", annotation(e)))),
              row.names = TRUE)
    ch <- grep("characteristics", names(pd), value = TRUE)
    for (cf in ch) {
      uv <- unique(as.character(pd[[cf]]))
      cat(sprintf("   %s [%d]: %s\n", cf, length(uv), paste(utils::head(uv, 14), collapse = " | ")))
    }
    m <- exprs(e)
    q <- quantile(m, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
    cat(sprintf("   vals: min=%.3f p25=%.3f med=%.3f p75=%.3f max=%.3f negfrac=%.4f int=%s\n",
                q[1], q[2], q[3], q[4], q[5], mean(m < 0, na.rm = TRUE),
                all(m == round(m), na.rm = TRUE)))
  }
}
cat("DONE\n")
