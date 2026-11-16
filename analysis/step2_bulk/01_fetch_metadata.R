#!/usr/bin/env Rscript
# 01_fetch_metadata.R — download series matrices + sample metadata for step2 cohorts
suppressPackageStartupMessages({
  library(GEOquery)
})
OUT <- "analysis/step2_bulk"
dir.create(file.path(OUT, "meta"), showWarnings = FALSE, recursive = TRUE)
dl <- file.path(OUT, "meta", "geo_dl"); dir.create(dl, showWarnings = FALSE, recursive = TRUE)

gses <- c("GSE338426","GSE198950","GSE136825","GSE36830","GSE194282","GSE158277",
          "GSE72713","GSE179265","GSE189690","GSE267331","GSE69093","GSE313186",
          "GSE70898","GSE70900","GSE174681","GSE150846","GSE330761","GSE213809",
          "GSE270953")

for (g in gses) {
  cat("=====", g, "=====\n")
  f <- file.path(dl, paste0(g, ".rds"))
  ok <- FALSE
  if (!file.exists(f)) {
    es <- try(getGEO(g, destdir = dl, GSEMatrix = TRUE, AnnotGPL = TRUE), silent = TRUE)
    if (inherits(es, "try-error")) { cat("getGEO FAILED:", attr(es, "condition")$message, "\n"); next }
    saveRDS(es, f)
  } else {
    es <- readRDS(f)
  }
  # handle multiple platforms -> list
  plats <- sapply(es, function(x) annotation(x))
  cat("platforms:", paste(plats, collapse = ", "), "\n")
  for (j in seq_along(es)) {
    e <- es[[j]]
    cat(sprintf("  [%s] %s: %d samples x %d features\n", g, annotation(e), ncol(e), nrow(e)))
    pd <- pData(e)
    keep <- grep("characteristics|source_name|title|description|organism|tissue", names(pd), value = TRUE)
    write.csv(pd[, intersect(keep, names(pd)), drop = FALSE],
              file.path(OUT, "meta", sprintf("%s_%s_pheno.csv", g, gsub("\\W", "_", annotation(e)))),
              row.names = TRUE)
    cat("   metadata fields:", paste(intersect(keep, names(pd)), collapse = "; "), "\n")
    ch <- grep("characteristics", names(pd), value = TRUE)
    if (length(ch)) {
      cat("   --- unique values per characteristics field ---\n")
      for (cf in ch) {
        uv <- unique(as.character(pd[[cf]]))
        cat(sprintf("   %s [%d unique]: %s\n", cf, length(uv),
                    paste(utils::head(uv, 12), collapse = " | ")))
      }
    }
    # expression value summary for scale check
    m <- exprs(e)
    q <- quantile(m, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
    cat(sprintf("   value range: min=%.2f p25=%.2f med=%.2f p75=%.2f max=%.2f negfrac=%.3f\n",
                q[1], q[2], q[3], q[4], q[5], mean(m < 0, na.rm = TRUE)))
    fd <- fData(e)
    cat("   fData cols:", paste(head(names(fd), 10), collapse = ", "), "\n")
  }
}
cat("DONE\n")
