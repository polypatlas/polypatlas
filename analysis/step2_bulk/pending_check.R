#!/usr/bin/env Rscript
# pending_check.R — file-level metadata checks for the 11 pending series.
suppressPackageStartupMessages(library(GEOquery))
OUT <- "analysis/step2_bulk"
dl <- file.path(OUT, "meta", "matrix_pending")
dir.create(dl, showWarnings = FALSE, recursive = TRUE)
gdir <- function(gse) { n <- nchar(gse); if (n == 8) paste0("GSE", substr(gse,4,5), "nnn") else paste0("GSE", substr(gse,4,6), "nnn") }
gses <- c("GSE10406","GSE23552","GSE179268","GSE39957","GSE255544","GSE276503",
          "GSE295436","GSE270312","GSE279240","GSE279289","GSE315521")
for (g in gses) {
  cat("=====", g, "=====\n")
  dirurl <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", gdir(g), g)
  lst <- system(sprintf('curl -s --max-time 60 "%s"', dirurl), intern = TRUE)
  pat <- paste0(g, "_series_matrix[^\"]*\\.txt\\.gz")
  fns <- unique(unlist(regmatches(lst, gregexpr(pat, lst))))
  if (!length(fns)) { cat("no matrix files\n"); next }
  for (fn in fns) {
    dest <- file.path(dl, fn)
    if (!file.exists(dest)) system(sprintf('curl -s --max-time 900 -o "%s" "%s%s"', dest, dirurl, fn))
    e <- try(getGEO(filename = dest, AnnotGPL = FALSE), silent = TRUE)
    if (inherits(e, "try-error")) { cat("parse FAILED", fn, "\n"); next }
    es <- if (is.list(e) && !inherits(e, "ExpressionSet")) e else list(e)
    for (j in seq_along(es)) {
      ee <- es[[j]]
      cat(sprintf("  [%s] %d samples x %d features\n", annotation(ee), ncol(ee), nrow(ee)))
      pd <- pData(ee)
      cat("   title[1:3]:", paste(head(pd$title, 3), collapse = " | "), "\n")
      cat("   source[uniq]:", paste(unique(pd$source_name_ch1), collapse = " | "), "\n")
      if ("organism_ch1" %in% names(pd)) cat("   organism:", paste(unique(pd$organism_ch1), collapse = "|"), "\n")
      for (cf in grep("characteristics", names(pd), value = TRUE)) {
        uv <- unique(as.character(pd[[cf]]))
        cat(sprintf("   %s [%d]: %s\n", cf, length(uv), paste(utils::head(uv, 10), collapse = " | ")))
      }
    }
  }
}
cat("DONE\n")
