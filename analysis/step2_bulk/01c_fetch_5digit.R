#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(GEOquery))
OUT <- "analysis/step2_bulk"
dl <- file.path(OUT, "meta", "matrix")
gdir <- function(gse) { n <- nchar(gse); if (n == 8) paste0("GSE", substr(gse,4,5), "nnn") else paste0("GSE", substr(gse,4,6), "nnn") }
for (g in c("GSE36830","GSE72713","GSE69093","GSE70898","GSE70900")) {
  cat("=====", g, gdir(g), "=====\n")
  dirurl <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", gdir(g), g)
  lst <- system(sprintf('curl -s --max-time 60 "%s"', dirurl), intern = TRUE)
  pat <- paste0(g, "_series_matrix[^\"]*\\.txt\\.gz")
  fns <- unique(unlist(regmatches(lst, gregexpr(pat, lst))))
  cat("files:", paste(fns, collapse = ", "), "\n")
  es <- list()
  for (fn in fns) {
    dest <- file.path(dl, fn)
    if (!file.exists(dest)) system(sprintf('curl -s --max-time 600 -o "%s" "%s%s"', dest, dirurl, fn))
    e <- try(getGEO(filename = dest, AnnotGPL = FALSE), silent = TRUE)
    if (inherits(e, "try-error")) { cat("parse FAILED", fn, "\n"); next }
    es <- c(es, if (is.list(e) && !inherits(e, "ExpressionSet")) e else list(e))
  }
  saveRDS(es, file.path(dl, paste0(g, ".rds")))
  for (j in seq_along(es)) {
    e <- es[[j]]
    cat(sprintf("  [%s] %s: %d x %d\n", g, annotation(e), ncol(e), nrow(e)))
    pd <- pData(e)
    keep <- grep("characteristics|source_name|title$|description|organism", names(pd), value = TRUE)
    write.csv(pd[, intersect(keep, names(pd)), drop = FALSE],
              file.path(OUT, "meta", sprintf("%s_%s_pheno.csv", g, gsub("\\W", "_", annotation(e)))),
              row.names = TRUE)
    for (cf in grep("characteristics", names(pd), value = TRUE)) {
      uv <- unique(as.character(pd[[cf]]))
      cat(sprintf("   %s [%d]: %s\n", cf, length(uv), paste(utils::head(uv, 14), collapse = " | ")))
    }
    m <- exprs(e); q <- quantile(m, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
    cat(sprintf("   vals: min=%.2f p25=%.2f med=%.2f p75=%.2f max=%.2f negfrac=%.3f\n",
                q[1], q[2], q[3], q[4], q[5], mean(m < 0, na.rm = TRUE)))
  }
}
cat("DONE\n")
