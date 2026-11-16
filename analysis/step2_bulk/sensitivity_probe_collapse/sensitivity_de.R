#!/usr/bin/env Rscript
# sensitivity_de.R — probe-collapse sensitivity analysis for PolypAtlas step2.
# IDENTICAL pipeline to ../03_de.R (same prep inputs, same designs, same limma
# fits) with ONE variable changed: probe->gene collapse rule.
#   main rule  : collapse_sym()      = per symbol, probe with MINIMUM p-value
#   sens rule  : collapse_maxvar()   = per symbol, probe with MAXIMUM VARIANCE
#                (variance across the samples entering that contrast's fit,
#                 computed on the prep expression scale; pre-specified, does
#                 not use any test statistic)
# Scope: microarray cohorts only (GPL570 / Illumina HT-12 prep units):
#   GSE36830, GSE194282, GSE267331, GSE69093, GSE70898_70900, GSE150846
# The min-p collapse is recomputed in-memory and ASSERTED equal to the frozen
# de/*.csv outputs; any mismatch stops the run (pipeline-identity check).
suppressPackageStartupMessages({ library(limma) })
BASE <- "analysis/step2_bulk"
OUT  <- file.path(BASE, "sensitivity_probe_collapse")
dir.create(file.path(OUT, "de_maxvar"), showWarnings = FALSE, recursive = TRUE)

## ---- collapse rules -------------------------------------------------------
collapse_sym <- function(res, symbols) {           # verbatim from 03_de.R
  res$symbol <- symbols[rownames(res)]
  res <- res[!is.na(res$symbol) & res$symbol != "", ]
  res <- res[order(res$P.Value), ]
  res <- res[!duplicated(res$symbol), ]
  res[order(res$P.Value), ]
}
row_vars <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowVars(m)
  else apply(m, 1, var)
}
collapse_maxvar <- function(res, expr, samples, symbols) {
  res$symbol <- symbols[rownames(res)]
  res <- res[!is.na(res$symbol) & res$symbol != "", ]
  v <- row_vars(expr[rownames(res), samples, drop = FALSE])
  res$probe_var <- v[rownames(res)]
  res <- res[order(-res$probe_var, rownames(res)), ]   # deterministic tie-break
  res <- res[!duplicated(res$symbol), ]
  res[order(res$symbol), ]
}

## ---- fit functions (copied from 03_de.R; return fit samples too) ----------
run_limma <- function(mat, grp, cov = NULL) {
  mat <- mat[, grp %in% c("case", "control")]; grp <- factor(grp[grp %in% c("case", "control")], levels = c("control", "case"))
  if (!is.null(cov)) { cov <- cov[names(grp)]; design <- model.matrix(~ grp + cov) }
  else design <- model.matrix(~ grp)
  colnames(design) <- make.names(colnames(design))
  fit <- eBayes(lmFit(mat, design))
  tt <- topTable(fit, coef = grep("grpcase", colnames(design)), number = Inf, sort.by = "none")
  list(res = data.frame(logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
                        P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, row.names = rownames(tt)),
       samples = colnames(mat))
}
run_limma_paired <- function(mat, meta, arm) {
  m <- meta[meta$treatment == arm & meta$time %in% c("Pre-Treatment", "Post-Treatment"), ]
  tab <- table(m$subject, m$time)
  oks <- rownames(tab)[tab[, "Pre-Treatment"] > 0 & tab[, "Post-Treatment"] > 0]
  m <- m[m$subject %in% oks, ]
  mat <- mat[, rownames(m)]
  m$subject <- factor(m$subject); m$time <- factor(m$time, levels = c("Pre-Treatment", "Post-Treatment"))
  design <- model.matrix(~ subject + time, data = m)
  colnames(design) <- make.names(colnames(design))
  fit <- eBayes(lmFit(mat, design))
  tt <- topTable(fit, coef = "timePost.Treatment", number = Inf, sort.by = "none")
  list(res = data.frame(logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
                        P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, row.names = rownames(tt)),
       samples = rownames(m))
}

## ---- driver ----------------------------------------------------------------
comp_rows <- list()
do_contrast_sens <- function(sp, gse, cname, case_ids, ctrl_ids, meta_extra = NULL) {
  mat <- sp$expr; sym <- sp$symbols
  if (anyDuplicated(rownames(mat))) {
    rownames(mat) <- make.unique(rownames(mat)); sym <- unname(sym); names(sym) <- rownames(mat)
  }
  method <- sp$method
  if (method == "limma_paired") {
    fr <- run_limma_paired(mat, sp$meta, meta_extra$paired_arm)
  } else {
    ids <- c(case_ids, ctrl_ids)
    grp <- ifelse(colnames(mat) %in% case_ids, "case", "control")
    names(grp) <- colnames(mat); grp <- grp[ids]
    if (method == "limma_cov") {
      cov <- sp$meta[ids, "celltype"]; names(cov) <- ids
      if (length(unique(cov)) < 2) fr <- run_limma(mat[, ids], grp) else fr <- run_limma(mat[, ids], grp, cov)
    } else fr <- run_limma(mat[, ids], grp)
  }
  res <- fr$res
  ## assertion: min-p collapse must reproduce the frozen de/ output exactly
  main <- collapse_sym(res, sym)
  frozen <- read.csv(file.path(BASE, "de", sprintf("%s_%s.csv", gse, cname)), stringsAsFactors = FALSE)
  stopifnot(setequal(main$symbol, frozen$symbol))
  mm <- merge(data.frame(symbol = main$symbol, logFC = main$logFC, FDR = main$adj.P.Val),
              frozen, by = "symbol", suffixes = c("_new", "_frz"))
  stopifnot(nrow(mm) == nrow(frozen),
            max(abs(mm$logFC_new - mm$logFC_frz), na.rm = TRUE) < 1e-6,
            max(abs(mm$FDR_new - mm$FDR), na.rm = TRUE) < 1e-6)
  ## sensitivity collapse
  sv <- collapse_maxvar(res, mat, fr$samples, sym)
  out <- data.frame(symbol = sv$symbol, logFC = sv$logFC, AveExpr = sv$AveExpr,
                    t = sv$t, p = sv$P.Value, FDR = sv$adj.P.Val)
  write.csv(out, file.path(OUT, "de_maxvar", sprintf("%s_%s.csv", gse, cname)), row.names = FALSE)
  deg_main <- main$symbol[!is.na(main$adj.P.Val) & main$adj.P.Val < 0.05]
  deg_sv   <- out$symbol[!is.na(out$FDR) & out$FDR < 0.05]
  comp_rows[[length(comp_rows) + 1]] <<- data.frame(
    gse = gse, contrast = cname, n_genes_main = nrow(frozen), n_genes_maxvar = nrow(out),
    deg_main_minp = length(deg_main), deg_maxvar = length(deg_sv),
    deg_intersection = length(intersect(deg_main, deg_sv)),
    deg_union = length(union(deg_main, deg_sv)), stringsAsFactors = FALSE)
  cat(sprintf("[sens] %s %s: genes %d->%d | DEG minp=%d maxvar=%d overlap=%d (minp refit == frozen de/ OK)\n",
              gse, cname, nrow(frozen), nrow(out), length(deg_main), length(deg_sv),
              length(intersect(deg_main, deg_sv))))
  invisible(NULL)
}

specs <- list(
  list(gse = "GSE36830", primary = list(name = "CRSwNP_vs_UTctrl"), extras = c("CRSsNP_vs_UTctrl")),
  list(gse = "GSE194282", primary = list(name = "NP_vs_ctrl"), extras = c("elderly_vs_nonelderly")),
  list(gse = "GSE267331", primary = list(name = "dupilumab_post_vs_pre", paired_arm = "Dupilumab"),
       extras = list(list(name = "placebo_post_vs_pre", paired_arm = "Placebo"))),
  list(gse = "GSE69093", primary = list(name = "case_vs_ctrl_celltypeAdj"),
       extras = c("epi_case_vs_ctrl", "fib_case_vs_ctrl")),
  list(gse = "GSE70898_70900", primary = list(name = "polypT_vs_nasalT"), extras = NULL),
  list(gse = "GSE150846", primary = list(name = "NP_ILC2_vs_blood_ILC2"), extras = c("type2stim_vs_baseline"))
)
for (spc in specs) {
  gse <- spc$gse
  sp <- readRDS(file.path(BASE, "prep", paste0(gse, ".rds")))
  pr <- sp$primary
  if (!is.null(pr$paired_arm)) do_contrast_sens(sp, gse, pr$name, NULL, NULL, meta_extra = pr)
  else do_contrast_sens(sp, gse, pr$name,
                        colnames(sp$expr)[sp$meta$group == "case"],
                        colnames(sp$expr)[sp$meta$group == "control"])
  exs <- sp$extras; want <- spc$extras
  if (!is.null(want) && length(exs)) {
    for (w in want) {
      wname <- if (is.list(w)) w$name else w
      hit <- Filter(function(e) e$name == wname, exs)
      if (!length(hit)) { cat(sprintf("[warn] %s extra %s not in prep\n", gse, wname)); next }
      e <- hit[[1]]
      if (!is.null(e$paired_arm)) do_contrast_sens(sp, gse, e$name, NULL, NULL, meta_extra = e)
      else do_contrast_sens(sp, gse, e$name, e$case_vals, e$ctrl_vals)
    }
  }
  rm(sp); gc(verbose = FALSE)
}
comp <- do.call(rbind, comp_rows)
write.csv(comp, file.path(OUT, "deg_comparison.csv"), row.names = FALSE)
cat("ALL SENSITIVITY DE DONE\n")
