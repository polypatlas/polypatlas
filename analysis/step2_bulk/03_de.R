#!/usr/bin/env Rscript
# 03_de.R — generic DE driver over prep/<GSE>.rds.
# Writes de/<GSE>_<contrast>.csv + runlog.csv (appended).
suppressPackageStartupMessages({ library(limma); library(DESeq2) })
OUT <- "analysis/step2_bulk"
dir.create(file.path(OUT, "de"), showWarnings = FALSE, recursive = TRUE)
RL <- file.path(OUT, "runlog.csv")
if (!file.exists(RL)) {
  writeLines("gse,contrast,n_case,n_control,scale_verdict,method,status,note", RL)
}
done <- if (file.exists(RL)) {
  d <- read.csv(RL); paste(d$gse, d$contrast)
} else character(0)

collapse_sym <- function(res, symbols) {
  # res: data.frame with rownames = features, cols logFC,P.Value,adj.P.Val,(AveExpr|baseMean|t)
  res$symbol <- symbols[rownames(res)]
  res <- res[!is.na(res$symbol) & res$symbol != "", ]
  res <- res[order(res$P.Value), ]
  res <- res[!duplicated(res$symbol), ]
  res[order(res$P.Value), ]
}

run_deseq2 <- function(mat, grp) {
  mat <- mat[, grp %in% c("case", "control")]; grp <- grp[grp %in% c("case", "control")]
  keep <- rowSums(mat) >= 10
  mat <- mat[keep, ]
  dds <- DESeqDataSetFromMatrix(countData = round(mat),
                                colData = data.frame(group = factor(grp, levels = c("control", "case"))),
                                design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("group", "case", "control")))
  data.frame(logFC = r$log2FoldChange, AveExpr = log2(r$baseMean + 1), t = r$stat,
             P.Value = r$pvalue, adj.P.Val = r$padj, row.names = rownames(r))
}

run_limma <- function(mat, grp, cov = NULL) {
  mat <- mat[, grp %in% c("case", "control")]; grp <- factor(grp[grp %in% c("case", "control")], levels = c("control", "case"))
  if (!is.null(cov)) { cov <- cov[names(grp)]; design <- model.matrix(~ grp + cov) }
  else design <- model.matrix(~ grp)
  colnames(design) <- make.names(colnames(design))
  fit <- eBayes(lmFit(mat, design))
  tt <- topTable(fit, coef = grep("grpcase", colnames(design)), number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, row.names = rownames(tt))
}

run_limma_paired <- function(mat, meta, arm) {
  m <- meta[meta$treatment == arm & meta$time %in% c("Pre-Treatment", "Post-Treatment"), ]
  # keep only subjects with both timepoints
  tab <- table(m$subject, m$time)
  oks <- rownames(tab)[tab[, "Pre-Treatment"] > 0 & tab[, "Post-Treatment"] > 0]
  m <- m[m$subject %in% oks, ]
  mat <- mat[, rownames(m)]
  m$subject <- factor(m$subject); m$time <- factor(m$time, levels = c("Pre-Treatment", "Post-Treatment"))
  design <- model.matrix(~ subject + time, data = m)
  colnames(design) <- make.names(colnames(design))
  fit <- eBayes(lmFit(mat, design))
  tt <- topTable(fit, coef = "timePost.Treatment", number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, row.names = rownames(tt))
}

run_voom <- function(mat, grp) {
  mat <- mat[, grp %in% c("case", "control")]; grp <- factor(grp[grp %in% c("case", "control")], levels = c("control", "case"))
  design <- model.matrix(~ grp); colnames(design) <- make.names(colnames(design))
  v <- voom(mat, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, row.names = rownames(tt))
}

do_contrast <- function(sp, gse, cname, case_ids, ctrl_ids, method, meta_extra = NULL) {
  key <- paste(gse, cname)
  outfile <- file.path(OUT, "de", sprintf("%s_%s.csv", gse, cname))
  if (file.exists(outfile)) { cat(sprintf("[skip] %s exists\n", key)); return(invisible(NULL)) }
  if (key %in% done) { cat(sprintf("[skip] %s done\n", key)); return(invisible(NULL)) }
  mat <- sp$expr; sym <- sp$symbols
  if (anyDuplicated(rownames(mat))) {  # e.g. GSE72713 gene-name rownames: uniquify
    rownames(mat) <- make.unique(rownames(mat))
    sym <- unname(sym); names(sym) <- rownames(mat)
  }
  if (method == "limma_paired") {
    arm <- meta_extra$paired_arm
    res <- run_limma_paired(mat, sp$meta, arm)
    n_case <- sum(sp$meta$treatment == arm & sp$meta$time == "Post-Treatment")
    n_ctrl <- sum(sp$meta$treatment == arm & sp$meta$time == "Pre-Treatment")
  } else {
    ids <- c(case_ids, ctrl_ids)
    grp <- ifelse(colnames(mat) %in% case_ids, "case", "control")
    names(grp) <- colnames(mat); grp <- grp[ids]
    if (method == "deseq2") res <- run_deseq2(mat[, ids], grp)
    else if (method == "limma_cov") {
      cov <- sp$meta[ids, "celltype"]; names(cov) <- ids
      if (length(unique(cov)) < 2) res <- run_limma(mat[, ids], grp)  # single-level: drop covariate
      else res <- run_limma(mat[, ids], grp, cov)
    }
    else if (method == "limma_voom") res <- run_voom(mat[, ids], grp)
    else res <- run_limma(mat[, ids], grp)
    n_case <- length(case_ids); n_ctrl <- length(ctrl_ids)
  }
  res <- collapse_sym(res, sym)
  status <- if (n_case < 3 || n_ctrl < 3) "LOW_POWER" else "OK"
  out <- data.frame(symbol = res$symbol, logFC = res$logFC, AveExpr = res$AveExpr,
                    t = res$t, p = res$P.Value, FDR = res$adj.P.Val)
  write.csv(out, file.path(OUT, "de", sprintf("%s_%s.csv", gse, cname)), row.names = FALSE)
  sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("[de] %s %s: n=%d genes, %d FDR<0.05 | case=%d ctrl=%d | %s\n",
              gse, cname, nrow(out), sig, n_case, n_ctrl, status))
  line <- sprintf("%s,%s,%d,%d,%s,%s,%s,\"%s (%d FDR<0.05)\"", gse, cname, n_case, n_ctrl,
                  sp$scale_verdict, method, status, sp$note, sig)
  cat(line, "\n", file = RL, append = TRUE)
  invisible(NULL)
}

specs <- list(
  list(gse = "GSE338426", primary = list(name = "CRSwNP_vs_ctrl"),
       extras = c("eCRSwNP_vs_ctrl", "neCRSwNP_vs_ctrl")),
  list(gse = "GSE198950", primary = list(name = "CRSwNP_vs_NonCRS"),
       extras = c("CRSsNP_vs_NonCRS", "EOS_high_vs_low")),
  list(gse = "GSE136825", primary = list(name = "polyp_vs_ctrl"), extras = NULL),
  list(gse = "GSE36830", primary = list(name = "CRSwNP_vs_UTctrl"),
       extras = c("CRSsNP_vs_UTctrl")),
  list(gse = "GSE194282", primary = list(name = "NP_vs_ctrl"),
       extras = c("elderly_vs_nonelderly")),
  list(gse = "GSE72713", primary = list(name = "CRSwNP_vs_ctrl"), extras = c("e_vs_ne")),
  list(gse = "GSE179265", primary = list(name = "CRSwNP_vs_ctrl"),
       extras = c("US_vs_Jpn_polyps")),
  list(gse = "GSE189690", primary = list(name = "NERD_vs_nonNERD"), extras = NULL),
  list(gse = "GSE267331", primary = list(name = "dupilumab_post_vs_pre", paired_arm = "Dupilumab"),
       extras = list(list(name = "placebo_post_vs_pre", paired_arm = "Placebo"))),
  list(gse = "GSE69093", primary = list(name = "case_vs_ctrl_celltypeAdj"),
       extras = c("epi_case_vs_ctrl", "fib_case_vs_ctrl")),
  list(gse = "GSE313186", primary = list(name = "exposed_vs_ctrl"), extras = NULL),
  list(gse = "GSE70898_70900", primary = list(name = "polypT_vs_nasalT"), extras = NULL),
  list(gse = "GSE174681", primary = list(name = "NP_ASC_vs_PBLnaive"),
       extras = c("NP_naive_vs_PBLnaive")),
  list(gse = "GSE150846", primary = list(name = "NP_ILC2_vs_blood_ILC2"),
       extras = c("type2stim_vs_baseline")),
  list(gse = "GSE330761", primary = list(name = "CRSwNP_MCTC_vs_CRSsNP_MCTC"), extras = NULL),
  list(gse = "GSE213809", primary = list(name = "NP_vs_PB_stimulated"),
       extras = c("EP_vs_EN_NPstim"))
)

for (spc in specs) {
  gse <- spc$gse
  sp <- readRDS(file.path(OUT, "prep", paste0(gse, ".rds")))
  method <- sp$method
  # primary
  pr <- sp$primary
  if (!is.null(pr$paired_arm)) {
    do_contrast(sp, gse, pr$name, NULL, NULL, method, meta_extra = pr)
  } else {
    case_ids <- colnames(sp$expr)[sp$meta$group == "case"]
    ctrl_ids <- colnames(sp$expr)[sp$meta$group == "control"]
    do_contrast(sp, gse, pr$name, case_ids, ctrl_ids, method)
  }
  # extras: resolve from sp$extras by name
  exs <- sp$extras
  want <- spc$extras
  if (!is.null(want) && length(exs)) {
    for (w in want) {
      wname <- if (is.list(w)) w$name else w
      hit <- Filter(function(e) e$name == wname, exs)
      if (!length(hit)) { cat(sprintf("[warn] %s extra %s not in prep\n", gse, wname)); next }
      e <- hit[[1]]
      if (!is.null(e$paired_arm)) do_contrast(sp, gse, e$name, NULL, NULL, method, meta_extra = e)
      else do_contrast(sp, gse, e$name, e$case_vals, e$ctrl_vals, method)
    }
  }
}
cat("ALL DE DONE\n")
