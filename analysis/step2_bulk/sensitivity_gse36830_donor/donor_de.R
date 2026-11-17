#!/usr/bin/env Rscript
# donor_de.R — GSE36830 polyp-only DE rerun (sensitivity: donor duplication).
# Frozen primary CRSwNP_vs_UTctrl treats 12 case samples as independent:
#   6 CRSwNP uncinate tissue (UT-1..6, GSM902692-697) + 6 CRSwNP nasal polyp
#   (polyp-1..6, GSM902698-703), parallel numbering suggesting the same 6 donors,
#   vs 6 healthy control UT (GSM902680-685).
# This script reruns the SAME pipeline as 03_de.R (run_limma + collapse_sym,
# verbatim logic) with cases restricted to the 6 polyp samples only
# (one sample per putative donor). No frozen file is modified.
# It first replicates the frozen 12-vs-6 contrast and asserts exact agreement
# with the frozen de table, so the polyp-only difference is attributable
# solely to the sample set.
suppressPackageStartupMessages(library(limma))
BASE <- "analysis/step2_bulk"
OUT  <- file.path(BASE, "sensitivity_gse36830_donor")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- verbatim copies from 03_de.R ----
collapse_sym <- function(res, symbols) {
  res$symbol <- symbols[rownames(res)]
  res <- res[!is.na(res$symbol) & res$symbol != "", ]
  res <- res[order(res$P.Value), ]
  res <- res[!duplicated(res$symbol), ]
  res[order(res$P.Value), ]
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
run_contrast <- function(sp, case_ids, ctrl_ids) {
  mat <- sp$expr; sym <- sp$symbols
  ids <- c(case_ids, ctrl_ids)
  grp <- ifelse(colnames(mat) %in% case_ids, "case", "control")
  names(grp) <- colnames(mat); grp <- grp[ids]
  res <- run_limma(mat[, ids], grp)
  res <- collapse_sym(res, sym)
  data.frame(symbol = res$symbol, logFC = res$logFC, AveExpr = res$AveExpr,
             t = res$t, p = res$P.Value, FDR = res$adj.P.Val)
}

## ---- inputs ----
sp <- readRDS(file.path(BASE, "prep", "GSE36830.rds"))
stopifnot(sp$scale_verdict == "linear_normalized_log2t", sp$method == "limma")
ph <- read.csv(file.path(BASE, "meta", "GSE36830_GPL570_pheno.csv"), row.names = 1)
stopifnot(all(colnames(sp$expr) %in% rownames(ph)))

ctrl_ids   <- rownames(ph)[ph$title %in% paste0("Control UT-", 1:6)]
ut_ids     <- rownames(ph)[ph$title %in% paste0("CRSwNP UT-", 1:6)]
polyp_ids  <- rownames(ph)[ph$title %in% paste0("CRSwNP polyp-", 1:6)]
stopifnot(length(ctrl_ids) == 6, length(ut_ids) == 6, length(polyp_ids) == 6)
# tissue field cross-check
stopifnot(all(ph[ut_ids, "characteristics_ch1.2"] == "tissue: uncinate tissue (UT)"))
stopifnot(all(ph[polyp_ids, "characteristics_ch1.2"] == "tissue: nasal polyp (NP) tissue"))
stopifnot(all(sp$meta[ctrl_ids, "group"] == "control"),
          all(sp$meta[c(ut_ids, polyp_ids), "group"] == "case"))

inv <- data.frame(
  gsm = c(ctrl_ids, ut_ids, polyp_ids),
  title = ph[c(ctrl_ids, ut_ids, polyp_ids), "title"],
  tissue = ph[c(ctrl_ids, ut_ids, polyp_ids), "characteristics_ch1.2"],
  frozen_primary_role = c(rep("control", 6), rep("case", 12)),
  polyp_only_role = c(rep("control", 6), rep("excluded_UT_of_CRSwNP_donor", 6), rep("case", 6)))
write.csv(inv, file.path(OUT, "gse36830_sample_inventory.csv"), row.names = FALSE)
cat("inventory written: ctrl=6, CRSwNP UT=6, CRSwNP polyp=6\n")

## ---- step 1: replicate frozen contrast (12 case vs 6 control) ----
rep <- run_contrast(sp, c(ut_ids, polyp_ids), ctrl_ids)
frozen <- read.csv(file.path(BASE, "de", "GSE36830_CRSwNP_vs_UTctrl.csv"),
                   colClasses = c(symbol = "character"))
m <- merge(frozen, rep, by = "symbol", suffixes = c("_frozen", "_rep"))
stopifnot(nrow(m) == nrow(frozen), nrow(m) == nrow(rep))
d <- max(abs(m$logFC_frozen - m$logFC_rep), abs(m$t_frozen - m$t_rep),
         abs(m$p_frozen - m$p_rep), na.rm = TRUE)
cat(sprintf("frozen replication: n=%d genes, max abs diff (logFC/t/p) = %.3g\n", nrow(m), d))
stopifnot(d < 1e-8)
cat("REPLICATION OK: pipeline + prep reproduce the frozen DE table exactly\n")

## ---- step 2: polyp-only contrast (6 polyp vs 6 control UT) ----
new <- run_contrast(sp, polyp_ids, ctrl_ids)
write.csv(new, file.path(OUT, "GSE36830_CRSwNP_polypOnly_vs_UTctrl.csv"), row.names = FALSE)
cat(sprintf("polyp-only DE: n=%d genes, %d FDR<0.05 (frozen 12v6: %d FDR<0.05)\n",
            nrow(new), sum(new$FDR < 0.05, na.rm = TRUE),
            sum(frozen$FDR < 0.05, na.rm = TRUE)))
cat("wrote GSE36830_CRSwNP_polypOnly_vs_UTctrl.csv\n")
