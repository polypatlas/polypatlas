#!/usr/bin/env Rscript
# 04_controls.R — positive-control panel behavior across cohorts.
# Panel (type-2 / barrier genes with strong prior): CLC CCL26 CCL13 IL33 TSLP POSTN SERPINB2 ALOX15
# Writes controls_panel.csv: rows = genes, cols = <GSE>_logFC, <GSE>_FDR; plus incoherence flags.
OUT <- "analysis/step2_bulk"
em <- readRDS(file.path(OUT, "evidence_matrices.rds"))
panel <- c("CLC","CCL26","CCL13","IL33","TSLP","POSTN","SERPINB2","ALOX15")
have <- panel[panel %in% rownames(em$logFC)]
cat("panel genes present in matrix:", paste(have, collapse = ", "), "\n")
miss <- setdiff(panel, have)
if (length(miss)) cat("panel genes ABSENT from all matrices:", paste(miss, collapse = ", "), "\n")
lfc <- em$logFC[have, , drop = FALSE]
fdr <- em$FDR[have, , drop = FALSE]
out <- data.frame(gene = have, stringsAsFactors = FALSE)
for (g in colnames(lfc)) {
  out[[paste0(g, "_logFC")]] <- round(lfc[, g], 3)
  out[[paste0(g, "_FDR")]] <- signif(fdr[, g], 3)
}
write.csv(out, file.path(OUT, "controls_panel.csv"), row.names = FALSE)
# incoherence check: for each cohort, count panel genes with significant DE in the
# expected direction. Expected (nasal polyp vs control, from prior work):
# up: CLC CCL26 CCL13 IL33 TSLP POSTN SERPINB2 ALOX15 (all canonical up in CRSwNP)
exp_up <- panel
flags <- data.frame(gse = colnames(lfc), n_sig_up = NA_integer_, n_measured = NA_integer_,
                    flag = "", stringsAsFactors = FALSE)
for (j in seq_along(colnames(lfc))) {
  g <- colnames(lfc)[j]
  m <- !is.na(lfc[, g])
  flags$n_measured[j] <- sum(m)
  flags$n_sig_up[j] <- sum(m & lfc[, g] > 0 & fdr[, g] < 0.05, na.rm = TRUE)
  if (sum(m) >= 5 && sum(m & lfc[, g] > 0 & fdr[, g] < 0.05, na.rm = TRUE) <= 1)
    flags$flag[j] <- "INVESTIGATE: panel incoherent"
}
write.csv(flags, file.path(OUT, "controls_panel_flags.csv"), row.names = FALSE)
print(out)
print(flags)
cat("DONE\n")
