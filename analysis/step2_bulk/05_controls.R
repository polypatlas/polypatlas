#!/usr/bin/env Rscript
# 05_controls.R — positive-control panel behavior across cohorts (primary contrasts).
# Panel: CLC CCL26 CCL13 IL33 TSLP POSTN SERPINB2 ALOX15 (canonical type-2 up in CRSwNP).
# Expected direction depends on contrast type:
#   disease_vs_ctrl -> UP ; treatment_reversal (dupilumab post vs pre) -> DOWN ;
#   subtype (no directional prior) -> no flag.
OUT <- "analysis/step2_bulk"
em <- readRDS(file.path(OUT, "evidence_matrices.rds"))
panel <- c("CLC", "CCL26", "CCL13", "IL33", "TSLP", "POSTN", "SERPINB2", "ALOX15")
have <- panel[panel %in% rownames(em$logFC)]
cat("panel genes present:", paste(have, collapse = ", "), "\n")
miss <- setdiff(panel, have)
if (length(miss)) cat("panel genes ABSENT:", paste(miss, collapse = ", "), "\n")
# contrast type per cohort
ctype <- c(
  GSE338426 = "disease_vs_ctrl", GSE198950 = "disease_vs_ctrl",
  GSE136825 = "disease_vs_ctrl", GSE36830 = "disease_vs_ctrl",
  GSE194282 = "disease_vs_ctrl", GSE72713 = "disease_vs_ctrl",
  GSE179265 = "disease_vs_ctrl", GSE189690 = "subtype",
  GSE267331 = "treatment_reversal", GSE69093 = "disease_vs_ctrl",
  GSE313186 = "exposure_vs_ctrl", GSE70898_70900 = "disease_vs_ctrl",
  GSE174681 = "cell_subset", GSE150846 = "cell_subset",
  GSE330761 = "subtype", GSE213809 = "cell_subset"
)
lfc <- em$logFC[have, , drop = FALSE]
fdr <- em$FDR[have, , drop = FALSE]
out <- data.frame(gene = have, stringsAsFactors = FALSE)
for (g in colnames(lfc)) {
  out[[paste0(g, "_logFC")]] <- round(lfc[, g], 3)
  out[[paste0(g, "_FDR")]] <- signif(fdr[, g], 3)
}
write.csv(out, file.path(OUT, "controls_panel.csv"), row.names = FALSE)
flags <- data.frame(gse = colnames(lfc), contrast_type = ctype[colnames(lfc)],
                    n_sig_expected = NA_integer_, n_measured = NA_integer_,
                    flag = "", stringsAsFactors = FALSE)
for (j in seq_along(colnames(lfc))) {
  g <- colnames(lfc)[j]
  m <- !is.na(lfc[, g])
  flags$n_measured[j] <- sum(m)
  ct <- ctype[[g]]
  if (ct %in% c("disease_vs_ctrl")) {
    s <- sum(m & lfc[, g] > 0 & fdr[, g] < 0.05, na.rm = TRUE)
  } else if (ct == "treatment_reversal") {
    s <- sum(m & lfc[, g] < 0 & fdr[, g] < 0.05, na.rm = TRUE)
  } else if (ct == "exposure_vs_ctrl") {
    # weaker prior: flag only if panel is significantly REVERSED
    s <- sum(m & lfc[, g] > 0 & fdr[, g] < 0.05, na.rm = TRUE)
    rev <- sum(m & lfc[, g] < 0 & fdr[, g] < 0.05, na.rm = TRUE)
    if (!is.na(rev) && rev >= 3) flags$flag[j] <- "INVESTIGATE: panel reversed in exposure contrast"
    flags$n_sig_expected[j] <- s; next
  } else s <- NA_integer_
  flags$n_sig_expected[j] <- s
  if (!is.na(s) && sum(m) >= 5 && s <= 1)
    flags$flag[j] <- "INVESTIGATE: panel incoherent with contrast expectation"
}
write.csv(flags, file.path(OUT, "controls_panel_flags.csv"), row.names = FALSE)
print(flags)
cat("DONE\n")
