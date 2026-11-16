#!/usr/bin/env Rscript
# 04_assemble_evidence.R — protein-coding evidence matrices (genes x cohorts), PRIMARY contrasts only.
suppressPackageStartupMessages({ library(org.Hs.eg.db) })
OUT <- "analysis/step2_bulk"
primaries <- c(
  GSE338426 = "GSE338426_CRSwNP_vs_ctrl.csv",
  GSE198950 = "GSE198950_CRSwNP_vs_NonCRS.csv",
  GSE136825 = "GSE136825_polyp_vs_ctrl.csv",
  GSE36830 = "GSE36830_CRSwNP_vs_UTctrl.csv",
  GSE194282 = "GSE194282_NP_vs_ctrl.csv",
  GSE72713 = "GSE72713_CRSwNP_vs_ctrl.csv",
  GSE179265 = "GSE179265_CRSwNP_vs_ctrl.csv",
  GSE189690 = "GSE189690_NERD_vs_nonNERD.csv",
  GSE267331 = "GSE267331_dupilumab_post_vs_pre.csv",
  GSE69093 = "GSE69093_case_vs_ctrl_celltypeAdj.csv",
  GSE313186 = "GSE313186_exposed_vs_ctrl.csv",
  GSE70898_70900 = "GSE70898_70900_polypT_vs_nasalT.csv",
  GSE174681 = "GSE174681_NP_ASC_vs_PBLnaive.csv",
  GSE150846 = "GSE150846_NP_ILC2_vs_blood_ILC2.csv",
  GSE330761 = "GSE330761_CRSwNP_MCTC_vs_CRSsNP_MCTC.csv",
  GSE213809 = "GSE213809_NP_vs_PB_stimulated.csv"
)
fs <- file.path(OUT, "de", primaries)
names(fs) <- names(primaries)
missing <- fs[!file.exists(fs)]
if (length(missing)) { cat("MISSING:", paste(names(missing), collapse = ","), "\n"); quit(status = 1) }
tab <- do.call(rbind, lapply(names(fs), function(g) {
  d <- read.csv(fs[[g]], stringsAsFactors = FALSE)
  d$gse <- g; d
}))
tab$gene <- toupper(trimws(tab$symbol))
tab <- tab[!is.na(tab$gene) & tab$gene != "", ]
tab <- tab[order(tab$gse, tab$FDR), ]
tab <- tab[!duplicated(paste(tab$gse, tab$gene)), ]
cat(sprintf("rows: %d, cohorts: %d, genes: %d\n", nrow(tab), length(unique(tab$gse)), length(unique(tab$gene))))
pc <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(tab$gene), keytype = "SYMBOL",
                            columns = "GENETYPE")
pcgenes <- unique(pc$SYMBOL[pc$GENETYPE == "protein-coding"])
cat("protein-coding symbols in universe:", length(pcgenes), "\n")
tab <- tab[tab$gene %in% pcgenes, ]
mat_lfc <- with(tab, tapply(logFC, list(gene, gse), function(x) x[1]))
mat_fdr <- with(tab, tapply(FDR, list(gene, gse), function(x) x[1]))
write.csv(mat_lfc, file.path(OUT, "evidence_matrix_logFC.csv"))
write.csv(mat_fdr, file.path(OUT, "evidence_matrix_FDR.csv"))
cat(sprintf("matrices: %d genes x %d cohorts\n", nrow(mat_lfc), ncol(mat_lfc)))
cat("non-NA coverage per cohort:\n")
print(colSums(!is.na(mat_lfc)))
saveRDS(list(logFC = mat_lfc, FDR = mat_fdr), file.path(OUT, "evidence_matrices.rds"))
cat("DONE\n")
