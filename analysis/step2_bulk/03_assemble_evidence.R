#!/usr/bin/env Rscript
# 03_assemble_evidence.R — protein-coding evidence matrices (genes x cohorts) from primary contrasts.
# Reads de/<GSE>_DE.csv, restricts to protein-coding symbols via org.Hs.eg.db,
# dedups duplicate symbols (keep min-FDR row), writes evidence_matrix_logFC.csv / _FDR.csv.
suppressPackageStartupMessages({ library(org.Hs.eg.db) })
OUT <- "analysis/step2_bulk"
fs <- list.files(file.path(OUT, "de"), pattern = "_DE\\.csv$", full.names = TRUE)
cat("cohort DE files:", length(fs), "\n")
if (!length(fs)) quit(status = 1)
tab <- do.call(rbind, lapply(fs, function(f) {
  g <- sub("_DE\\.csv$", "", basename(f))
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$gse <- g
  d
}))
tab$gene <- toupper(trimws(tab$gene))
tab <- tab[!is.na(tab$gene) & tab$gene != "", ]
# dedup: keep row with min FDR per (gse, gene)
tab <- tab[order(tab$gse, tab$FDR), ]
tab <- tab[!duplicated(paste(tab$gse, tab$gene)), ]
cat(sprintf("rows: %d, cohorts: %d, genes: %d\n", nrow(tab), length(unique(tab$gse)), length(unique(tab$gene))))
# protein-coding filter
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
