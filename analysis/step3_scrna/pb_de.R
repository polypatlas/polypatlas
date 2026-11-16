#!/usr/bin/env Rscript
# Pseudobulk DE with edgeR QL F-test. Usage:
# Rscript pb_de.R counts.csv meta.csv groupA groupB out.csv
# counts.csv: donor x gene counts (first column = donor id)
# meta.csv: donor,group
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
counts_file <- args[1]; meta_file <- args[2]; gA <- args[3]; gB <- args[4]; out_file <- args[5]
suppressPackageStartupMessages({library(edgeR); library(limma)})
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
counts <- counts[meta$donor, ]
keep <- meta$group %in% c(gA, gB)
counts <- counts[keep, ]; meta <- meta[keep, ]
if (nrow(meta) < 6) { write.csv(data.frame(gene=character(), logFC=numeric(), PValue=numeric(), FDR=numeric()), out_file, row.names = FALSE); quit(save = "no") }
y <- DGEList(counts = t(as.matrix(counts)), group = factor(meta$group, levels = c(gB, gA)))
keep2 <- filterByExpr(y)
y <- y[keep2, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0 + y$samples$group); colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
con <- makeContrasts(contrasts = paste0(gA, "-", gB), levels = design)
qlf <- glmQLFTest(fit, contrast = con)
tt <- topTags(qlf, n = Inf)$table
out <- data.frame(gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM,
                  PValue = tt$PValue, FDR = tt$FDR, row.names = NULL)
write.csv(out, out_file, row.names = FALSE)
