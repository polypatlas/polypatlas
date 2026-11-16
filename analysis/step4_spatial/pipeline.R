suppressMessages({
  library(GEOquery); library(openxlsx); library(limma)
})
dir <- "analysis/step4_spatial"
lg <- file(file.path(dir, "log/pipeline.log"), open="wt")
logm <- function(...) { m <- paste0(...); cat(m, "\n", file=lg, append=TRUE); cat(m, "\n") }

# ---------- 1. ROI metadata ----------
g <- getGEO(filename=file.path(dir, "meta/GSE235714_series_matrix.txt.gz"))
p <- pData(phenoData(g))
meta <- data.frame(
  GSM = rownames(p),
  title = p$title,
  source = p$source_name_ch1,
  tissue_ch = p$characteristics_ch1,
  stringsAsFactors = FALSE
)
# parse "CRSwNP_4, CD45+ ROI1"
mm <- regmatches(meta$title, regexec("^(CRSwNP|CTRL)_([0-9]+), (CD45|PanCK)\\+ ROI([0-9]+)$", meta$title))
ok <- sapply(mm, length) == 5
logm("titles parsed: ", sum(ok), "/", nrow(meta))
meta$donor <- sapply(mm, function(x) if (length(x)==5) paste0(x[2], "_", x[3]) else NA)
meta$group  <- sapply(mm, function(x) if (length(x)==5) x[2] else NA)
meta$segment<- sapply(mm, function(x) if (length(x)==5) x[4] else NA)
meta$roi_no <- sapply(mm, function(x) if (length(x)==5) x[5] else NA)

seg <- read.xlsx(file.path(dir, "meta/GSE235714_CRS_Q3_Normalization.xlsx"), sheet="SegmentProperties")
# match TargetCountMatrix columns to series titles by parsing both sides
parse_col <- function(x) sub("^([A-Za-z]+)_([0-9]+)_([0-9]+)_([A-Za-z45]+)$", "\\1_\\2:\\4:ROI\\3", x)
# series titles -> key
meta$key <- paste0(meta$donor, ":", meta$segment, ":ROI", meta$roi_no)
# seg SampleName -> key
sm <- regmatches(seg$SampleName, regexec("^(CRSwNP|CTRL)_([0-9]+), (CD45|PanCK)\\+ ROI([0-9]+)$", seg$SampleName))
seg$key <- sapply(sm, function(x) if (length(x)==5) paste0(x[2], "_", x[3], ":", x[4], ":ROI", x[5]) else NA)
logm("seg keys parsed: ", sum(!is.na(seg$key)), "/", nrow(seg))

tc <- read.xlsx(file.path(dir, "meta/GSE235714_CRS_Q3_Normalization.xlsx"), sheet="TargetCountMatrix")
tp <- read.xlsx(file.path(dir, "meta/GSE235714_CRS_Q3_Normalization.xlsx"), sheet="TargetProperties")
tc_keys <- parse_col(colnames(tc)[-1])
# match tc columns to seg keys: seg$SampleName "CRSwNP_4, CD45+ ROI1" -> "CRSwNP_4:CD45:ROI1"; tc "CRSwNP_4_1_CD45" -> "CRSwNP_4:CD45:ROI1"
segkey_for_tc <- sapply(tc_keys, function(k) {
  mm2 <- regmatches(k, regexec("^([A-Za-z]+_[0-9]+):(CD45|PanCK):ROI([0-9]+)$", k))[[1]]
  if (length(mm2)==4) paste0(mm2[2], ":", mm2[3], ":ROI", mm2[4]) else NA
})
match_ok <- segkey_for_tc %in% seg$key
logm("tc columns matched to seg: ", sum(match_ok), "/", length(tc_keys))

meta2 <- merge(meta[, c("GSM","title","source","tissue_ch","donor","group","segment","roi_no","key")],
               seg[, c("key","QCFlags","RawReads","DeduplicatedReads","AOINucleiCount","Group")],
               by="key", all.x=TRUE)
meta2 <- meta2[order(meta2$donor, meta2$segment, meta2$roi_no), ]
rownames(meta2) <- NULL
# order tc columns to match meta2
ord <- match(meta2$key, segkey_for_tc)
logm("order mapping ok: ", sum(!is.na(ord)), "/", nrow(meta2))
stopifnot(all(!is.na(ord)))
write.csv(meta2, file.path(dir, "roi_metadata.csv"), row.names=FALSE)
logm("roi_metadata rows: ", nrow(meta2), " | donors: ", length(unique(meta2$donor)),
     " | segments: ", paste(table(meta2$segment), collapse=", "))

# ---------- 2. expression matrix ----------
genes <- tc$TargetName
is_spike <- tp$AnalyteType != "RNA" | grepl("NEG", genes, ignore.case=TRUE)
logm("total targets: ", length(genes), " | RNA targets: ", sum(!is_spike), " | non-RNA: ", sum(is_spike))
E <- as.matrix(tc[, -1][!is_spike, ])
rownames(E) <- genes[!is_spike]
E <- E[, ord]
colnames(E) <- meta2$key
expr_df <- data.frame(gene = rownames(E), E, check.names=FALSE)
write.csv(expr_df, file.path(dir, "roi_expression.csv"), row.names=FALSE)
logm("roi_expression: ", nrow(E), " genes x ", ncol(E), " ROIs")

# QC flags summary
logm("QCFlags table: ", paste(names(table(meta2$QCFlags, useNA="ifany")), collapse=", "), " -> ",
     paste(as.integer(table(meta2$QCFlags, useNA="ifany")), collapse=", "))
lowseg <- meta2[meta2$QCFlags=="Low Surface Area", ]
logm("Low Surface Area ROI: ", paste(lowseg$key, collapse="; "))

# ---------- 3. DE per segment ----------
El <- log2(E + 1)
donors <- meta2$donor
grp <- factor(ifelse(meta2$group=="CRSwNP", "CRSwNP", "CTRL"), levels=c("CTRL","CRSwNP"))
des <- model.matrix(~grp)
de_all <- list()
for (sg in c("CD45","PanCK")) {
  idx <- meta2$segment == sg
  y <- El[, idx]
  d2 <- des[idx, , drop=FALSE]; dn <- donors[idx]
  dup <- duplicateCorrelation(y, d2, block=dn)
  fit <- lmFit(y, d2, block=dn, correlation=dup$consensus.correlation)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef="grpCRSwNP", number=nrow(y), sort.by="none")
  out <- data.frame(gene=rownames(y), logFC=tt$logFC, AveExpr=tt$AveExpr, t=tt$t, p=tt$P.Value, FDR=tt$adj.P.Val,
                    cor=dup$consensus.correlation, stringsAsFactors=FALSE)
  fn <- file.path(dir, "de", paste0("GSE235714_", sg, "_CRSwNP_vs_CTRL.csv"))
  write.csv(out, fn, row.names=FALSE)
  de_all[[sg]] <- out
  logm("DE ", sg, ": FDR<0.05 = ", sum(out$FDR < 0.05, na.rm=TRUE), " / ", nrow(out),
       " | donor-consensus-cor = ", round(dup$consensus.correlation, 3))
}
saveRDS(de_all, file.path(dir, "de_results.rds"))

# ---------- 4. Alignment ----------
logm("--- alignment ---")
lf <- read.csv("analysis/step2_bulk/evidence_matrix_logFC.csv", check.names=FALSE)
ff <- read.csv("analysis/step2_bulk/evidence_matrix_FDR.csv", check.names=FALSE)
colnames(lf)[1] <- colnames(ff)[1] <- "gene"
bulk_cohorts <- colnames(lf)[-1]
logm("bulk matrix: ", nrow(lf), " genes x ", length(bulk_cohorts), " cohorts")

sc <- read.csv("analysis/step3_scrna/celltype_evidence_matrix.csv", check.names=FALSE)
colnames(sc)[1] <- "gene"
logm("scRNA matrix: ", nrow(sc), " genes x ", (ncol(sc)-1), " cols")
logm("scRNA col sample: ", paste(colnames(sc)[2:7], collapse=" | "))

immune_ct <- c("CD4_T","CD8_T","Treg","NK","NKT","ILC2","B","plasma","mast","macrophage","monocyte","cDC","pDC","neutrophil","eosinophil")
epi_ct <- c("basal","ciliated","secretory","goblet","glandular","ionocyte","tuft","fibroblast","endothelial","pericyte")

# build alignment table: gene x [spatial CD45/PanCK logFC+FDR] + [bulk cohorts] + [scRNA matched]
matched_sc_cd45 <- grep(paste0("GSE235711.*(", paste(immune_ct, collapse="|"), ")"), colnames(sc), value=TRUE)
matched_sc_panck <- grep(paste0("GSE235711.*(", paste(epi_ct, collapse="|"), ")"), colnames(sc), value=TRUE)
logm("matched scRNA cols: CD45-linked ", length(matched_sc_cd45), " | PanCK-linked ", length(matched_sc_panck))

align <- data.frame(gene = de_all$CD45$gene, stringsAsFactors=FALSE)
align$sp_CD45_logFC <- de_all$CD45$logFC; align$sp_CD45_FDR <- de_all$CD45$FDR
align$sp_PanCK_logFC <- de_all$PanCK$logFC; align$sp_PanCK_FDR <- de_all$PanCK$FDR
# bulk
lf2 <- lf[match(align$gene, lf$gene), ]; ff2 <- ff[match(align$gene, ff$gene), ]
for (cc in bulk_cohorts) { align[[paste0("bulk_", cc, "_logFC")]] <- lf2[[cc]]; align[[paste0("bulk_", cc, "_FDR")]] <- ff2[[cc]] }
# scRNA: average logFC across matched cell-type contrasts (simple mean of available logFC cols)
scL <- colnames(sc)[grepl("logFC", colnames(sc))]
scF <- colnames(sc)[grepl("FDR", colnames(sc))]
mean_of <- function(cols) {
  m <- as.matrix(sc[match(align$gene, sc$gene), cols, drop=FALSE])
  suppressWarnings(rowMeans(m, na.rm=TRUE)); 
}
if (length(matched_sc_cd45)) {
  lc <- intersect(matched_sc_cd45, scL); fc <- sub("logFC","FDR",lc); fc <- intersect(fc, scF)
  align$sc_immune_meanlogFC <- mean_of(lc)
}
if (length(matched_sc_panck)) {
  lc <- intersect(matched_sc_panck, scL); fc <- sub("logFC","FDR",lc); fc <- intersect(fc, scF)
  align$sc_epi_meanlogFC <- mean_of(lc)
}
write.csv(align, file.path(dir, "alignment_table.csv"), row.names=FALSE)
logm("alignment_table: ", nrow(align), " genes x ", ncol(align), " cols")

# concordance: for spatial-significant genes, sign agreement with each bulk cohort (where FDR<0.05 there too) + with scRNA mean
conc <- data.frame()
for (sg in c("CD45","PanCK")) {
  sL <- paste0("sp_", sg, "_logFC"); sF <- paste0("sp_", sg, "_FDR")
  sig <- which(align[[sF]] < 0.05)
  for (cc in bulk_cohorts) {
    bL <- paste0("bulk_", cc, "_logFC"); bF <- paste0("bulk_", cc, "_FDR")
    both <- sig[align[[bF]][sig] < 0.05]
    rate <- if (length(both) >= 3) mean(sign(align[[sL]][both]) == sign(align[[bL]][both]), na.rm=TRUE) else NA
    conc <- rbind(conc, data.frame(segment=sg, cohort=cc, n_spatial_sig=length(sig),
                                  n_overlap_sig=length(both), sign_agreement=rate, stringsAsFactors=FALSE))
  }
  scm <- paste0("sc_", if (sg=="CD45") "immune" else "epi", "_meanlogFC")
  if (scm %in% colnames(align)) {
    both <- sig[!is.na(align[[scm]][sig]) & align[[scm]][sig] != 0]
    rate <- if (length(both) >= 3) mean(sign(align[[sL]][both]) == sign(align[[scm]][both]), na.rm=TRUE) else NA
    conc <- rbind(conc, data.frame(segment=sg, cohort="scRNA_same_study_mean", n_spatial_sig=length(sig),
                                  n_overlap_sig=length(both), sign_agreement=rate, stringsAsFactors=FALSE))
  }
}
write.csv(conc, file.path(dir, "concordance_summary.csv"), row.names=FALSE)
logm("concordance rows: ", nrow(conc))

# marker behavior
markers <- c("CLC","CCL26","CCL13","IL33","TSLP","POSTN","SERPINB2","ALOX15")
mk <- align[align$gene %in% markers, ]
mcols <- c("gene","sp_CD45_logFC","sp_CD45_FDR","sp_PanCK_logFC","sp_PanCK_FDR")
mb <- merge(mk[, mcols], lf[lf$gene %in% markers, c("gene", bulk_cohorts)], by="gene", all.x=TRUE)
write.csv(mb, file.path(dir, "markers_alignment.csv"), row.names=FALSE)
logm("markers present in panel: ", nrow(mk), "/8")

# panel coverage of bulk headline genes (top100 |logFC| per cohort union)
topg <- unique(unlist(lapply(bulk_cohorts, function(cc) {
  o <- order(abs(lf[[cc]]), decreasing=TRUE); lf$gene[o[1:min(100, nrow(lf))]]
})))
logm("bulk top100-union genes: ", length(topg), " | covered by panel: ", sum(topg %in% align$gene),
     " (", round(100*sum(topg %in% align$gene)/length(topg), 1), "%)")
# step3 headline genes presence
scgenes <- c("EPAS1","PTGS2","BIRC3","CXCR6","CRTAM","ANKRD37","POU2F3","TRPM5","FOXI1","S100A8")
logm("step3 headline genes in panel: ", paste(paste0(scgenes, "=", scgenes %in% align$gene), collapse=", "))

close(lg)
cat("PIPELINE DONE\n")
