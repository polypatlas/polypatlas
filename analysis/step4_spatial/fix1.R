suppressMessages({ library(GEOquery); library(openxlsx); library(limma) })
dir <- "analysis/step4_spatial"
lg <- file(file.path(dir, "log/fix1.log"), open="wt")
logm <- function(...) { m <- paste0(...); cat(m, "\n", file=lg, append=TRUE); cat(m, "\n") }

key_of <- function(x) sub("^([A-Za-z]+_[0-9]+)_([0-9]+)_([A-Za-z45]+)$", "\\1:\\3:ROI\\2", x)

g <- getGEO(filename=file.path(dir, "meta/GSE235714_series_matrix.txt.gz"))
p <- pData(phenoData(g))
meta <- data.frame(GSM=rownames(p), title=p$title, source=p$source_name_ch1, tissue_ch=p$characteristics_ch1, stringsAsFactors=FALSE)
mm <- regmatches(meta$title, regexec("^(CRSwNP|CTRL)_([0-9]+), (CD45|PanCK)\\+ ROI([0-9]+)$", meta$title))
meta$donor <- sapply(mm, function(x) if (length(x)==5) paste0(x[2],"_",x[3]) else NA)
meta$group  <- sapply(mm, function(x) if (length(x)==5) x[2] else NA)
meta$segment<- sapply(mm, function(x) if (length(x)==5) x[4] else NA)
meta$key <- paste0(meta$donor, ":", meta$segment, ":ROI", sapply(mm, function(x) if (length(x)==5) x[5] else NA))

seg <- read.xlsx(file.path(dir, "meta/GSE235714_CRS_Q3_Normalization.xlsx"), sheet="SegmentProperties")
seg$key <- key_of(seg$SampleName)
tc <- read.xlsx(file.path(dir, "meta/GSE235714_CRS_Q3_Normalization.xlsx"), sheet="TargetCountMatrix")
tp <- read.xlsx(file.path(dir, "meta/GSE235714_CRS_Q3_Normalization.xlsx"), sheet="TargetProperties")

# merge QC by tc-style key
tc_keys <- key_of(colnames(tc)[-1])
stopifnot(all(tc_keys %in% seg$key))
qc <- seg[match(tc_keys, seg$key), c("key","SampleName","QCFlags","RawReads","DeduplicatedReads","AOINucleiCount","Group",
  "LOQ.(Cancer.Transcriptome.Atlas)","NormalizationFactor")]
m2 <- meta[match(qc$key, meta$key), ]
stopifnot(!any(is.na(m2$key)))
meta_full <- cbind(m2[, c("GSM","title","donor","group","segment","source","tissue_ch")], qc)
meta_full <- meta_full[order(meta_full$donor, meta_full$segment, meta_full$key), ]
write.csv(meta_full, file.path(dir, "roi_metadata.csv"), row.names=FALSE)
logm("roi_metadata: ", nrow(meta_full), " rows")
lowsa <- meta_full[meta_full$QCFlags=="Low Surface Area", ]
logm("Low Surface Area ROI: GSM=", lowsa$GSM, " key=", lowsa$key, " group=", lowsa$Group,
     " RawReads=", lowsa$RawReads, " Dedup=", lowsa$DeduplicatedReads, " Nuclei=", lowsa$AOINucleiCount)

# expression
genes <- tc$TargetName
is_spike <- tp$AnalyteType != "RNA" | grepl("NEG", genes, ignore.case=TRUE)
E <- as.matrix(tc[, -1][!is_spike, ]); rownames(E) <- genes[!is_spike]
E <- E[, match(meta_full$key, tc_keys)]
colnames(E) <- meta_full$key
write.csv(data.frame(gene=rownames(E), E, check.names=FALSE), file.path(dir, "roi_expression.csv"), row.names=FALSE)
logm("roi_expression: ", nrow(E), " x ", ncol(E))

# donor-level means
mk <- paste0(meta_full$donor, ":", meta_full$segment)
uniq <- unique(mk)
Ed <- sapply(uniq, function(k) rowMeans(E[, mk==k, drop=FALSE]))
colnames(Ed) <- uniq
dinf <- data.frame(key=uniq, donor=sub(":.*","",uniq), segment=sub(".*:","",uniq), stringsAsFactors=FALSE)
dinf$group <- ifelse(grepl("^CRSwNP", dinf$donor), "CRSwNP", "CTRL")
logm("donor-level matrix: ", nrow(Ed), " x ", ncol(Ed), " (6 CRSwNP + 4 CTRL donors per segment)")

res <- list()
for (sg in c("CD45","PanCK")) {
  cols <- which(dinf$segment==sg)
  Y <- log2(Ed[, cols] + 1)
  grp <- factor(ifelse(dinf$group[cols]=="CRSwNP","CRSwNP","CTRL"), levels=c("CTRL","CRSwNP"))
  des <- model.matrix(~grp)
  fit <- eBayes(lmFit(Y, des))
  tt <- topTable(fit, coef="grpCRSwNP", number=nrow(Y), sort.by="none")
  out <- data.frame(gene=rownames(Y), logFC=tt$logFC, AveExpr=tt$AveExpr, t=tt$t, p=tt$P.Value, FDR=tt$adj.P.Val, stringsAsFactors=FALSE)
  write.csv(out, file.path(dir, "de", paste0("GSE235714_", sg, "_CRSwNP_vs_CTRL_donorlevel.csv")), row.names=FALSE)
  res[[sg]] <- out
  logm("donor-level DE ", sg, " (6v4): FDR<0.05 = ", sum(out$FDR<0.05, na.rm=TRUE))
}
saveRDS(res, file.path(dir, "de_donorlevel.rds"))

# ---------- alignment ----------
lf <- read.csv("analysis/step2_bulk/evidence_matrix_logFC.csv", check.names=FALSE)
ff <- read.csv("analysis/step2_bulk/evidence_matrix_FDR.csv", check.names=FALSE)
colnames(lf)[1] <- colnames(ff)[1] <- "gene"
bulk_cohorts <- colnames(lf)[-1]
sc <- read.csv("analysis/step3_scrna/celltype_evidence_matrix.csv", check.names=FALSE)
colnames(sc)[1] <- "gene"
scL <- colnames(sc)[grepl("logFC", colnames(sc))]
immune_ct <- c("CD4_T","CD8_T","Treg","NK","NKT","ILC2","B","plasma","mast","macrophage","monocyte","cDC","pDC","neutrophil","eosinophil")
epi_ct <- c("basal","ciliated","secretory","goblet","glandular","ionocyte","tuft","fibroblast","endothelial","pericyte")
pick <- function(pat, cols) cols[grepl(paste0("GSE235711.*(", paste(pat, collapse="|"), ")"), cols)]
lci <- pick(immune_ct, scL); lce <- pick(epi_ct, scL)
logm("scRNA matched cols: immune ", length(lci), " | epithelial ", length(lce))
mean_of <- function(cols) suppressWarnings(rowMeans(as.matrix(sc[match(res$CD45$gene, sc$gene), cols, drop=FALSE]), na.rm=TRUE))

align <- data.frame(gene=res$CD45$gene, stringsAsFactors=FALSE)
align$sp_CD45_logFC <- res$CD45$logFC; align$sp_CD45_FDR <- res$CD45$FDR
align$sp_PanCK_logFC <- res$PanCK$logFC; align$sp_PanCK_FDR <- res$PanCK$FDR
lf2 <- lf[match(align$gene, lf$gene), ]; ff2 <- ff[match(align$gene, ff$gene), ]
for (cc in bulk_cohorts) { align[[paste0("bulk_", cc, "_logFC")]] <- lf2[[cc]]; align[[paste0("bulk_", cc, "_FDR")]] <- ff2[[cc]] }
align$sc_immune_meanlogFC <- mean_of(lci); align$sc_epi_meanlogFC <- mean_of(lce)
write.csv(align, file.path(dir, "alignment_table.csv"), row.names=FALSE)
logm("alignment_table: ", nrow(align), " x ", ncol(align))

conc <- data.frame()
for (sg in c("CD45","PanCK")) {
  sL <- paste0("sp_", sg, "_logFC"); sF <- paste0("sp_", sg, "_FDR")
  sig <- which(align[[sF]] < 0.05)
  for (cc in bulk_cohorts) {
    bL <- paste0("bulk_", cc, "_logFC"); bF <- paste0("bulk_", cc, "_FDR")
    both <- sig[align[[bF]][sig] < 0.05]
    rate <- if (length(both) >= 5) mean(sign(align[[sL]][both]) == sign(align[[bL]][both]), na.rm=TRUE) else NA
    conc <- rbind(conc, data.frame(segment=sg, cohort=cc, n_spatial_sig=length(sig),
      n_overlap_sig=length(both), sign_agreement=round(rate,3), stringsAsFactors=FALSE))
  }
  scm <- paste0("sc_", if (sg=="CD45") "immune" else "epi", "_meanlogFC")
  both <- sig[!is.na(align[[scm]][sig])]
  rate <- if (length(both) >= 5) mean(sign(align[[sL]][both]) == sign(align[[scm]][both]), na.rm=TRUE) else NA
  conc <- rbind(conc, data.frame(segment=sg, cohort="scRNA_GSE235711_mean", n_spatial_sig=length(sig),
    n_overlap_sig=length(both), sign_agreement=round(rate,3), stringsAsFactors=FALSE))
}
write.csv(conc, file.path(dir, "concordance_summary.csv"), row.names=FALSE)

markers <- c("CLC","CCL26","CCL13","IL33","TSLP","POSTN","SERPINB2","ALOX15")
logm("markers in CTA panel: ", paste(markers[markers %in% align$gene], collapse=", "), " (", sum(markers %in% align$gene), "/8)")
mk <- align[align$gene %in% markers, c("gene","sp_CD45_logFC","sp_CD45_FDR","sp_PanCK_logFC","sp_PanCK_FDR")]
print(mk)
close(lg); cat("FIX DONE\n")
