# CellChat ligand-receptor analysis, GSE235711 (CRSwNP flagship scRNA)
# NORMALIZED RERUN (2026-11-16) of ../cellchat_run.R.
# Single methodological change vs the original: library-size normalization
# (CellChat::normalizeData, scale.factor=10000, log1p) inserted after
# subsetData() and before identifyOverExpressedGenes().
#
# Implementation note (CellChat 1.6.1): normalizeData() in this version takes
# a data matrix (data.raw, scale.factor=10000, do.log=TRUE), NOT a CellChat
# object (no normalizeData.CellChat method exists). The official v1 workflow
# applies it to the expression matrix; here it is applied to the exact matrix
# every downstream step reads, i.e. cc@data.signaling as set by subsetData()
# (verified in the installed package source: identifyOverExpressedGenes,
# identifyOverExpressedInteractions and computeCommunProb reference
# @data.signaling only). The input matrix was already restricted to
# CellChatDB genes (1,075 genes), so this is the post-subsetData equivalent
# of the official pre-createCellChat normalization.
#
# Repairs vs the saved ../cellchat_run.R (documented, non-methodological):
#  - lines 54/62 of the original contain a syntax error
#    (m[order(-abs(m$delta)), -m$prob_dis), ]) and cannot parse; the ordering
#    expression is restored to the form used by ../cellchat_compare.R and
#    ../cellchat_snp.R, which produced the frozen CSVs.
#  - comparison tables are built from subsetCommunication() (as in
#    cellchat_compare.R / cellchat_snp.R) instead of @net$df.net, which is
#    NULL in the saved 1.6.1 objects; subsetCommunication row counts equal
#    sum(@net$count) in the frozen objects.
# Inputs are read from the original directory by absolute path; all outputs
# are written to this rerun directory with a _norm suffix. Frozen files are
# never overwritten.

suppressMessages({library(CellChat); library(Matrix); library(dplyr)})

cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")
cat("rerun started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
t0 <- Sys.time()

orig <- "analysis/step8_expansion/lr"
outd <- "analysis/step8_expansion/lr/cellchat_rerun_norm"
setwd(outd)
counts <- readMM(file.path(orig, "counts_full.mtx"))
genes <- readLines(file.path(orig, "genes_full.txt"))
rownames(counts) <- genes
meta <- read.csv(file.path(orig, "meta_full.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$cell_id
colnames(counts) <- meta$cell_id
counts <- counts[, rownames(meta)]
cat("counts:", nrow(counts), "x", ncol(counts), "\n")
print(table(meta$group))

mk <- function(cond) {
  m <- meta[meta$group == cond, ]
  ct <- counts[, rownames(m)]
  cc <- createCellChat(object = ct, meta = m, group.by = "celltype")
  cc <- setIdent(cc, ident.use = "celltype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  ## === THE FIX: normalize (was absent in the original run) ===
  cc@data.signaling <- normalizeData(cc@data.signaling)
  ## ============================================================
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

for (cond in c("control", "CRSwNP", "CRSsNP")) {
  cat("=== building", cond, "===\n")
  cc <- mk(cond)
  saveRDS(cc, file = paste0("cellchat_", cond, "_norm.rds"))
  cat("interactions:", nrow(subsetCommunication(cc)),
      "| net count sum:", sum(cc@net$count),
      "| pathways:", length(cc@netP$pathways), "\n")
  rm(cc); gc()
}

# comparison: CRSwNP vs control (primary)
ctrl <- readRDS("cellchat_control_norm.rds"); dis <- readRDS("cellchat_CRSwNP_norm.rds")
merged <- mergeCellChat(list(CTRL = ctrl, CRSwNP = dis), add.names = c("CTRL","CRSwNP"))
gg1 <- compareInteractions(merged, show.legend = FALSE, group = c(1,2))
ggsave("cellchat_compare_interactions_norm.png", gg1, width = 7, height = 5, dpi = 300)

dfc <- subsetCommunication(ctrl); dfd <- subsetCommunication(dis)
key <- c("source","target","ligand","receptor")
m <- merge(dfc[, c(key, "prob")], dfd[, c(key, "prob")], by = key, suffixes = c("_ctrl","_dis"), all = TRUE)
m$prob_ctrl[is.na(m$prob_ctrl)] <- 0; m$prob_dis[is.na(m$prob_dis)] <- 0
m$delta <- m$prob_dis - m$prob_ctrl
m <- m[order(-abs(m$delta), -m$prob_dis), ]
write.csv(m, "lr_interactions_CRSwNP_vs_control_norm.csv", row.names = FALSE)

# secondary: CRSsNP vs control (logic merged from ../cellchat_snp.R)
snp <- readRDS("cellchat_CRSsNP_norm.rds")
dfs <- subsetCommunication(snp)
m2 <- merge(dfc[, c(key, "prob")], dfs[, c(key, "prob")], by = key, suffixes = c("_ctrl","_sn"), all = TRUE)
m2$prob_ctrl[is.na(m2$prob_ctrl)] <- 0; m2$prob_sn[is.na(m2$prob_sn)] <- 0
m2$delta <- m2$prob_sn - m2$prob_ctrl
m2 <- m2[order(-abs(m2$delta), -m2$prob_sn), ]
write.csv(m2, "lr_interactions_CRSsNP_vs_control_norm.csv", row.names = FALSE)

cat("CRSwNP_vs_ctrl interactions:", nrow(m), "| CRSsNP_vs_ctrl:", nrow(m2), "\n")
cat("top10 CRSwNP gainers:\n"); print(head(m[m$delta > 0, ], 10))
cat("top10 CRSwNP losers:\n"); print(head(m[m$delta < 0, ], 10))
cat("elapsed mins:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("rerun finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
