# CellChat ligand-receptor analysis, GSE235711 (CRSwNP flagship scRNA)
# Input: counts_full.mtx (CellChatDB genes x annotated cells), meta_full.csv, genes_full.txt
# Raw counts rebuilt from per-donor 10x tars; annotated celltypes/groups from the step-3 h5ad.
# Disease vs control aggregated; CRSsNP vs control as secondary.

suppressMessages({library(CellChat); library(Matrix); library(dplyr)})

wd <- "analysis/step8_expansion/lr"
setwd(wd)
counts <- readMM("counts_full.mtx")
genes <- readLines("genes_full.txt")
rownames(counts) <- genes
meta <- read.csv("meta_full.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$cell_id
colnames(counts) <- meta$cell_id
counts <- counts[, rownames(meta)]
cat("counts:", nrow(counts), "x", ncol(counts), "\n")
print(as.data.frame(table(meta$group, meta$celltype))[1:6, ])

mk <- function(cond) {
  m <- meta[meta$group == cond, ]
  ct <- counts[, rownames(m)]
  cc <- createCellChat(object = ct, meta = m, group.by = "celltype")
  cc <- setIdent(cc, ident.use = "celltype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
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
  saveRDS(cc, file = paste0("cellchat_", cond, ".rds"))
  cat("interactions:", nrow(cc@net$df.net), "pathways:", length(cc@netP$pathways), "\n")
}

# comparison: CRSwNP vs control (primary)
ctrl <- readRDS("cellchat_control.rds"); dis <- readRDS("cellchat_CRSwNP.rds")
merged <- mergeCellChat(list(CTRL = ctrl, CRSwNP = dis), add.names = c("CTRL","CRSwNP"))
gg1 <- compareInteractions(merged, show.legend = FALSE, group = c(1,2))
ggsave("cellchat_compare_interactions.png", gg1, width = 7, height = 5, dpi = 300)

dfc <- ctrl@net$df.net; dfd <- dis@net$df.net
key <- c("source","target","ligand","receptor")
m <- merge(dfc[, c(key, "prob")], dfd[, c(key, "prob")], by = key, suffixes = c("_ctrl","_dis"), all = TRUE)
m$prob_ctrl[is.na(m$prob_ctrl)] <- 0; m$prob_dis[is.na(m$prob_dis)] <- 0
m$delta <- m$prob_dis - m$prob_ctrl
m <- m[order(-abs(m$delta)), -m$prob_dis), ]
write.csv(m, "lr_interactions_CRSwNP_vs_control.csv", row.names = FALSE)

# secondary: CRSsNP vs control
snp <- readRDS("cellchat_CRSsNP.rds")
m2 <- merge(dfc[, c(key, "prob")], snp@net$df.net[, c(key, "prob")], by = key, suffixes = c("_ctrl","_sn"), all = TRUE)
m2$prob_ctrl[is.na(m2$prob_ctrl)] <- 0; m2$prob_sn[is.na(m2$prob_sn)] <- 0
m2$delta <- m2$prob_sn - m2$prob_ctrl
m2 <- m2[order(-abs(m2$delta)), -m2$prob_sn), ]
write.csv(m2, "lr_interactions_CRSsNP_vs_control.csv", row.names = FALSE)

cat("CRSwNP_vs_ctrl interactions:", nrow(m), "| CRSsNP_vs_ctrl:", nrow(m2), "\n")
cat("top10 CRSwNP gainers:\n"); print(head(m[m$delta > 0, ], 10))
cat("top10 CRSwNP losers:\n"); print(head(m[m$delta < 0, ], 10))
