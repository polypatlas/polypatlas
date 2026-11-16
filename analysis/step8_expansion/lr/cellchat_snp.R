# CellChat for CRSsNP only (fresh process to avoid OOM). Compares against saved control object.
suppressMessages({library(CellChat); library(Matrix)})
wd <- "analysis/step8_expansion/lr"
setwd(wd)
counts <- readMM("counts_full.mtx")
meta <- read.csv("meta_full.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$cell_id
colnames(counts) <- meta$cell_id
genes <- readLines("genes_full.txt")
rownames(counts) <- genes
m <- meta[meta$group == "CRSsNP", ]
ct <- counts[, rownames(m)]
cc <- createCellChat(object = ct, meta = m, group.by = "celltype")
cc@DB <- CellChatDB.human
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1)
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
saveRDS(cc, "cellchat_CRSsNP.rds")
cat("CRSsNP interactions: pathways:", length(cc@netP$pathways), "\n")
rm(ct, counts); gc()
# comparison vs control
ctrl <- readRDS("cellchat_control.rds")
dfc <- subsetCommunication(ctrl); dfs <- subsetCommunication(cc)
key <- c("source","target","ligand","receptor")
m2 <- merge(dfc[, c(key, "prob")], dfs[, c(key, "prob")], by = key, suffixes = c("_ctrl","_sn"), all = TRUE)
m2$prob_ctrl[is.na(m2$prob_ctrl)] <- 0; m2$prob_sn[is.na(m2$prob_sn)] <- 0
m2$delta <- m2$prob_sn - m2$prob_ctrl
m2 <- m2[order(-abs(m2$delta), -m2$prob_sn), ]
write.csv(m2, "lr_interactions_CRSsNP_vs_control.csv", row.names = FALSE)
cat("CRSsNP_vs_ctrl interactions:", nrow(m2), "\n")
