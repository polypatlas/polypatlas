# CellChat comparison, GSE235711: CRSwNP vs control (primary). Fresh process, loads saved objects only.
suppressMessages({library(CellChat); library(dplyr)})
wd <- "analysis/step8_expansion/lr"
setwd(wd)
ctrl <- readRDS("cellchat_control.rds"); dis <- readRDS("cellchat_CRSwNP.rds")
# subsetCommunication returns the interaction data.frame (not the object)
dfc <- subsetCommunication(ctrl); dfd <- subsetCommunication(dis)
key <- c("source","target","ligand","receptor")
m <- merge(dfc[, c(key, "prob")], dfd[, c(key, "prob")], by = key, suffixes = c("_ctrl","_dis"), all = TRUE)
m$prob_ctrl[is.na(m$prob_ctrl)] <- 0; m$prob_dis[is.na(m$prob_dis)] <- 0
m$delta <- m$prob_dis - m$prob_ctrl
m <- m[order(-abs(m$delta), -m$prob_dis), ]
write.csv(m, "lr_interactions_CRSwNP_vs_control.csv", row.names = FALSE)
cat("interactions:", nrow(m), "\n")
cat("--- top 10 gained in CRSwNP ---\n"); print(as.data.frame(head(m[m$delta>0, ],10))[, c(key,"prob_ctrl","prob_dis","delta")])
cat("--- top 10 lost in CRSwNP ---\n"); print(as.data.frame(head(m[m$delta<0, ],10))[, c(key,"prob_ctrl","prob_dis","delta")])
cat("--- IL33/TSLP/IL13/CCL26/CCL13 axis ---\n")
ax <- m[m$ligand %in% c("IL33","TSLP","IL13","CCL26","CCL13","POSTN","IL5","CSF2"), ]
print(as.data.frame(ax[order(-abs(ax$delta)), ])[, c(key,"prob_ctrl","prob_dis","delta")]))
