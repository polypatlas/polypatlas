#!/usr/bin/env Rscript
# F4 step A: per-cohort Spearman top-30 positive partners, with checkpoint files
# (survives VM restarts). Checkpoints: analysis/step9_features/checkpoints/f4_pairs_<GSE>.rds
suppressMessages(library(data.table))

univ <- fread("analysis/step2_bulk/evidence_matrix_logFC.csv", select=1)[[1]]
cohorts <- c("GSE136825","GSE179265","GSE189690","GSE198950","GSE213809",
             "GSE267331","GSE313186","GSE338426","GSE36830","GSE69093",
             "GSE70898_70900")
cpdir <- "analysis/step9_features/checkpoints"
dir.create(cpdir, showWarnings=FALSE)
topN <- 30L
bsize <- 1500L

used <- data.table(cohort=character(), n_samples=integer(), n_genes=integer())
for (gse in cohorts) {
  cp <- file.path(cpdir, paste0("f4_pairs_", gse, ".rds"))
  if (file.exists(cp)) { cat(gse, "checkpoint exists, skipping\n"); next }
  t0 <- Sys.time()
  x <- readRDS(file.path("analysis/step2_bulk/prep", paste0(gse, ".rds")))
  M <- x$expr; sym <- x$symbols
  keep <- !is.na(sym) & sym != ""
  M <- M[keep, , drop=FALSE]; sym <- sym[keep]
  vr <- apply(M, 1, var); ord <- order(sym, -vr)
  M <- M[ord, , drop=FALSE]; sym <- sym[ord]
  M <- M[!duplicated(sym), , drop=FALSE]; sym <- sym[!duplicated(sym)]
  rownames(M) <- sym
  M <- M[intersect(sym, univ), , drop=FALSE]
  M <- t(M)
  zv <- apply(M, 2, function(v) var(v, na.rm=TRUE) == 0 | all(is.na(v)))
  M <- M[, !zv, drop=FALSE]
  ng <- ncol(M); ns <- nrow(M)
  cat(gse, "samples=", ns, "genes=", ng, " prep:", round(difftime(Sys.time(), t0, units="mins"),1), "min\n", sep="")
  used <- rbind(used, data.table(cohort=gse, n_samples=ns, n_genes=ng))
  R <- apply(M, 2, rank, na.last="keep")
  genes <- colnames(R); nb <- ceiling(ng/bsize)
  use_pairwise <- anyNA(R)
  cat(" ", gse, "anyNA in ranks:", use_pairwise, "\n")
  blk_list <- list()
  for (b in seq_len(nb)) {
    bt <- Sys.time()
    idx <- ((b-1)*bsize+1):min(b*bsize, ng)
    C <- if (use_pairwise) cor(R[, idx, drop=FALSE], R, use="pairwise.complete.obs")
         else cor(R[, idx, drop=FALSE], R)   # fast BLAS path: no NAs present
    C[cbind(seq_len(nrow(C)), idx)] <- -Inf
    C[is.na(C) | C <= 0] <- -Inf
    haspos <- rowSums(C > -Inf) > 0
    Cb <- C[haspos, , drop=FALSE]; gb <- genes[idx][haspos]
    rows <- vector("list", nrow(Cb))
    for (j in seq_len(nrow(Cb))) {
      v <- Cb[j, ]; o <- order(v, decreasing=TRUE)
      o <- o[seq_len(min(topN, sum(v > -Inf)))]
      rows[[j]] <- list(gene=gb[j], partner=genes[o], r=v[o])
    }
    bb <- rbindlist(lapply(rows, function(x) data.table(gene=x$gene, partner=x$partner, r=x$r)))
    bb <- bb[is.finite(r) & r > 0]; bb[, cohort := gse]
    blk_list[[b]] <- bb
    cat(" ", gse, "block", b, "/", nb, "pairs:", nrow(bb),
        "(", round(difftime(Sys.time(), bt, units="secs"),0), "s)\n", sep="")
    rm(C, Cb, rows, bb); gc()
  }
  dt <- rbindlist(blk_list)
  saveRDS(dt, cp)
  cat(gse, "CHECKPOINTED pairs:", nrow(dt), " total:",
      round(difftime(Sys.time(), t0, units="mins"),1), "min\n")
  rm(x, M, R, dt, blk_list); gc()
}
fwrite(used, file.path(cpdir, "f4_coexpr_cohorts.csv"))
cat("done step A\n")
