#!/usr/bin/env Rscript
# PolypAtlas enrichment layer: ORA + GSEA across 71 contrasts
# Rules per task brief:
#  - symbol->Entrez via AnnotationDbi::select(org.Hs.eg.db); drop NA, dedupe by smallest p
#  - ORA up/down = FDR<0.05 & logFC>0 / <0; skip direction if <10 genes (no padding)
#  - enrichGO BP/CC/MF, enrichPathway, enrichKEGG (BH, pvalueCutoff 0.05, qvalueCutoff 0.2)
#  - GSEA (fgsea) on sign(logFC)*(-log10(p)), p floored at 1e-300; H, C5 BP, C2 REACTOME, C2 KEGG;
#    minSize=15, maxSize=500, nperm=10000; powered contrasts only
# Resumable via checkpoint files in chk_ora/ and chk_gsea/.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(clusterProfiler); library(ReactomePA); library(msigdbr)
  library(org.Hs.eg.db); library(fgsea); library(AnnotationDbi)
})

OUT <- "analysis/step8_expansion/enrichment"
BULK <- "analysis/step2_bulk/de"
SCR  <- "analysis/step3_scrna/pseudobulk_de"
dir.create(file.path(OUT, "chk_ora"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "chk_gsea"), showWarnings = FALSE, recursive = TRUE)

logmsg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console() }

## ---------- 1. load contrasts ----------
bulk_files <- sort(list.files(BULK, pattern = "\\.csv$", full.names = TRUE))
scr_files  <- sort(list.files(SCR,  pattern = "\\.csv$", full.names = TRUE))
contrast_id <- function(f) sub("\\.csv$", "", basename(f))

load_one <- function(f, layer) {
  d <- read.csv(f)
  if (layer == "bulk") {
    d2 <- data.frame(symbol = as.character(d$symbol), logFC = d$logFC, p = d$p, FDR = d$FDR)
  } else {
    d2 <- data.frame(symbol = as.character(d$gene), logFC = d$logFC, p = d$PValue, FDR = d$FDR)
  }
  d2 <- d2[!is.na(d2$symbol) & d2$symbol != "" & !is.na(d2$p) & !is.na(d2$FDR) & !is.na(d2$logFC), ]
  d2
}

## ---------- 2. symbol -> entrez, dedupe by smallest p ----------
logmsg("mapping symbols to Entrez")
alldat <- list()
mapinfo <- list()
bf <- setNames(bulk_files, rep("bulk", length(bulk_files)))
sf <- setNames(scr_files, rep("scrna", length(scr_files)))
all_syms <- unique(unlist(mapply(function(f, ly) load_one(f, ly)$symbol, c(bf, sf), names(c(bf, sf)), SIMPLIFY = FALSE)))
mp <- AnnotationDbi::select(org.Hs.eg.db, keys = all_syms, keytype = "SYMBOL", columns = "ENTREZID")
mp <- mp[!is.na(mp$ENTREZID), ]
mp$ENTREZID <- as.character(mp$ENTREZID)

load_map <- function(f, layer) {
  cid <- contrast_id(f)
  d <- load_one(f, layer)
  m <- merge(d, mp, by.x = "symbol", by.y = "SYMBOL")
  n_unmapped <- nrow(d) - length(unique(d$symbol[d$symbol %in% mp$SYMBOL]))
  # dedupe by smallest p
  m <- m[order(m$p), ]
  m <- m[!duplicated(m$ENTREZID), ]
  m <- m[order(m$p), ]
  e2s <- setNames(m$symbol, m$ENTREZID)
  list(cid = cid, layer = layer, dat = m, e2s = e2s, n_genes = nrow(m),
       n_unmapped = n_unmapped,
       n_up = sum(m$FDR < 0.05 & m$logFC > 0), n_down = sum(m$FDR < 0.05 & m$logFC < 0))
}
allc <- c(lapply(bulk_files, load_map, "bulk"), lapply(scr_files, load_map, "scrna"))
names(allc) <- sapply(allc, `[[`, "cid")
saveRDS(allc, file.path(OUT, "chk_mapped.rds"))
logmsg("mapped", length(allc), "contrasts")

## ---------- 3. ORA ----------
ora_row <- function(res, cid, layer, direction, db) {
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  e2s <- allc[[cid]]$e2s
  sym_of <- function(gid) {
    ids <- strsplit(gid, "/", fixed = TRUE)[[1]]
    paste(e2s[ids], collapse = ";")
  }
  data.frame(
    contrast = cid, layer = layer, direction = direction, database = db,
    term_id = df$ID, term_name = df$Description,
    set_size = as.integer(sub("/.*", "", df$BgRatio)),
    count = df$Count, pvalue = df$pvalue, FDR = df$p.adjust,
    gene_ratio = df$GeneRatio, genes = vapply(df$geneID, sym_of, character(1)),
    stringsAsFactors = FALSE)
}

run_ora_dir <- function(cid, layer, genes_entrez, universe, direction) {
  out <- list()
  safe <- function(expr, db) {
    tryCatch({ r <- expr; ora_row(r, cid, layer, direction, db) },
             error = function(e) { logmsg("ORA FAIL", cid, direction, db, conditionMessage(e)); NULL })
  }
  r <- safe(enrichGO(gene = genes_entrez, universe = universe, OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     readable = FALSE), "GO_BP")
  if (!is.null(r)) out <- c(out, list(r))
  r <- safe(enrichGO(gene = genes_entrez, universe = universe, OrgDb = org.Hs.eg.db,
                     ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     readable = FALSE), "GO_CC")
  if (!is.null(r)) out <- c(out, list(r))
  r <- safe(enrichGO(gene = genes_entrez, universe = universe, OrgDb = org.Hs.eg.db,
                     ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     readable = FALSE), "GO_MF")
  if (!is.null(r)) out <- c(out, list(r))
  r <- safe(enrichPathway(gene = genes_entrez, universe = universe, organism = "human",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                          readable = FALSE), "Reactome")
  if (!is.null(r)) out <- c(out, list(r))
  r <- safe(enrichKEGG(gene = genes_entrez, universe = universe, organism = "hsa",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2), "KEGG")
  if (!is.null(r)) out <- c(out, list(r))
  out
}

ora_status <- list()
for (cid in names(allc)) {
  chk <- file.path(OUT, "chk_ora", paste0(cid, ".rds"))
  if (file.exists(chk)) { ora_status[[cid]] <- readRDS(chk)$status; next }
  a <- allc[[cid]]; m <- a$dat
  univ <- m$ENTREZID
  up_e <- m$ENTREZID[m$FDR < 0.05 & m$logFC > 0]
  dn_e <- m$ENTREZID[m$FDR < 0.05 & m$logFC < 0]
  rows <- list(); st <- c(up = NA, down = NA)
  if (length(up_e) >= 10) { rows <- c(rows, run_ora_dir(cid, a$layer, up_e, univ, "up")); st["up"] <- "done" }
  else st["up"] <- "skipped_insufficient_DE"
  if (length(dn_e) >= 10) { rows <- c(rows, run_ora_dir(cid, a$layer, dn_e, univ, "down")); st["down"] <- "done" }
  else st["down"] <- "skipped_insufficient_DE"
  ora <- if (length(rows)) do.call(rbind, rows) else NULL
  saveRDS(list(rows = ora, status = st), chk)
  ora_status[[cid]] <- st
  nrows <- if (is.null(ora)) 0L else nrow(ora)
  rm(rows, ora); gc(verbose = FALSE)
  vmrss <- as.numeric(gsub("[^0-9]", "", system(paste0("grep VmRSS /proc/", Sys.getpid(), "/status"), intern = TRUE)[1])) / 1024
  logmsg("ORA done", cid, "up:", length(up_e), st["up"], "| down:", length(dn_e), st["down"],
         "| rows:", nrows, "| VmRSS(MB):", round(vmrss))
}
logmsg("ORA phase complete")

## ---------- 4. GSEA ----------
# powered bulk contrasts from runlog
rl <- read.csv("analysis/step2_bulk/runlog.csv")
pow_bulk <- paste0(rl$gse[rl$n_case >= 5 & rl$n_control >= 5], "_",
                   rl$contrast[rl$n_case >= 5 & rl$n_control >= 5])
# powered scRNA contrasts (per task brief, from step-3 report donor counts >=3/arm)
pow_scr <- c(paste0("GSE235711_", rep(c("B","CD4_T","CD8_T","NK","basal","ciliated","endothelial",
  "fibroblast","glandular","goblet","mast","monocyte","myeloid","neutrophil","pDC","plasma",
  "secretory","secretory_club"), each = 2), "_", rep(c("CRSsNP_vs_control","CRSwNP_vs_control"), 18)),
  "GSE324445_ILC2_NP_vs_PB", "GSE324445_CD4_T_NP_vs_PB",
  "GSE249159_CD4_T_Polyp_vs_PBMC",
  "GSE270953_B_NP_vs_tonsil", "GSE270953_plasma_NP_vs_tonsil")
powered <- unique(c(pow_bulk, pow_scr))
writeLines(powered, file.path(OUT, "powered_contrasts.txt"))

# gene sets
logmsg("fetching msigdbr collections")
msig <- list(
  Hallmark = msigdbr(species = "Homo sapiens", collection = "H"),
  GO_BP    = msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"),
  Reactome = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  KEGG     = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")) # msigdbr>=26 renamed CP:KEGG -> CP:KEGG_LEGACY
gs_info <- do.call(rbind, lapply(names(msig), function(nm) {
  d <- msig[[nm]]
  data.frame(collection = nm, n_sets = length(unique(d$gs_name)),
             version = paste0("msigdbr_", as.character(packageVersion("msigdbr"))),
             stringsAsFactors = FALSE)
}))
write.csv(gs_info, file.path(OUT, "msigdb_collections.csv"), row.names = FALSE)
# msigdbr 26.1.1 dropped the entrez_gene column; use gene_symbol and rank stats by symbol
pw <- lapply(msig, function(d) split(as.character(d$gene_symbol), d$gs_name))
saveRDS(pw, file.path(OUT, "chk_pathways.rds"))
logmsg("collections ready")

gsea_status <- list()
for (cid in names(allc)) {
  chk <- file.path(OUT, "chk_gsea", paste0(cid, ".rds"))
  if (file.exists(chk)) { gsea_status[[cid]] <- readRDS(chk)$status; next }
  a <- allc[[cid]]; m <- a$dat
  if (!(cid %in% powered)) { saveRDS(list(rows = NULL, status = "skipped_not_powered"), chk); gsea_status[[cid]] <- "skipped_not_powered"; next }
  p <- pmax(m$p, 1e-300)
  stat <- setNames(sign(m$logFC) * (-log10(p)), m$symbol)
  stat <- stat[!is.na(names(stat)) & names(stat) != ""]          # drop missing symbols
  o <- order(abs(stat), decreasing = TRUE)                       # dedup: keep max |stat| per symbol
  stat <- stat[o[!duplicated(names(stat)[o])]]
  stat <- sort(stat, decreasing = TRUE)
  rows <- list()
  st <- "done"
  for (nm in names(pw)) {
    rr <- tryCatch(
      fgsea(pathways = pw[[nm]], stats = stat, minSize = 15, maxSize = 500, nPermSimple = 10000),
      error = function(e) { logmsg("GSEA FAIL", cid, nm, conditionMessage(e)); NULL })
    if (is.null(rr)) { st <- paste0("failed_", nm); next }
    if (nrow(rr) > 0) {
      rr <- as.data.frame(rr)
      rows[[nm]] <- data.frame(
        contrast = cid, layer = a$layer, collection = nm,
        term_id = rr$pathway, term_name = rr$pathway,
        NES = rr$NES, pvalue = rr$pval, FDR = rr$padj, size = rr$size,
        leading_edge_n = vapply(rr$leadingEdge, length, integer(1)),
        stringsAsFactors = FALSE)
    }
  }
  gsea <- if (length(rows)) do.call(rbind, rows) else NULL
  saveRDS(list(rows = gsea, status = st), chk)
  gsea_status[[cid]] <- st
  nrows <- if (is.null(gsea)) 0L else nrow(gsea)
  rm(rows, gsea, stat); gc(verbose = FALSE)
  logmsg("GSEA done", cid, st, "| rows:", nrows)
}
logmsg("GSEA phase complete")

## ---------- 5. assemble outputs ----------
ora_all <- do.call(rbind, lapply(names(allc), function(cid) readRDS(file.path(OUT, "chk_ora", paste0(cid, ".rds")))$rows))
gsea_all <- do.call(rbind, lapply(names(allc), function(cid) readRDS(file.path(OUT, "chk_gsea", paste0(cid, ".rds")))$rows))
if (is.null(ora_all)) ora_all <- data.frame(contrast=character(),layer=character(),direction=character(),database=character(),term_id=character(),term_name=character(),set_size=integer(),count=integer(),pvalue=numeric(),FDR=numeric(),gene_ratio=character(),genes=character())
if (is.null(gsea_all)) gsea_all <- data.frame(contrast=character(),layer=character(),collection=character(),term_id=character(),term_name=character(),NES=numeric(),pvalue=numeric(),FDR=numeric(),size=integer(),leading_edge_n=integer())
ora_all <- ora_all[order(ora_all$contrast, ora_all$direction, ora_all$database, ora_all$FDR), ]
gsea_all <- gsea_all[order(gsea_all$contrast, gsea_all$collection, gsea_all$FDR), ]
write.csv(ora_all, file.path(OUT, "enrichment_ora.csv"), row.names = FALSE)
write.csv(gsea_all, file.path(OUT, "enrichment_gsea.csv"), row.names = FALSE)

rl_out <- do.call(rbind, lapply(names(allc), function(cid) {
  a <- allc[[cid]]
  data.frame(contrast = cid, layer = a$layer, n_genes = a$n_genes,
             n_up_FDR005 = a$n_up, n_down_FDR005 = a$n_down, n_unmapped = a$n_unmapped,
             ora_up_status = ora_status[[cid]]["up"], ora_down_status = ora_status[[cid]]["down"],
             gsea_status = gsea_status[[cid]],
             notes = if (!(cid %in% powered)) "GSEA skipped: low power (bulk n_case<5 or n_control<5; scRNA <3 donors/arm)" else "",
             stringsAsFactors = FALSE, row.names = NULL)
}))
write.csv(rl_out, file.path(OUT, "enrichment_runlog.csv"), row.names = FALSE)
logmsg("WROTE: ora rows =", nrow(ora_all), "| gsea rows =", nrow(gsea_all), "| runlog rows =", nrow(rl_out))
