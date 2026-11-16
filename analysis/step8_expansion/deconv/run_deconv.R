#!/usr/bin/env Rscript
# Step 8 expansion ITEM 3: bulk immune deconvolution, MCP-counter-STYLE marker scoring (PolypAtlas)
# CRITICAL LABELING: every output/legend says "MCP-counter-style marker scoring";
# values are "abundance estimates" -- NEVER presented as official MCPcounter package output.
# Formula (per task): score_pop(s) = mean_over_markers( log2( x_lin + 1 ) )
#   = log2 of the geometric mean of (x_lin + 1), where x_lin is linear-scale
#   expression: raw counts used as-is; log2-scale matrices back-transformed via 2^expr.
suppressPackageStartupMessages({ library(utils) })

BASE  <- "analysis/step8_expansion/deconv"
PREP  <- "analysis/step2_bulk/prep"
MCP_URL <- "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
MCP_FILE <- "/tmp/mcp_genes.txt"

## ---- load marker set ----
if (!file.exists(MCP_FILE)) stop("marker file missing; download failed")
mcp <- read.table(MCP_FILE, sep = "\t", header = TRUE, quote = "\"",
                  stringsAsFactors = FALSE, check.names = FALSE)
colnames(mcp) <- c("symbol", "population", "entrez", "ensembl")
mcp$symbol <- gsub('"', "", mcp$symbol); mcp$population <- gsub('"', "", mcp$population)
marker_src <- paste0("MCP-counter published signature file (", MCP_URL, "); ",
                     nrow(mcp), " marker genes / ", length(unique(mcp$population)), " populations; downloaded 2026-10-18")
pops <- sort(unique(mcp$population))
cat("Populations:", paste(pops, collapse = " | "), "\n")
cat(marker_src, "\n")

pop_markers <- split(mcp$symbol, mcp$population)
pop_info <- data.frame(population = pops,
                       n_markers_total = sapply(pop_markers[pops], length),
                       method = "MCP-counter-style marker scoring (abundance estimates; NOT official MCPcounter output)",
                       stringsAsFactors = FALSE)
write.csv(pop_info, file.path(BASE, "marker_set_info.csv"), row.names = FALSE)

## ---- cohort group rules (same primary contrasts as programs) ----
rules <- list(
  GSE136825     = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE150846     = list(drop = "stim",      case = "case", ctrl = "control", paired = FALSE),
  GSE174681     = list(drop = "other",     case = "case", ctrl = "control", paired = FALSE),
  GSE179265     = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE189690     = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE194282     = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE198950     = list(drop = "other",     case = "case", ctrl = "control", paired = FALSE),
  GSE213809     = list(drop = "other",     case = "case", ctrl = "control", paired = FALSE),
  GSE267331     = list(drop = NULL,        case = "case", ctrl = "control", paired = TRUE,  arm = "Dupilumab"),
  GSE313186     = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE330761     = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE338426     = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE36830      = list(drop = "other",     case = "case", ctrl = "control", paired = FALSE),
  GSE69093      = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE),
  GSE70898_70900= list(drop = "blood_ref", case = "case", ctrl = "control", paired = FALSE),
  GSE72713      = list(drop = NULL,        case = "case", ctrl = "control", paired = FALSE)
)

scores_all <- list(); summ_all <- list(); runlog <- list()

for (co in names(rules)) {
  t0 <- Sys.time()
  r <- rules[[co]]
  d <- readRDS(file.path(PREP, paste0(co, ".rds")))
  e <- as.matrix(d$expr)
  m <- d$meta
  scalev <- as.character(d$scale_verdict)

  ## map rownames to gene symbols via d$symbols (some cohorts store Ensembl/probe IDs in rownames)
  sym <- d$symbols
  mapped <- FALSE
  if (!is.null(sym) && length(sym) == nrow(e)) {
    ns <- names(sym)
    if (!is.null(ns) && all(ns %in% rownames(e))) sym <- sym[rownames(e)]
    good <- !is.na(sym) & nzchar(as.character(sym))
    e <- e[good, , drop = FALSE]
    rownames(e) <- as.character(sym[good])
    mapped <- TRUE
  }
  map_note <- paste0("symbol_mapped=", mapped)

  ## collapse duplicated symbols (max-mean), drop all-zero rows
  if (any(duplicated(rownames(e)))) {
    mmx <- rowMeans(e, na.rm = TRUE)
    keep <- tapply(seq_len(nrow(e)), rownames(e), function(ix) ix[which.max(mmx[ix])])
    e <- e[unlist(keep), , drop = FALSE]
  }
  e <- e[rowSums(e, na.rm = TRUE) > 0, , drop = FALSE]

  ## linear-scale expression for the log2(x+1) formula
  e_lin <- if (scalev %in% c("log2_normalized", "linear_normalized_log2t")) 2^e else e

  pop_used <- list(); S <- matrix(NA_real_, nrow = ncol(e), ncol = length(pops))
  colnames(S) <- pops; rownames(S) <- colnames(e)
  for (p in pops) {
    g <- intersect(pop_markers[[p]], rownames(e_lin))
    if (length(g) >= 1) {
      S[, p] <- colMeans(log2(e_lin[g, , drop = FALSE] + 1))
      pop_used[[p]] <- length(g)
    }
  }
  keep_pops <- names(pop_used)
  skip_pops <- setdiff(pops, keep_pops)

  sc <- NULL
  if (length(keep_pops) > 0) {
    smp <- rownames(S)
    sc <- do.call(rbind, lapply(keep_pops, function(p) data.frame(
      cohort = co, sample = smp,
      group = as.character(m$group[match(smp, rownames(m))]),
      method = "MCP-counter-style marker scoring",
      population = p, score = S[, p], n_markers_used = pop_used[[p]],
      stringsAsFactors = FALSE)))
    colnames(sc) <- c("cohort","sample","group","method","population",
                      "abundance_estimate","n_markers_used")
    write.csv(sc, file.path(BASE, "checkpoints", paste0(co, "_deconv_scores.csv")), row.names = FALSE)
    scores_all[[co]] <- sc
  }

  ## summary on primary-contrast samples
  n_primary <- NA_integer_
  if (!is.null(sc)) {
    mm2 <- m
    if (!is.null(r$arm)) mm2 <- mm2[mm2$treatment == r$arm & mm2$time %in% c("Pre-Treatment","Post-Treatment"), ]
    if (!is.null(r$drop)) mm2 <- mm2[!mm2$group %in% r$drop, ]
    n_primary <- nrow(mm2)
    rows <- list()
    for (p in keep_pops) {
      v <- sc[sc$sample %in% rownames(mm2), ]
      xc <- v$abundance_estimate[v$population == p & v$group == r$case]
      xn <- v$abundance_estimate[v$population == p & v$group == r$ctrl]
      if (r$paired) {
        vv <- v[v$population == p, ]
        key <- mm2$subject[match(vv$sample, rownames(mm2))]
        kc <- key[vv$group == r$case]; kn <- key[vv$group == r$ctrl]
        zc <- vv$abundance_estimate[vv$group == r$case]; zn <- vv$abundance_estimate[vv$group == r$ctrl]
        oc <- order(kc); on <- order(kn)
        stopifnot(all(kc[oc] == kn[on]))
        diffs <- zc[oc] - zn[on]
        pv_p <- wilcox.test(diffs)$p.value; tst <- "paired signed-rank (post-minus-pre per subject)"
      } else {
        pv_p <- wilcox.test(xc, xn)$p.value; tst <- "Wilcoxon rank-sum"
      }
      rows[[p]] <- data.frame(cohort = co, population = p, n_case = length(xc), n_ctrl = length(xn),
        mean_case = mean(xc), mean_ctrl = mean(xn), diff = mean(xc) - mean(xn),
        wilcox_p = pv_p, test = tst, n_markers_used = pop_used[[p]],
        method = "MCP-counter-style marker scoring (abundance estimates; NOT official MCPcounter output)")
    }
    summ <- do.call(rbind, rows)
    summ$FDR <- p.adjust(summ$wilcox_p, method = "BH")
    summ_all[[co]] <- summ
    write.csv(summ, file.path(BASE, "checkpoints", paste0(co, "_deconv_summary.csv")), row.names = FALSE)
  }

  runlog[[co]] <- data.frame(
    cohort = co, marker_set_source = marker_src, scale_verdict = scalev,
    scale_handling = if (scalev %in% c("log2_normalized","linear_normalized_log2t"))
      "matrix already log2; back-transformed 2^expr before log2(x+1)" else "counts used as-is in log2(x+1)",
    n_samples_total = ncol(e), n_samples_primary = n_primary,
    populations_scored = paste(keep_pops, collapse = ";"),
    populations_skipped = if (length(skip_pops)) paste(skip_pops, collapse = ";") else "",
    qc_note = map_note,
    elapsed_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    stringsAsFactors = FALSE)
  cat(sprintf("[%s] scored %d populations | elapsed %.1f min\n",
      co, length(keep_pops), runlog[[co]]$elapsed_min), flush = TRUE)
}

deconv_scores  <- do.call(rbind, scores_all)
deconv_summary <- do.call(rbind, summ_all)
runlog_df      <- do.call(rbind, runlog)
write.csv(deconv_scores,  file.path(BASE, "deconv_scores.csv"),  row.names = FALSE)
write.csv(deconv_summary, file.path(BASE, "deconv_summary.csv"), row.names = FALSE)
write.csv(runlog_df,      file.path(BASE, "deconv_runlog.csv"),  row.names = FALSE)
cat("DONE. scores:", nrow(deconv_scores), "summary rows:", nrow(deconv_summary), "\n")
