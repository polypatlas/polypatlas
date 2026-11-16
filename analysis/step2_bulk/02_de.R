#!/usr/bin/env Rscript
# 02_de.R — parameterized DE driver for step2 bulk cohorts.
# Usage: Rscript 02_de.R <gse>
# Reads spec from cohort_spec.csv; writes de/<GSE>_DE.csv (+ extras), appends runlog.csv.
suppressPackageStartupMessages({
  library(limma); library(DESeq2); library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
GSE <- args[1]
OUT <- "analysis/step2_bulk"
SPEC <- read.csv(file.path(OUT, "cohort_spec.csv"), stringsAsFactors = FALSE, colClasses = "character")
sp <- SPEC[SPEC$gse == GSE, ]
stopifnot(nrow(sp) == 1)
sp <- as.list(sp)

logf <- file.path(OUT, "log", paste0("02_de_", GSE, ".log"))
)
sink(logf, split = TRUE)
on.exit(sink())
cat("===== DE driver:", GSE, "=====\n")

# ---------- load expression ----------
es <- readRDS(file.path(OUT, "meta", "geo_dl", paste0(GSE, ".rds")))
plats <- sapply(es, annotation)
pi <- as.integer(sp$platform_idx)
cat("platforms:", paste(plats, collapse = ", "), "-> using idx", pi, plats[pi], "\n")
e <- es[[pi]]
mat <- exprs(e)
pd <- pData(e)
fd <- fData(e)
cat(sprintf("matrix: %d features x %d samples\n", nrow(mat), ncol(mat)))

# ---------- scale check (three-way) ----------
q <- quantile(mat, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
negfrac <- mean(mat < 0, na.rm = TRUE)
is_int <- all(mat == round(mat), na.rm = TRUE)
cat(sprintf("scale: min=%.3f p25=%.3f med=%.3f p75=%.3f max=%.3f negfrac=%.4f integer=%s\n",
            q[1], q[2], q[3], q[4], q[5], negfrac, is_int))
pub_note <- sp$pub_note
scale_verdict <- if (is_int && q[1] >= 0) "raw_counts"
                 else if (negfrac > 0.01) "log2_normalized"
                 else if (q[5] > 50) "linear_normalized"
                 else "log2_normalized?"
cat("publication note:", pub_note, "\n")
cat("SCALE VERDICT:", scale_verdict, "\n")

# ---------- groups ----------
gcol <- sp$group_col
case_v <- trimws(strsplit(sp$case_values, ";")[[1]])
ctrl_v <- trimws(strsplit(sp$control_values, ";")[[1]])
stopifnot(gcol %in% names(pd))
grp <- as.character(pd[[gcol]])
cat("group column:", gcol, "| unique:", paste(unique(grp), collapse = " / "), "\n")
keep <- grp %in% c(case_v, ctrl_v)
grp <- factor(ifelse(grp[keep] %in% case_v, "case", "control"), levels = c("control", "case"))
mat <- mat[, keep, drop = FALSE]; pd <- pd[keep, , drop = FALSE]
n_case <- sum(grp == "case"); n_ctrl <- sum(grp == "control")
cat(sprintf("primary contrast: case(n=%d) vs control(n=%d)\n", n_case, n_ctrl))
status <- "OK"
if (n_case < 3 || n_ctrl < 3) { status <- "LOW_POWER"; cat("!! LOW_POWER: n<3 on a side\n") }

method <- if (scale_verdict == "raw_counts") "DESeq2" else "limma"
cat("method:", method, "\n")

# ---------- gene symbols ----------
symcol <- sp$symbol_col
syms <- if (symcol %in% names(fd)) as.character(fd[[symcol]]) else rep(NA, nrow(mat))
if (all(is.na(syms)) && "ENTREZID" %in% names(fd)) {
  ez <- suppressWarnings(as.character(mapIds(org.Hs.eg.db, keys = as.character(fd$ENTREZID),
                                             column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")))
  syms <- ifelse(is.na(syms), ez, syms)
}
# split multi-symbol probes, take first
syms <- trimws(sapply(strsplit(syms, " /// | // | \\|"), function(x) x[1]))
syms[syms %in% c("", "NA", "---")] <- NA
cat(sprintf("symbol coverage: %d/%d features\n", sum(!is.na(syms)), length(syms)))

# ---------- DE ----------
paired <- !is.na(sp$paired_subject) && sp$paired_subject != "" && sp$paired_subject %in% names(pd)
if (method == "DESeq2") {
  dds <- DESeqDataSetFromMatrix(countData = round(mat), colData = data.frame(grp = grp),
                                design = ~grp)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("grp", "case", "control")))
  out <- data.frame(gene = syms, logFC = res$log2FoldChange, p = res$pvalue, FDR = res$padj)
} else {
  dmat <- if (paired) {
    subj <- factor(pd[[sp$paired_subject]])
    cat("paired design, subject col:", sp$paired_subject, "\n")
    model.matrix(~0 + grp + subj)
  } else model.matrix(~0 + grp)
  fit <- lmFit(mat, dmat)
  cm <- if (paired) makeContrasts("grpcase - grpcontrol", levels = dmat)
         else makeContrasts("grpcase - grpcontrol", levels = dmat)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  out <- data.frame(gene = syms, logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val)
}
out <- out[order(out$FDR), ]
dir.create(file.path(OUT, "de"), showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(OUT, "de", paste0(GSE, "_DE.csv")), row.names = FALSE)
cat(sprintf("wrote primary DE table: %d genes\n", nrow(out)))

# ---------- extra contrasts ----------
xc <- sp$extra_contrasts
if (!is.na(xc) && xc != "") {
  # format: "name|col|caseA;caseB|ctrlA;ctrlB ## name2|..." — simplified: same col, new groupings
  for (blk in strsplit(xc, " ## ")[[1]]) {
    parts <- strsplit(blk, "\\|")[[1]]
    xn <- parts[1]; xcase <- trimws(strsplit(parts[2], ";")[[1]]); xctrl <- trimws(strsplit(parts[3], ";")[[1]])
    allg <- as.character(pd[[gcol]])
    kk <- allg %in% c(xcase, xctrl)
    if (sum(allg[kk] %in% xcase) < 2 || sum(allg[kk] %in% xctrl) < 2) {
      cat(sprintf("extra contrast %s skipped: n too small\n", xn)); next
    }
    xm <- mat[, kk, drop = FALSE]
    xg <- factor(ifelse(allg[kk] %in% xcase, "case", "control"), levels = c("control", "case"))
    xd <- model.matrix(~0 + xg)
    xf <- eBayes(contrasts.fit(lmFit(xm, xd), makeContrasts("xgcase - xgcontrol", levels = xd)))
    xt <- topTable(xf, number = Inf, sort.by = "none")
    xo <- data.frame(gene = syms, logFC = xt$logFC, p = xt$P.Value, FDR = xt$adj.P.Val)
    xo <- xo[order(xo$FDR), ]
    write.csv(xo, file.path(OUT, "de", sprintf("%s_DE_%s.csv", GSE, xn)), row.names = FALSE)
    cat(sprintf("wrote extra contrast %s: n=%d/%d\n", xn, sum(xg == "case"), sum(xg == "control")))
  }
}

# ---------- runlog ----------
rl <- data.frame(gse = GSE, n_case = n_case, n_control = n_ctrl,
                 groups_used = paste0(sp$case_values, " VS ", sp$control_values),
                 scale_verdict = scale_verdict, method = method, status = status,
                 note = sp$note, stringsAsFactors = FALSE)
rlf <- file.path(OUT, "runlog.csv")
write.table(rl, rlf, sep = ",", row.names = FALSE, col.names = !file.exists(rlf), append = file.exists(rlf))
cat("runlog appended, status:", status, "\nDONE\n")
