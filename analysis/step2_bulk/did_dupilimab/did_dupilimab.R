#!/usr/bin/env Rscript
# did_dupilimab.R — treatment x time (difference-in-differences) analysis of
# GSE267331 (SINUS-52 nasal brushings, GPL570), the analysis the main pipeline
# did NOT run (main = within-arm paired limma per arm).
# Design: complete-pair subjects only (Dupilumab 47, Placebo 19).
#   1) probes -> genes by the pre-specified MAX-VARIANCE rule (variance across
#      the 132 paired samples; same rule as the probe-collapse sensitivity
#      analysis; no test statistic involved in probe choice)
#   2) per-subject paired difference  d = post - pre  (log2 scale)
#   3) limma on d with design ~ arm; coefficient = mean(d|D) - mean(d|P)
#      = DiD estimate; BH FDR across genes.
suppressPackageStartupMessages({ library(limma) })
BASE <- "analysis/step2_bulk"
OUT  <- file.path(BASE, "did_dupilimab")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

sp <- readRDS(file.path(BASE, "prep", "GSE267331.rds"))
meta <- sp$meta
meta <- meta[meta$time %in% c("Pre-Treatment", "Post-Treatment"), ]
expr <- sp$expr[, rownames(meta)]
stopifnot(nrow(meta) == 132)
tab <- table(meta$subject, meta$time)
oks <- rownames(tab)[tab[, "Pre-Treatment"] > 0 & tab[, "Post-Treatment"] > 0]
meta <- meta[meta$subject %in% oks, ]
expr <- expr[, rownames(meta)]
arm_of <- tapply(as.character(meta$treatment), meta$subject, function(x) x[1])
cat("subjects:", length(oks), " Dupilumab:", sum(arm_of == "Dupilumab"),
    " Placebo:", sum(arm_of == "Placebo"), "\n")

## probes -> genes, max-variance rule
sym <- sp$symbols[rownames(expr)]
keep <- !is.na(sym) & sym != ""
expr_k <- expr[keep, ]; sym_k <- sym[keep]
v <- if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowVars(expr_k) else apply(expr_k, 1, var)
o <- order(-v, rownames(expr_k))
expr_k <- expr_k[o, ]; sym_k <- sym_k[o]
first <- !duplicated(sym_k)
gexpr <- expr_k[first, ]; rownames(gexpr) <- sym_k[first]
cat("gene-level matrix:", nrow(gexpr), "genes x", ncol(gexpr), "samples\n")

## per-subject paired differences
pre_id  <- rownames(meta)[meta$time == "Pre-Treatment"]
post_id <- rownames(meta)[meta$time == "Post-Treatment"]
subj_pre  <- meta[pre_id, "subject"];  subj_post <- meta[post_id, "subject"]
stopifnot(all(subj_pre == subj_post))            # prep ordering check
D <- gexpr[, post_id] - gexpr[, pre_id]
colnames(D) <- subj_post
arm <- factor(arm_of[colnames(D)], levels = c("Placebo", "Dupilumab"))

design <- model.matrix(~ arm)
fit <- eBayes(lmFit(D, design))
tt <- topTable(fit, coef = "armDupilumab", number = Inf, sort.by = "none")
out <- data.frame(symbol = rownames(tt), DiD_logFC = tt$logFC,
                  AveExpr = rowMeans(gexpr)[rownames(tt)],
                  t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val)
out <- out[order(out$FDR, -abs(out$DiD_logFC)), ]
write.csv(out, file.path(OUT, "did_de.csv"), row.names = FALSE)
write.csv(head(out, 30), file.path(OUT, "did_top30.csv"), row.names = FALSE)
n_sig <- sum(out$FDR < 0.05, na.rm = TRUE)
cat("DiD FDR<0.05 genes:", n_sig, "of", nrow(out), "\n")

## overlap with within-arm dupilumab DEGs (frozen main result)
main_d <- read.csv(file.path(BASE, "de", "GSE267331_dupilumab_post_vs_pre.csv"))
main_deg <- main_d$symbol[!is.na(main_d$FDR) & main_d$FDR < 0.05]
did_deg <- out$symbol[out$FDR < 0.05]
cat("within-D-arm DEGs:", length(main_deg), "| DiD DEGs among them:",
    length(intersect(main_deg, did_deg)), "\n")

## 8 pre-defined type-2 markers (step2 panel; step6 validation set)
panel <- c("CLC", "CCL26", "CCL13", "IL33", "TSLP", "POSTN", "SERPINB2", "ALOX15")
pl <- read.csv(file.path(BASE, "de", "GSE267331_placebo_post_vs_pre.csv"))
mk <- data.frame(gene = panel,
                 D_within_logFC = main_d$logFC[match(panel, main_d$symbol)],
                 D_within_FDR   = main_d$FDR[match(panel, main_d$symbol)],
                 P_within_logFC = pl$logFC[match(panel, pl$symbol)],
                 P_within_FDR   = pl$FDR[match(panel, pl$symbol)],
                 DiD_logFC      = out$DiD_logFC[match(panel, out$symbol)],
                 DiD_FDR        = out$FDR[match(panel, out$symbol)])
write.csv(mk, file.path(OUT, "markers_did.csv"), row.names = FALSE)
print(mk)

## ---- type-2 program score DiD (scores = frozen step8 ssGSEA output) -------
sc <- read.csv(file.path(BASE, "..", "step8_expansion", "programs", "program_scores.csv"),
               stringsAsFactors = FALSE)
sc <- sc[sc$cohort == "GSE267331" & sc$program == "type2_custom", ]
m2 <- meta; m2$sample <- rownames(m2)
sc <- merge(sc, m2[, c("sample", "subject", "time", "treatment")], by = "sample")
wide <- reshape(sc[, c("subject", "treatment", "time", "score")],
                idvar = c("subject", "treatment"), timevar = "time", direction = "wide")
wide$diff <- wide[["score.Post-Treatment"]] - wide[["score.Pre-Treatment"]]
wide <- wide[!is.na(wide$diff), ]
write.csv(wide, file.path(OUT, "program_score_subject_diffs.csv"), row.names = FALSE)
dD <- wide$diff[wide$treatment == "Dupilumab"]; dP <- wide$diff[wide$treatment == "Placebo"]
wt <- wilcox.test(dD, dP, conf.int = TRUE, exact = FALSE)
r_rb <- 2 * unname(wt$statistic) / (length(dD) * length(dP)) - 1
prog <- data.frame(program = "type2_custom", n_D = length(dD), n_P = length(dP),
                   median_diff_D = median(dD), median_diff_P = median(dP),
                   HL_estimate_DminusP = unname(wt$estimate),
                   HL_lo = wt$conf.int[1], HL_hi = wt$conf.int[2],
                   wilcox_p = wt$p.value, rank_biserial_r = r_rb)
write.csv(prog, file.path(OUT, "program_score_did.csv"), row.names = FALSE)
print(prog)
cat("DID DONE\n")
