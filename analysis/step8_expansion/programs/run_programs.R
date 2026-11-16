#!/usr/bin/env Rscript
# Step 8 expansion ITEM 2: ssGSEA program scores per cohort (PolypAtlas)
# All numbers from produced files. ssGSEA is rank-based, so log2 vs count
# scale matrices are both acceptable (documented).
suppressPackageStartupMessages({ library(GSVA); library(msigdbr) })

BASE <- "analysis/step8_expansion/programs"
PREP <- "analysis/step2_bulk/prep"

## ---- program definitions (a priori panel) ----
h <- msigdbr(species = "human", collection = "H")
hall <- function(nm) unique(h$gene_symbol[h$gs_name == nm])
programs <- list(
  type2_custom      = c("CLC","CCL26","CCL13","IL4","IL5","IL13","IL1RL1","IL33","TSLP","POSTN","SERPINB2","ALOX15","PTGDR2","HPGDS","ALOX5","CYSLTR1","IL13RA2","EPX"),
  il13_response     = c("POSTN","SERPINB2","CLCA1","MUC5AC","CCL26","IL13RA2","ALOX15","DUOX1","ITLN1","TFF2"),
  interferon        = hall("HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  hypoxia           = hall("HALLMARK_HYPOXIA"),
  proliferation     = hall("HALLMARK_E2F_TARGETS"),
  barrier           = c("CDH1","OCLN","CLDN1","CLDN4","TJP1","DSC2","DSG3","CTNNA1","CTNNB1"),
  ciliated          = c("FOXJ1","DNAH5","DNAH11","CCDC39","RSPH1","HYDIN","TCTN1","CFAP221","SPAG6"),
  goblet_secretory  = c("MUC5AC","MUC5B","SPDEF","FOXA3","AGR2","TFF3","SCGB1A1","CLCA1"),
  neutrophil        = c("CSF3R","S100A8","S100A9","ELANE","MPO","CXCR2","FCGR3B"),
  mast              = c("TPSAB1","TPSB2","CPA3","KIT","HDC","MS4A2"),
  ilc2              = c("GATA3","IL1RL1","IL13","KLRG1","PTGDR2","IL17RB","ICOS"),
  eosinophil        = c("CLC","EPX","IL5RA","PRG2","PRG3","SIGLEC8","CCR3"),
  fibroblast_ecm    = c("COL1A1","COL3A1","DCN","LUM","FBLN1","POSTN","VCAN","MMP2")
)
sources <- c(
  type2_custom     = "literature-defined CRSwNP type-2 signature (canonical type-2 mediators + project marker panel)",
  il13_response    = "IL-13-induced airway epithelial response genes (literature-curated)",
  interferon       = "MSigDB Hallmark INTERFERON_GAMMA_RESPONSE (msigdbr)",
  hypoxia          = "MSigDB Hallmark HYPOXIA (msigdbr)",
  proliferation    = "MSigDB Hallmark E2F_TARGETS (msigdbr)",
  barrier          = "epithelial junction/barrier (canonical)",
  ciliated         = "ciliated cell canonical markers",
  goblet_secretory = "goblet/secretory canonical markers",
  neutrophil       = "neutrophil canonical markers",
  mast             = "mast cell canonical markers",
  ilc2             = "ILC2 canonical markers",
  eosinophil       = "eosinophil canonical markers",
  fibroblast_ecm   = "fibroblast/ECM canonical markers"
)
prog_def <- data.frame(program = names(programs), n_genes_total = sapply(programs, length),
                       source = sources[names(programs)], row.names = NULL)
write.csv(prog_def, file.path(BASE, "programs_def.csv"), row.names = FALSE)

## ---- cohort group rules (primary contrast per step2_bulk/runlog.csv) ----
## drop = groups excluded from the primary contrast (kept in scores file but unused in summary)
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
rule_text <- c(
  GSE136825="polyp (n=75) vs control inferior turbinate (n=28)",
  GSE150846="sorted NP-ILC2 (n=3) vs PB-ILC2 (n=3); stim arm excluded",
  GSE174681="sorted NP antibody-secreting B (n=3) vs peripheral naive B (n=4); NP naive B excluded",
  GSE179265="CRSwNP (n=17) vs healthy control (n=7)",
  GSE189690="N-ERD (n=13) vs non-N-ERD CRSwNP (n=13); no healthy control",
  GSE194282="NP (n=7) vs control uncinate tissue (n=7)",
  GSE198950="CRSwNP ethmoid (n=5) vs non-CRS inferior turbinate (n=4); CRSsNP excluded; TISSUE-CONFOUNDED",
  GSE213809="stimulated NP-ILC2 (n=14) vs stimulated PB-ILC2 (n=14); unstimulated excluded",
  GSE267331="dupilumab post vs pre, PAIRED within-subject (47 pairs, SINUS-52); placebo arm excluded",
  GSE313186="deployment-exposed CRS (n=10) vs unexposed control (n=10)",
  GSE330761="CRSwNP MCTC (n=3) vs CRSsNP MCTC (n=3); sorted mast cells, no healthy control",
  GSE338426="CRSwNP (n=14) vs control (n=8)",
  GSE36830="CRSwNP (n=12) vs healthy UT control (n=6); CRSsNP excluded",
  GSE69093="case (n=25) vs control (n=23) nasal epi+fibroblast primary cells; weak signal",
  GSE70898_70900="polyp T cells (n=12) vs normal nasal mucosa T cells (n=6); blood T ref excluded",
  GSE72713="CRSwNP (n=6, ECR+nonECR) vs control (n=3)"
)

scores_all <- list(); summ_all <- list(); runlog <- list()

for (co in names(rules)) {
  t0 <- Sys.time()
  r <- rules[[co]]
  d <- readRDS(file.path(PREP, paste0(co, ".rds")))
  e <- as.matrix(d$expr)
  m <- d$meta
  note0 <- paste0("dup_rows=", sum(duplicated(rownames(e))))

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
  note0 <- paste0(note0, "; symbol_mapped=", mapped)

  ## collapse duplicated symbols (max-mean), drop all-zero rows
  if (any(duplicated(rownames(e)))) {
    mm <- rowMeans(e, na.rm = TRUE)
    keep <- tapply(seq_len(nrow(e)), rownames(e), function(ix) ix[which.max(mm[ix])])
    e <- e[unlist(keep), , drop = FALSE]
  }
  e <- e[rowSums(e, na.rm = TRUE) > 0, , drop = FALSE]
  note0 <- paste0(note0, "; zero_rows_dropped; genes_after_qc=", nrow(e))

  ## program gene coverage
  cov <- sapply(programs, function(g) sum(g %in% rownames(e)))
  keep_p <- names(cov)[cov >= 3]
  skip_p <- names(cov)[cov < 3]

  sc <- NULL
  if (length(keep_p) > 0) {
    sets <- lapply(programs[keep_p], function(g) intersect(g, rownames(e)))
    esc <- gsva(ssgseaParam(e, sets), verbose = FALSE)   # rank-based; count or log2 scale both OK
    smp <- colnames(e)
    sc <- do.call(rbind, lapply(keep_p, function(p) data.frame(
      cohort = co, sample = smp,
      group = as.character(m$group[match(smp, rownames(m))]),
      program = p, score = as.numeric(esc[p, smp]), n_genes_used = cov[p],
      stringsAsFactors = FALSE)))
    write.csv(sc, file.path(BASE, "checkpoints", paste0(co, "_scores.csv")), row.names = FALSE)
    scores_all[[co]] <- sc
  }

  ## summary on primary-contrast samples
  n_primary <- NA_integer_
  if (!is.null(sc)) {
    mm <- m
    if (!is.null(r$arm)) mm <- mm[mm$treatment == r$arm & mm$time %in% c("Pre-Treatment","Post-Treatment"), ]
    if (!is.null(r$drop)) mm <- mm[!mm$group %in% r$drop, ]
    n_primary <- nrow(mm)
    sco <- sc[sc$sample %in% rownames(mm), ]
    rows <- list()
    for (p in keep_p) {
      v <- sco[sco$program == p, ]
      xc <- v$score[v$group == r$case]; xn <- v$score[v$group == r$ctrl]
      if (r$paired) {
        key <- mm$subject[match(v$sample, rownames(mm))]
        dfr <- data.frame(sample = v$sample, group = v$group, score = v$score, subject = key)
        pw <- dfr[dfr$group == r$case, ]; pv <- dfr[dfr$group == r$ctrl, ]
        pw <- pw[order(pw$subject), ]; pv <- pv[order(pv$subject), ]
        stopifnot(all(pw$subject == pv$subject))
        diffs <- pw$score - pv$score
        pv_p <- wilcox.test(diffs)$p.value
        rows[[p]] <- data.frame(cohort = co, program = p, n_case = length(xc), n_ctrl = length(xn),
          mean_case = mean(xc), mean_ctrl = mean(xn), diff = mean(xc) - mean(xn),
          wilcox_p = pv_p, test = "paired signed-rank (post-minus-pre per subject)")
      } else {
        pv_p <- wilcox.test(xc, xn)$p.value
        rows[[p]] <- data.frame(cohort = co, program = p, n_case = length(xc), n_ctrl = length(xn),
          mean_case = mean(xc), mean_ctrl = mean(xn), diff = mean(xc) - mean(xn),
          wilcox_p = pv_p, test = "Wilcoxon rank-sum")
      }
    }
    summ <- do.call(rbind, rows)
    summ$FDR <- p.adjust(summ$wilcox_p, method = "BH")
    summ$note <- rule_text[co]
    summ$note[co == "GSE198950"] <- paste0(rule_text["GSE198950"], " -> do NOT interpret as clean disease effect")
    summ$note[co == "GSE69093"]  <- paste0(rule_text["GSE69093"], " -> primary contrast was celltype-adjusted; program diff is unadjusted")
    summ_all[[co]] <- summ
    write.csv(summ, file.path(BASE, "checkpoints", paste0(co, "_summary.csv")), row.names = FALSE)
  }

  runlog[[co]] <- data.frame(
    cohort = co, group_rule = rule_text[co], n_samples_total = ncol(e),
    n_samples_primary = n_primary,
    programs_scored = paste(keep_p, collapse = ";"),
    programs_skipped = if (length(skip_p)) paste(skip_p, collapse = ";") else "",
    qc_note = note0,
    elapsed_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    stringsAsFactors = FALSE)
  cat(sprintf("[%s] scored %d programs, skipped %d | %s\n",
      co, length(keep_p), length(skip_p), runlog[[co]]$qc_note), flush = TRUE)
}

program_scores  <- do.call(rbind, scores_all)
program_summary <- do.call(rbind, summ_all)
runlog_df       <- do.call(rbind, runlog)
write.csv(program_scores,  file.path(BASE, "program_scores.csv"),  row.names = FALSE)
write.csv(program_summary, file.path(BASE, "program_summary.csv"), row.names = FALSE)
write.csv(runlog_df,       file.path(BASE, "programs_runlog.csv"), row.names = FALSE)
cat("DONE. scores:", nrow(program_scores), "summary rows:", nrow(program_summary), "\n")
