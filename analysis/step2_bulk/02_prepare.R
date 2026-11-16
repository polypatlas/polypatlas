#!/usr/bin/env Rscript
# 02_prepare.R — per-cohort data prep: clean expression + group metadata + symbols.
# Writes prep/<GSE>.rds = list(expr, meta, symbols, scale_verdict, pub_note, method,
#                              primary, extras, note)
# meta: data.frame(sample, group[case/control], cov1..) aligned to expr columns.
suppressPackageStartupMessages({ library(org.Hs.eg.db) })
OUT <- "analysis/step2_bulk"
SUP <- file.path(OUT, "meta", "suppl")
dir.create(file.path(OUT, "prep"), showWarnings = FALSE, recursive = TRUE)

save_prep <- function(gse, expr, meta, symbols, scale_verdict, pub_note, method,
                      primary, extras = list(), note = "") {
  stopifnot(ncol(expr) == nrow(meta), length(symbols) == nrow(expr))
  names(symbols) <- rownames(expr)
  meta$sample <- colnames(expr)
  p <- list(expr = expr, meta = meta, symbols = symbols, scale_verdict = scale_verdict,
            pub_note = pub_note, method = method, primary = primary, extras = extras, note = note)
  saveRDS(p, file.path(OUT, "prep", paste0(gse, ".rds")))
  cat(sprintf("[prep] %s: %d genes x %d samples | scale=%s | method=%s | primary %s n=%d vs %d | %s\n",
              gse, nrow(expr), ncol(expr), scale_verdict, method,
              primary$name, sum(meta$group == "case"), sum(meta$group == "control"), note))
}

ens2sym <- function(ids) {
  ids <- gsub('^"|"$', "", trimws(as.character(ids)))
  ids <- sub("\\.[0-9]+$", "", ids)
  suppressWarnings(mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
                          keytype = "ENSEMBL", multiVals = "first"))
}
affy_sym <- function(probes, pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  suppressWarnings(mapIds(get(pkg), keys = as.character(probes), column = "SYMBOL",
                          keytype = "PROBEID", multiVals = "first"))
}
ilmn_sym <- function(probes) {
  suppressPackageStartupMessages(library(illuminaHumanv4.db))
  suppressWarnings(mapIds(illuminaHumanv4.db, keys = as.character(probes), column = "SYMBOL",
                          keytype = "PROBEID", multiVals = "first"))
}

## ---------------- 1. GSE338426 (counts csv) ----------------
{
  gse <- "GSE338426"
  d <- read.csv(gzfile(file.path(SUP, gse, "GSE338426_Expression.csv.gz")),
                check.names = FALSE)
  rc <- grep(":read count$", names(d), value = TRUE)
  samp <- sub(":read count$", "", rc)
  keepc <- samp != "" & samp != "read_count" & !is.na(samp)
  rc <- rc[keepc]; samp <- samp[keepc]
  mat <- as.matrix(d[, rc]); rownames(mat) <- d$Gene_ID; colnames(mat) <- samp
  storage.mode(mat) <- "integer"
  ph <- read.csv(file.path(OUT, "meta", "GSE338426_GPL24676_pheno.csv"), row.names = 1)
  grp <- ifelse(ph$title %in% paste0("EosNP_", 1:9) | ph$title %in% paste0("Non_EP_", 1:5), "case",
         ifelse(grepl("^Control", ph$title), "control", NA))
  meta <- data.frame(group = grp[match(samp, ph$title)], row.names = samp)
  stopifnot(!any(is.na(meta$group)))
  ex <- list(name = "e_vs_ctrl", case = "eCRSwNP", control = "Control",
             def = list(eCRSwNP = paste0("EosNP_", 1:9), Control = paste0("Control", 1:8)))
  nx <- list(name = "ne_vs_ctrl", case = "neCRSwNP", control = "Control",
             def = list(neCRSwNP = paste0("Non_EP_", 1:5), Control = paste0("Control", 1:8)))
  save_prep(gse, mat, meta, ens2sym(rownames(mat)), "raw_counts",
            "Deposited read-count columns (Expression.csv.gz); FPKM columns ignored",
            "deseq2", list(name = "CRSwNP_vs_ctrl", case = "case", control = "control"),
            extras = list(
              list(name = "eCRSwNP_vs_ctrl", case_vals = paste0("EosNP_", 1:9), ctrl_vals = paste0("Control", 1:8), title_map = ph$title),
              list(name = "neCRSwNP_vs_ctrl", case_vals = paste0("Non_EP_", 1:5), ctrl_vals = paste0("Control", 1:8), title_map = ph$title)),
            note = "endotype labels from deposited metadata; extras use title-based groups")
}

## ---------------- 2. GSE198950 (counts, _count cols) ----------------
{
  gse <- "GSE198950"
  con <- gzfile(file.path(SUP, gse, "GSE198950_CRS_GeneCount.txt.gz"), "rt")
  hdr <- strsplit(readLines(con, n = 1), "\t")[[1]]; close(con)
  cc <- grep("_count$", hdr, value = TRUE)
  sids <- sub("_count$", "", cc)
  d <- read.delim(gzfile(file.path(SUP, gse, "GSE198950_CRS_GeneCount.txt.gz")),
                  check.names = FALSE)
  mat <- as.matrix(d[, cc]); rownames(mat) <- d$GeneName; colnames(mat) <- sids
  storage.mode(mat) <- "integer"
  ph <- read.csv(file.path(OUT, "meta", "GSE198950_GPL20301_pheno.csv"), row.names = 1)
  # map sample id (description) -> GSM -> diagnosis
  desc2diag <- setNames(ph$characteristics_ch1.1, ph$description)
  diag <- desc2diag[sids]
  stopifnot(!any(is.na(diag)))
  meta <- data.frame(group = ifelse(grepl("CRSwNP", diag), "case",
                             ifelse(grepl("Non CRS", diag), "control", "other")),
                     diagnosis = diag, row.names = sids)
  eos <- setNames(ph$characteristics_ch1.2, ph$description)[sids]
  save_prep(gse, mat, meta, toupper(rownames(mat)), "raw_counts",
            "Deposited _count columns (gene-level); _rpkm columns ignored",
            "deseq2", list(name = "CRSwNP_vs_NonCRS", case = "case", control = "control"),
            extras = list(
              list(name = "CRSsNP_vs_NonCRS", case_vals = sids[grepl("CRSsNP", diag)],
                   ctrl_vals = sids[grepl("Non CRS", diag)]),
              list(name = "EOS_high_vs_low", case_vals = sids[grepl("Over 10", eos)],
                   ctrl_vals = sids[grepl("Under 10", eos)])),
            note = "tissue: 16 ethmoidal + 4 inferior turbinate (the 4 Non-CRS); eos split tissue EOS/HPF")
}

## ---------------- 3. GSE136825 (counts, RHP bam cols) ----------------
{
  gse <- "GSE136825"
  d <- read.delim(gzfile(file.path(SUP, gse, "GSE136825_genecounts_20190903.txt.gz")),
                  check.names = FALSE, quote = "")
  names(d) <- gsub('^"|"$', "", names(d))
  cc <- grep("\\.bam$", names(d), value = TRUE)
  sids <- sub("\\.bam$", "", cc)
  mat <- as.matrix(d[, cc]); rownames(mat) <- d$Geneid; colnames(mat) <- sids
  storage.mode(mat) <- "integer"
  ph <- read.csv(file.path(OUT, "meta", "GSE136825_GPL20301_pheno.csv"), row.names = 1)
  desc2tissue <- setNames(ph$characteristics_ch1, ph$description)
  tis <- desc2tissue[sids]
  stopifnot(!any(is.na(tis)))
  meta <- data.frame(group = ifelse(grepl("Nasal Polyp", tis), "case", "control"),
                     tissue = tis, row.names = sids)
  save_prep(gse, mat, meta, ens2sym(rownames(mat)), "raw_counts",
            "Deposited featureCounts gene counts (RHP*.bam columns)",
            "deseq2", list(name = "polyp_vs_ctrl", case = "case", control = "control"),
            note = "largest CRSwNP bulk cohort (75 polyp vs 28 control inferior turbinate)")
}

## ---------------- 4. GSE36830 (GPL570, linear -> log2) ----------------
{
  gse <- "GSE36830"
  es <- readRDS(file.path(OUT, "meta", "matrix", paste0(gse, ".rds")))[[1]]
  mat <- log2(exprs(es) + 1)
  ph <- read.csv(file.path(OUT, "meta", "GSE36830_GPL570_pheno.csv"), row.names = 1)
  grp <- ifelse(ph$characteristics_ch1.1 == "crs type: CRS with nasal polyps (CRSwNP)", "case",
         ifelse(ph$characteristics_ch1 == "tissue: uncinate tissue (UT)" & ph$characteristics_ch1.1 == "", "control", "other"))
  meta <- data.frame(group = grp, row.names = rownames(ph))
  save_prep(gse, mat, meta, affy_sym(rownames(mat), "hgu133plus2.db"),
            "linear_normalized_log2t",
            "Series matrix values centered ~1.0 (linear); log2(x+1) applied in prep",
            "limma", list(name = "CRSwNP_vs_UTctrl", case = "case", control = "control"),
            extras = list(list(name = "CRSsNP_vs_UTctrl",
                               case_vals = rownames(ph)[ph$characteristics_ch1.1 == "crs type: CRS without nasal polyps (CRSsNP)"],
                               ctrl_vals = rownames(ph)[grp == "control"])),
            note = "NP samples are CRSwNP only; control = 6 UT healthy")
}

## ---------------- 5. GSE194282 (GPL17692, log2) ----------------
{
  gse <- "GSE194282"
  es <- readRDS(file.path(OUT, "meta", "matrix", paste0(gse, ".rds")))[[1]]
  mat <- exprs(es)
  ph <- read.csv(file.path(OUT, "meta", "GSE194282_GPL17692_pheno.csv"), row.names = 1)
  grp <- ifelse(grepl("^NP", ph$title), "case", "control")
  age <- ifelse(grepl("elderly", ph$title) & !grepl("non-elderly", ph$title), "elderly", "non-elderly")
  meta <- data.frame(group = grp, age = age, row.names = rownames(ph))
  mp <- read.delim(file.path(OUT, "meta", "gpl", "GPL17692_id2symbol.tsv"), header = FALSE)
  sym <- mp$V2[match(rownames(mat), mp$V1)]
  save_prep(gse, mat, meta, sym, "log2_normalized",
            "Series matrix log2 (RMA-like, med~4)",
            "limma", list(name = "NP_vs_ctrl", case = "case", control = "control"),
            extras = list(list(name = "elderly_vs_nonelderly",
                               case_vals = rownames(ph)[age == "elderly"],
                               ctrl_vals = rownames(ph)[age == "non-elderly"])),
            note = "control = uncinate tissue incl. tumor-adjacent normal; age-stratified design")
}

## ---------------- 7. GSE72713 (FPKM-like -> log2) ----------------
{
  gse <- "GSE72713"
  d <- read.delim(gzfile(file.path(SUP, gse, "GSE72713_RNA-seq_processed_data.txt.gz")),
                  check.names = FALSE)
  sc <- grep("Ctrl|ECRSwNP|NECRSwNP|neCRSwNP", names(d), value = TRUE)
  mat <- log2(as.matrix(d[, sc]) + 1); rownames(mat) <- d$`Gene name`; colnames(mat) <- sc
  ph <- read.csv(file.path(OUT, "meta", "GSE72713_GPL11154_pheno.csv"), row.names = 1)
  # map count cols to GSM via title keywords
  diag <- ph$characteristics_ch1.1
  case_s <- rownames(ph)[grepl("chronic rhinosinusitis with nasal polyps", diag, ignore.case = TRUE)]
  ctrl_s <- rownames(ph)[grepl("Normal control", diag)]
  # column names in processed file: Ctrl_1..3, ECRSwNP_1..3, NECRSwNP_1..3 (verify order below)
  cat("[GSE72713] count cols:", paste(sc, collapse = ","), "\n")
  grp <- ifelse(grepl("ECRSwNP|NECRSwNP", sc, ignore.case = TRUE), "case",
         ifelse(grepl("Ctrl", sc), "control", NA))
  meta <- data.frame(group = grp, row.names = sc)
  e_s <- sc[grepl("^ECRSwNP", sc)]; ne_s <- sc[grepl("^non-ECRSwNP", sc, ignore.case = TRUE)]
  save_prep(gse, mat, meta, toupper(rownames(mat)), "linear_normalized_log2t",
            "Deposited processed values are FPKM-like (linear); log2(x+1) applied",
            "limma", list(name = "CRSwNP_vs_ctrl", case = "case", control = "control"),
            extras = list(list(name = "e_vs_ne", case_vals = e_s, ctrl_vals = ne_s)),
            note = "tiny n (3/3/3); lncRNA+mRNA mixed file, gene-name column used")
}

## ---------------- 8. GSE179265 (counts) ----------------
{
  gse <- "GSE179265"
  d <- read.delim(gzfile(file.path(SUP, gse, "GSE179265_readcount_genename.txt.gz")),
                  check.names = FALSE)
  smp <- names(d)[2:25]
  mat <- as.matrix(d[, smp]); rownames(mat) <- d$geneID; colnames(mat) <- smp
  storage.mode(mat) <- "integer"
  ph <- read.csv(file.path(OUT, "meta", "GSE179265_GPL24676_pheno.csv"), row.names = 1)
  # map: Control1..7, Jpn_NP_1..9, US_NP_1..8 (titles in pheno)
  diag <- ph$characteristics_ch1
  meta <- data.frame(group = ifelse(grepl("Healthy control", diag), "control", "case"),
                     ethnicity = ifelse(grepl("Caucasian", diag), "US",
                                 ifelse(grepl("Japanese", diag), "Jpn", "ctrl")),
                     row.names = smp)
  jp <- smp[grepl("^Jpn", smp)]; us <- smp[grepl("^US", smp)]; ct <- smp[grepl("^Control", smp)]
  save_prep(gse, mat, meta, ens2sym(rownames(mat)), "raw_counts",
            "Deposited read counts (geneID = Ensembl)",
            "deseq2", list(name = "CRSwNP_vs_ctrl", case = "case", control = "control"),
            extras = list(list(name = "US_vs_Jpn_polyps", case_vals = us, ctrl_vals = jp)),
            note = "ethnicity contrast available (US vs Japanese polyps)")
}

## ---------------- 9. GSE189690 (counts, N-ERD vs non) ----------------
{
  gse <- "GSE189690"
  d <- read.delim(gzfile(file.path(SUP, gse, "GSE189690_gene_counts.txt.gz")),
                  check.names = FALSE, comment.char = "#")
  cc <- grep("\\.bam$", names(d), value = TRUE)
  sids <- sub("\\.bam$", "", cc)
  mat <- as.matrix(d[, cc]); rownames(mat) <- d$Geneid; colnames(mat) <- sids
  storage.mode(mat) <- "integer"
  ph <- read.csv(file.path(OUT, "meta", "GSE189690_GPL20301_pheno.csv"), row.names = 1)
  pheno_id <- sub(".*, ", "", ph$title)  # "N-ERD, R8337" -> "R8337"
  pheno_grp <- sub(",.*", "", ph$title)
  m2 <- setNames(pheno_grp, pheno_id)
  pg <- m2[sids]; stopifnot(!any(is.na(pg)))
  meta <- data.frame(group = ifelse(pg == "N-ERD", "case", "control"), phenotype = pg, row.names = sids)
  save_prep(gse, mat, meta, ens2sym(rownames(mat)), "raw_counts",
            "Deposited featureCounts gene counts",
            "deseq2", list(name = "NERD_vs_nonNERD", case = "case", control = "control"),
            note = "NO healthy control: disease-subtype contrast only (N-ERD vs non-N-ERD CRSwNP)")
}

## ---------------- 10. GSE267331 (paired, log2) ----------------
{
  gse <- "GSE267331"
  es <- readRDS(file.path(OUT, "meta", "matrix", paste0(gse, ".rds")))[[1]]
  mat <- exprs(es)
  ph <- read.csv(file.path(OUT, "meta", "GSE267331_GPL570_pheno.csv"), row.names = 1)
  subj <- sub("subject id: ", "", ph$characteristics_ch1.1)
  time <- sub("time: ", "", ph$characteristics_ch1.2)
  trt <- sub("treatment: ", "", ph$characteristics_ch1)
  meta <- data.frame(subject = subj, time = time, treatment = trt, row.names = rownames(ph))
  # paired subjects with pre+post per arm
  dp_subj <- unique(subj[trt == "Dupilumab" & time %in% c("Pre-Treatment", "Post-Treatment")])
  dp_subj <- dp_subj[sapply(dp_subj, function(s) all(c("Pre-Treatment", "Post-Treatment") %in% time[subj == s & trt == "Dupilumab"]))]
  pl_subj <- unique(subj[trt == "Placebo" & time %in% c("Pre-Treatment", "Post-Treatment")])
  pl_subj <- pl_subj[sapply(pl_subj, function(s) all(c("Pre-Treatment", "Post-Treatment") %in% time[subj == s & trt == "Placebo"]))]
  cat(sprintf("[GSE267331] dupilumab pairs=%d placebo pairs=%d\n", length(dp_subj), length(pl_subj)))
  keep <- (trt == "Dupilumab" & subj %in% dp_subj) | (trt == "Placebo" & subj %in% pl_subj)
  mat <- mat[, keep]; meta <- meta[keep, ]
  meta$group <- ifelse(meta$time == "Post-Treatment", "case", "control")  # placeholder; paired handled via extras
  save_prep(gse, mat, meta, affy_sym(rownames(mat), "hgu133plus2.db"), "log2_normalized",
            "Series matrix RMA log2 (floor 2.31)",
            "limma_paired", list(name = "dupilumab_post_vs_pre", paired_arm = "Dupilumab"),
            extras = list(list(name = "placebo_post_vs_pre", paired_arm = "Placebo")),
            note = "SINUS-52 nasal brushings; paired within-subject; Early Termination excluded")
}

## ---------------- 11. GSE69093 (log2, celltype covariate) ----------------
{
  gse <- "GSE69093"
  es <- readRDS(file.path(OUT, "meta", "matrix", paste0(gse, ".rds")))[[1]]
  mat <- exprs(es)
  ph <- read.csv(file.path(OUT, "meta", "GSE69093_GPL10558_pheno.csv"), row.names = 1)
  grp <- ifelse(ph$characteristics_ch1 == "sample group: case", "case", "control")
  ct <- ifelse(grepl("Epithelial", ph$characteristics_ch1.1), "epi", "fib")
  meta <- data.frame(group = grp, celltype = ct, row.names = rownames(ph))
  save_prep(gse, mat, meta, ilmn_sym(rownames(mat)),
            "log2_normalized", "Series matrix log2; ID_REF are ILMN_ (HumanHT-12 v4) not Affy: platform misregistered as GPL10558",
            "limma_cov", list(name = "case_vs_ctrl_celltypeAdj", case = "case", control = "control"),
            extras = list(
              list(name = "epi_case_vs_ctrl", case_vals = rownames(ph)[grp == "case" & ct == "epi"],
                   ctrl_vals = rownames(ph)[grp == "control" & ct == "epi"]),
              list(name = "fib_case_vs_ctrl", case_vals = rownames(ph)[grp == "case" & ct == "fib"],
                   ctrl_vals = rownames(ph)[grp == "control" & ct == "fib"])),
            note = "primary nasal epithelial+fibroblast cells (culture status per series: primary component cells)")
}

## ---------------- 12. GSE313186 (counts, Exposed vs Control from title) ----------------
{
  gse <- "GSE313186"
  d <- read.delim(gzfile(file.path(SUP, gse, "GSE313186_Clinical_nasal_BPC_raw_featurecounts_matrix.txt.gz")),
                  check.names = FALSE)
  smp <- names(d)[8:27]
  mat <- as.matrix(d[, smp]); rownames(mat) <- d$GeneID; colnames(mat) <- smp
  storage.mode(mat) <- "integer"
  ph <- read.csv(file.path(OUT, "meta", "GSE313186_GPL34281_pheno.csv"), row.names = 1)
  # title: "Ethmoid sinus tissue, Exposed, BP001_S49_L008" -> key BP001 / BP_16 / UC_10 / UCSD0001
  tkey <- sub("_S[0-9]+_L[0-9]+$", "", sub("^.*, (Exposed|Control), ", "", ph$title))
  norm <- function(x) {
    x <- sub("^BP_0?", "BP", x); x <- sub("^UC_0?", "UCSD", x); x <- sub("^UCSD0*(?=[0-9])", "UCSD", x, perl = TRUE)
    x <- sub("^(BP|UCSD)([0-9]+)$", "\\1", x)
    x
  }
  # normalize to zero-padded form matching count columns
  norm2 <- function(x) {
    pre <- sub("[0-9]+$", "", x)
    dig <- regmatches(x, regexpr("[0-9]+$", x))
    dig <- sprintf("%03d", as.integer(dig))
    pre2 <- ifelse(pre %in% c("BP", "BP_"), "BP", ifelse(pre %in% c("UC_", "UCSD"), "UCSD", pre))
    paste0(pre2, dig)
  }
  tk <- norm2(tkey)
  col2gsm <- sapply(smp, function(s) { hit <- which(tk == s); if (!length(hit)) NA else rownames(ph)[hit[1]] })
  grp <- ifelse(grepl("Exposed", ph$title), "case", "control"); names(grp) <- rownames(ph)
  gg <- grp[col2gsm]; stopifnot(!any(is.na(gg)))
  meta <- data.frame(group = gg, row.names = smp)
  save_prep(gse, mat, meta, toupper(rownames(mat)), "raw_counts",
            "Deposited raw featureCounts matrix (gene symbols)",
            "deseq2", list(name = "exposed_vs_ctrl", case = "case", control = "control"),
            note = "deployment-exposed CRS vs unexposed controls; count cols mapped via title keys")
}

## ---------------- 13. GSE70898 + GSE70900 merged (linear -> log2) ----------------
{
  gse <- "GSE70898_70900"
  e1 <- readRDS(file.path(OUT, "meta", "matrix", "GSE70898.rds"))[[1]]
  e2 <- readRDS(file.path(OUT, "meta", "matrix", "GSE70900.rds"))[[1]]
  common <- intersect(rownames(e1), rownames(e2))
  mat <- log2(cbind(exprs(e1)[common, ], exprs(e2)[common, ]) + 1)
  p1 <- read.csv(file.path(OUT, "meta", "GSE70898_GPL10558_pheno.csv"), row.names = 1)
  p2 <- read.csv(file.path(OUT, "meta", "GSE70900_GPL10558_pheno.csv"), row.names = 1)
  g1 <- rep("case", nrow(p1))  # all polyp T cells
  g2 <- ifelse(grepl("Normal nasal", p2$title), "control",
        ifelse(grepl("Blood", p2$title), "blood_ref", NA))
  meta <- data.frame(group = c(g1, g2), row.names = c(rownames(p1), rownames(p2)))
  colnames(mat) <- rownames(meta)
  save_prep(gse, mat, meta, ilmn_sym(rownames(mat)),
            "linear_normalized_log2t",
            "Both series matrices linear-scale (med~100); merged same-study same-platform; log2(x+1)",
            "limma", list(name = "polypT_vs_nasalT", case = "case", control = "control"),
            note = "merged GSE70898 (cases) + GSE70900 (controls); blood T kept as reference arm only; short-term explant culture")
}

## ---------------- 14. GSE174681 (TMM linear -> log2) ----------------
{
  gse <- "GSE174681"
  d <- read.csv(gzfile(file.path(SUP, gse, "GSE174681_RNAseq_TMMnormalized.csv.gz")),
                check.names = FALSE)
  smp <- names(d)[2:11]
  mat <- log2(as.matrix(d[, smp]) + 1); rownames(mat) <- d$ENSEMBLE_ID; colnames(mat) <- smp
  tt <- smp  # TMM columns are the sample titles (NP_ASC_1789 etc.)
  meta <- data.frame(group = ifelse(grepl("^NP_ASC", tt), "case",
                             ifelse(grepl("^PBL_NAIVE", tt), "control", "other")),
                     subset = tt, row.names = smp)
  save_prep(gse, mat, meta, ens2sym(rownames(mat)), "linear_normalized_log2t",
            "Deposited TMM-normalized values (linear); log2(x+1) applied",
            "limma", list(name = "NP_ASC_vs_PBLnaive", case = "case", control = "control"),
            extras = list(list(name = "NP_naive_vs_PBLnaive",
                               case_vals = colnames(mat)[grepl("^NP_NAIVE", tt)],
                               ctrl_vals = colnames(mat)[grepl("^PBL_NAIVE", tt)])),
            note = "sorted B-cell subsets; atlas-type cohort, not a case-control cohort")
}

## ---------------- 15. GSE150846 (log2 matrix) ----------------
{
  gse <- "GSE150846"
  es <- readRDS(file.path(OUT, "meta", "matrix", paste0(gse, ".rds")))[[1]]
  mat <- exprs(es)
  ph <- read.csv(file.path(OUT, "meta", "GSE150846_GPL23126_pheno.csv"), row.names = 1)
  src <- ph$source_name_ch1
  meta <- data.frame(group = ifelse(grepl("Nasal polyp", src), "case",
                             ifelse(grepl("stimulated", src), "stim", "control")),
                     src = src, row.names = rownames(ph))
  save_prep(gse, mat, meta, affy_sym(rownames(mat), "clariomshumantranscriptcluster.db"),
            "log2_normalized", "Series matrix log2 (RMA-like)",
            "limma", list(name = "NP_ILC2_vs_blood_ILC2", case = "case", control = "control"),
            extras = list(list(name = "type2stim_vs_baseline",
                               case_vals = rownames(ph)[meta$group == "stim"],
                               ctrl_vals = rownames(ph)[meta$group == "control"])),
            note = "sorted ILC2 (3/3/3); LOW_POWER-adjacent tiny n")
}

## ---------------- 16. GSE330761 (txi counts -> voom) ----------------
{
  gse <- "GSE330761"
  e <- new.env(); load("/tmp/mctc_crs", envir = e)
  txi <- get("txi.CRS", envir = e)
  mat <- txi$counts
  smp <- colnames(mat)
  meta <- data.frame(group = ifelse(grepl("CRSwNP", smp), "case", "control"), row.names = smp)
  save_prep(gse, mat, meta, toupper(rownames(mat)), "tximport_counts",
            "Deposited tximport counts (non-integer estimates); gene symbols",
            "limma_voom", list(name = "CRSwNP_MCTC_vs_CRSsNP_MCTC", case = "case", control = "control"),
            note = "flow-sorted MCTC mast cells; no healthy control (CRSsNP as reference)")
}

## ---------------- 17. GSE213809 (per-sample counts in tar) ----------------
{
  gse <- "GSE213809"
  td <- file.path(SUP, gse, "x")
  unlink(td, recursive = TRUE)
  system(sprintf("mkdir -p %s && tar -xf %s -C %s", td,
                 file.path(SUP, gse, "GSE213809_RAW.tar"), td))
  fs <- list.files(td, full.names = TRUE)
  gsms <- sub("_.*$", "", basename(fs))
  gl <- lapply(fs, function(f) {
    d <- read.delim(gzfile(f), header = FALSE, check.names = FALSE)
    setNames(d$V2, d$V1)
  })
  genes <- Reduce(intersect, lapply(gl, names))
  mat <- do.call(cbind, lapply(gl, function(v) v[genes]))
  colnames(mat) <- gsms; storage.mode(mat) <- "integer"
  ph <- read.csv(file.path(OUT, "meta", "GSE213809_GPL15520_pheno.csv"), row.names = 1)
  tt <- ph[gsms, "title"]
  stim <- !grepl("before", tt)
  tissue <- ifelse(grepl(", NP,", tt), "NP", "PB")
  meta <- data.frame(group = ifelse(tissue == "NP" & stim, "case",
                             ifelse(tissue == "PB" & stim, "control", "other")),
                     tissue = tissue, stimulated = stim, title = tt, row.names = gsms)
  save_prep(gse, mat, meta, toupper(genes), "raw_counts",
            "Deposited per-sample count files (GSM*_countdata.tsv.gz in RAW.tar)",
            "deseq2", list(name = "NP_vs_PB_stimulated", case = "case", control = "control"),
            extras = list(list(name = "EP_vs_EN_NPstim",
                               case_vals = gsms[tissue == "NP" & stim & grepl(", EP", tt)],
                               ctrl_vals = gsms[tissue == "NP" & stim & grepl(", EN", tt)])),
            note = "sorted ILC2 +/- IL-33; tissue-of-origin x secretion response; donor #153 PB missing")
}
cat("ALL PREP DONE\n")
