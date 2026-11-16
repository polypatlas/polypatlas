#!/usr/bin/env Rscript
# rebuild_runlog.R — regenerate runlog.csv rows for already-computed de/ files.
OUT <- "analysis/step2_bulk"
# contrast -> prep gse + method + n_case/n_ctrl (from prep objects)
info <- list(
  GSE338426_CRSwNP_vs_ctrl = c("GSE338426", "deseq2"),
  GSE338426_eCRSwNP_vs_ctrl = c("GSE338426", "deseq2"),
  GSE338426_neCRSwNP_vs_ctrl = c("GSE338426", "deseq2"),
  GSE198950_CRSwNP_vs_NonCRS = c("GSE198950", "deseq2"),
  GSE198950_CRSsNP_vs_NonCRS = c("GSE198950", "deseq2"),
  GSE198950_EOS_high_vs_low = c("GSE198950", "deseq2"),
  GSE136825_polyp_vs_ctrl = c("GSE136825", "deseq2"),
  GSE36830_CRSwNP_vs_UTctrl = c("GSE36830", "limma"),
  GSE36830_CRSsNP_vs_UTctrl = c("GSE36830", "limma"),
  GSE194282_NP_vs_ctrl = c("GSE194282", "limma"),
  GSE194282_elderly_vs_nonelderly = c("GSE194282", "limma"),
  GSE72713_CRSwNP_vs_ctrl = c("GSE72713", "limma"),
  GSE72713_e_vs_ne = c("GSE72713", "limma"),
  GSE179265_CRSwNP_vs_ctrl = c("GSE179265", "deseq2"),
  GSE179265_US_vs_Jpn_polyps = c("GSE179265", "deseq2"),
  GSE189690_NERD_vs_nonNERD = c("GSE189690", "deseq2"),
  GSE267331_dupilumab_post_vs_pre = c("GSE267331", "limma_paired")
)
# n_case/n_ctrl per contrast (from prep logs)
nn <- list(
  GSE338426_CRSwNP_vs_ctrl = c(14, 8), GSE338426_eCRSwNP_vs_ctrl = c(9, 8),
  GSE338426_neCRSwNP_vs_ctrl = c(5, 8),
  GSE198950_CRSwNP_vs_NonCRS = c(5, 4), GSE198950_CRSsNP_vs_NonCRS = c(11, 4),
  GSE198950_EOS_high_vs_low = c(8, 7),
  GSE136825_polyp_vs_ctrl = c(75, 28),
  GSE36830_CRSwNP_vs_UTctrl = c(12, 6), GSE36830_CRSsNP_vs_UTctrl = c(6, 6),
  GSE194282_NP_vs_ctrl = c(7, 7), GSE194282_elderly_vs_nonelderly = c(8, 6),
  GSE72713_CRSwNP_vs_ctrl = c(6, 3), GSE72713_e_vs_ne = c(3, 3),
  GSE179265_CRSwNP_vs_ctrl = c(17, 7), GSE179265_US_vs_Jpn_polyps = c(8, 9),
  GSE189690_NERD_vs_nonNERD = c(13, 13),
  GSE267331_dupilumab_post_vs_pre = c(47, 47)
)
scales <- c(GSE338426 = "raw_counts", GSE198950 = "raw_counts", GSE136825 = "raw_counts",
            GSE36830 = "linear_normalized_log2t", GSE194282 = "log2_normalized",
            GSE72713 = "linear_normalized_log2t", GSE179265 = "raw_counts",
            GSE189690 = "raw_counts", GSE267331 = "log2_normalized")
RL <- file.path(OUT, "runlog.csv")
writeLines("gse,contrast,n_case,n_control,scale_verdict,method,status,note", RL)
for (key in names(info)) {
  f <- file.path(OUT, "de", paste0(key, ".csv"))
  if (!file.exists(f)) { cat("missing:", key, "\n"); next }
  d <- read.csv(f)
  sig <- sum(d$FDR < 0.05, na.rm = TRUE)
  gse <- info[[key]][1]; method <- info[[key]][2]
  n_case <- nn[[key]][1]; n_ctrl <- nn[[key]][2]
  cname <- sub(paste0(gse, "_"), "", key)
  sp <- readRDS(file.path(OUT, "prep", paste0(gse, ".rds")))
  status <- if (n_case < 3 || n_ctrl < 3) "LOW_POWER" else "OK"
  line <- sprintf('%s,%s,%d,%d,%s,%s,%s,"%s (%d FDR<0.05)"', gse, cname, n_case, n_ctrl,
                  scales[[gse]], method, status, sp$note, sig)
  cat(line, "\n", file = RL, append = TRUE)
}
cat("runlog rebuilt\n")
