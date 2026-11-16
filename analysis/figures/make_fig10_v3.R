# Fig 10 (v3 redraw): Per-sample program scores (ssGSEA).
# Conclusion: case-minus-control effects for 13 pre-defined programs across
# the 16 bulk cohorts; the type-2 program is up in five disease cohorts and
# down after dupilumab (GSE267331, FDR 2.4e-05). Asterisks: Wilcoxon BH-FDR.
# Heatmap with a single top annotation bar (contrast design), per v3 V12.
# Matrix construction replicates figE2_programs_fig10corrected.R exactly;
# palette, labels and framing follow v3 (row labels humanised, display only).
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(pheatmap); library(grid)})

d <- read.csv("analysis/step8_expansion/programs/program_summary.csv", stringsAsFactors = FALSE)
d$diff <- as.numeric(d$diff); d$FDR <- as.numeric(d$FDR)
ord_c <- c("GSE136825", "GSE179265", "GSE338426", "GSE36830", "GSE194282", "GSE72713",
           "GSE267331", "GSE198950", "GSE313186", "GSE69093",
           "GSE189690", "GSE330761", "GSE150846", "GSE174681", "GSE213809", "GSE70898_70900")
ord_p <- c("type2_custom", "il13_response", "eosinophil", "mast", "ilc2", "goblet_secretory",
           "fibroblast_ecm", "neutrophil", "interferon", "hypoxia", "ciliated", "barrier",
           "proliferation")
M <- matrix(NA, nrow = length(ord_p), ncol = length(ord_c),
            dimnames = list(ord_p, ord_c))
Sig <- M
for (i in seq_len(nrow(d))) {
  if (!d$cohort[i] %in% ord_c || !d$program[i] %in% ord_p) next
  M[d$program[i], d$cohort[i]] <- d$diff[i]
  Sig[d$program[i], d$cohort[i]] <- d$FDR[i] < 0.05
}
stopifnot(all(dim(M) == c(13, 16)))
fdr_t2_dupi <- d$FDR[d$cohort == "GSE267331" & d$program == "type2_custom"]
stopifnot(length(fdr_t2_dupi) == 1, round(-log10(fdr_t2_dupi), 1) == 4.6)  # FDR 2.4e-05
lim <- max(abs(M), na.rm = TRUE)

## display-only label cleaning (V8): program names humanised
prog_lab <- c(type2_custom = "Type 2", il13_response = "IL-13 response",
              eosinophil = "Eosinophil", mast = "Mast cell", ilc2 = "ILC2",
              goblet_secretory = "Goblet / secretory",
              fibroblast_ecm = "Fibroblast / ECM", neutrophil = "Neutrophil",
              interferon = "Interferon", hypoxia = "Hypoxia",
              ciliated = "Ciliated", barrier = "Barrier",
              proliferation = "Proliferation")
rownames(M) <- prog_lab[rownames(M)]
rownames(Sig) <- rownames(M)
ann <- data.frame(Design = c(rep("Disease vs control", 6), "Treatment",
                             "Tissue-confounded", "Exposure", "Weak signal",
                             rep("Subtype / sorted", 6)),
                  row.names = ord_c)
ann_cols <- list(Design = c("Disease vs control" = COL$disease2,
                            "Treatment" = COL$disease1,
                            "Tissue-confounded" = COL$neutral,
                            "Exposure" = COL$accent,
                            "Weak signal" = "#F39B7F",
                            "Subtype / sorted" = "#00A087"))
draw <- function() {
  ## pheatmap draws its own full page; capture it and redraw shifted 1% right
  ## so column labels never touch the left device edge (as the corrected
  ## submitted render did). Pure framing: no data/palette change.
  pheatmap(M, color = colorRampPalette(c(COL$down, "white", COL$up))(100),
           breaks = seq(-lim, lim, length.out = 101),
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = matrix(ifelse(Sig, "*", ""), nrow(Sig)),
           number_color = "black", fontsize_number = 10,
           annotation_col = ann, annotation_colors = ann_cols,
           annotation_names_col = FALSE,
           fontsize = 8, cellwidth = 17, cellheight = 22, angle_col = 45,
           gaps_col = c(6, 10), silent = FALSE)
  g <- grid.grab()
  grid.newpage()
  pushViewport(viewport(x = 0.51, y = 0.5, width = 0.98, height = 1, just = "center"))
  grid.draw(g)
}
w <- 180 / 25.4; h <- 135 / 25.4
grDevices::cairo_pdf("manuscript/figures/Fig10.pdf", width = w, height = h, family = "Arial")
draw(); dev.off()
grDevices::png("manuscript/figures/Fig10.png", width = w, height = h, units = "in", res = 300)
draw(); dev.off()
cat("Fig10 v3 done. NA cells:", sum(is.na(M)), "| type-2 dupilumab FDR:", fdr_t2_dupi, "\n")
