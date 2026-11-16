# Fig 3 (v3 redraw): Bulk differential-expression statistics.
# Conclusions: (a) DEG yield varies by orders of magnitude across contrasts;
# (b) yield does not track per-arm sample size (r = -0.11); (c-d) the two
# flagship cohorts show genome-wide, direction-balanced signal.
# Chart types: stacked direction bars (composition exception), scatter + fit,
# two volcano scatters. Data prep identical to the original script (same DE
# files, same runlog mapping); only presentation follows v3.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(ggrepel); library(patchwork)})

runlog <- read.csv("analysis/step2_bulk/runlog.csv", stringsAsFactors = FALSE)
de_files <- list.files("analysis/step2_bulk/de", pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(de_files) == 30, nrow(runlog) == 30)
file_for <- function(gse, contrast) {
  stems <- sub("\\.csv$", "", sub(paste0("^", gse, "_"), "", basename(de_files)))
  hit <- de_files[tolower(stems) == tolower(contrast)]
  if (length(hit) != 1) stop(paste("ambiguous/missing DE file for", gse, contrast))
  hit
}
deg <- do.call(rbind, lapply(seq_len(nrow(runlog)), function(i) {
  gse <- runlog$gse[i]; ct <- runlog$contrast[i]
  d <- read.csv(file_for(gse, ct), stringsAsFactors = FALSE)
  data.frame(gse = gse, contrast = ct,
             n_case = runlog$n_case[i], n_ctrl = runlog$n_control[i],
             up = sum(d$FDR < 0.05 & d$logFC > 0, na.rm = TRUE),
             down = sum(d$FDR < 0.05 & d$logFC < 0, na.rm = TRUE),
             total = sum(d$FDR < 0.05, na.rm = TRUE))
}))
deg <- deg[order(deg$total, decreasing = TRUE), ]
deg$harm_n <- 2 * deg$n_case * deg$n_ctrl / (deg$n_case + deg$n_ctrl)
stopifnot(deg$total[deg$gse == "GSE136825" & deg$contrast == "polyp_vs_ctrl"] == 5320,
          deg$total[deg$gse == "GSE179265" & deg$contrast == "CRSwNP_vs_ctrl"] == 5278)

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))

## (a) stacked up/down DEG bars (direction composition per contrast)
dA <- reshape(deg, varying = list(c("up", "down")), v.names = "n",
              timevar = "dir", times = c("Up", "Down"), direction = "long")
dA$dir <- factor(dA$dir, levels = c("Up", "Down"))
# display-only label cleaning (v3 V8): underscores in contrast names -> spaces
dA$short <- factor(paste0(dA$gse, " | ", gsub("_", " ", dA$contrast)),
                   levels = paste0(deg$gse, " | ", gsub("_", " ", deg$contrast)))
pA <- ggplot(dA, aes(x = short, y = n, fill = dir)) +
  geom_col(colour = NA) +
  coord_flip() +
  scale_fill_manual(values = c(Up = COL$up, Down = COL$down), name = "FDR < 0.05") +
  labs(y = "DEGs (FDR < 0.05)", x = NULL) + th +
  theme(axis.text.y = element_text(size = 7), legend.position = "top")

## (b) DEG yield vs harmonic-mean per-arm n
r <- cor(log10(deg$harm_n), log10(deg$total + 1), use = "complete.obs")
stopifnot(round(r, 2) == -0.11, nrow(deg) == 30)
pB <- ggplot(deg, aes(x = harm_n, y = total)) +
  geom_point(colour = COL$disease2, size = 2.2) +
  geom_smooth(method = "lm", se = FALSE, colour = COL$neutral,
              linewidth = 0.7, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x = min(deg$harm_n) * 1.05, y = max(deg$total) / 40, hjust = 0, vjust = 1,
           size = 2.7, label = sprintf("Pearson r (log10) = %.2f\nn = %d contrasts", r, nrow(deg))) +
  labs(x = "Harmonic mean samples per arm (log scale)",
       y = "Total DEGs at FDR < 0.05 (log scale)") + th

## (c-d) volcanoes, top 10 genes by FDR labelled
volc <- function(f, cohort_lab) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$mlp <- -log10(pmax(d$p, 1e-300))
  d$sig <- ifelse(d$FDR < 0.05, ifelse(d$logFC > 0, "Up", "Down"), "ns")
  lab <- d[order(d$FDR), ][1:10, ]
  ggplot(d, aes(x = logFC, y = mlp)) +
    geom_point(aes(colour = sig), size = 0.5, alpha = 0.6) +
    scale_colour_manual(values = c(ns = "#C9CED6", Up = COL$up, Down = COL$down),
                        guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = COL$neutral, linewidth = 0.5) +
    geom_text_repel(data = lab, aes(label = symbol), size = 2.5,
                    max.overlaps = 20, box.padding = 0.3) +
    annotate("text", x = min(d$logFC) * 0.98, y = -log10(0.05) + max(d$mlp) * 0.03, hjust = 0,
             vjust = 0, size = 3, fontface = "bold", label = cohort_lab) +
    labs(x = "logFC", y = "-log10 p-value") + th
}
pC <- volc("analysis/step2_bulk/de/GSE136825_polyp_vs_ctrl.csv", "GSE136825 · 75 vs 28")
pD <- volc("analysis/step2_bulk/de/GSE179265_CRSwNP_vs_ctrl.csv", "GSE179265 · 17 vs 7")

fig <- (pA / (pB | pC) / pD) +
  plot_layout(heights = c(1.35, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))
save_fig(fig, "manuscript/figures/Fig3", width_mm = 180, height_mm = 230)
cat(sprintf("Fig3 v3 done. Pearson r = %.2f; contrasts = %d\n", r, nrow(deg)))
