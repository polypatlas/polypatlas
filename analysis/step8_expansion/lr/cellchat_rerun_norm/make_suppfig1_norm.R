# Supplementary Fig 1 (v3 redraw, NORMALIZED rerun 2026-11-16): same layout as make_suppfig1_v3.R, data = CellChat rerun with normalizeData (see cellchat_run_norm.R).
# Conclusion: top gained interactions in CRSwNP are a tuft-cell-centred
# innate axis (APP-CD74, ANXA1-FPR1, MIF-CD74/CXCR4); classical type-2
# cytokine pairs are below reliable detection (caption).
# Chart type: diverging lollipop panels (ranking -> lollipop, v3 table;
# replaces the submitted base-R bars). Interaction selection replicates
# figE3_LR.R exactly (first 12 positive / negative rows in file order,
# which equal the 12 largest gains / losses - asserted).
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(patchwork)})

d <- read.csv("analysis/step8_expansion/lr/cellchat_rerun_norm/lr_interactions_CRSwNP_vs_control_norm.csv",
              stringsAsFactors = FALSE)
top_gain <- head(d[d$delta > 0, ], 12)
top_loss <- head(d[d$delta < 0, ], 12)
stopifnot(nrow(top_gain) == 12, nrow(top_loss) == 12,
          isTRUE(all.equal(sort(round(top_gain$delta, 6), decreasing = TRUE),
                           sort(round(d$delta[d$delta > 0], 6), decreasing = TRUE)[1:12])),
          isTRUE(all.equal(sort(round(top_loss$delta, 6)),
                           sort(round(d$delta[d$delta < 0], 6))[1:12])))
## display-only label cleaning (V8): arrows spaced; underscores in cell-type
## names -> spaces. Gene symbols (incl. the CD74_CXCR4 complex) are proper
## names and stay unchanged (registered exception).
mk_lab <- function(df) {
  paste0(df$ligand, " → ", df$receptor, "\n",
         gsub("_", " ", df$source), " → ", gsub("_", " ", df$target))
}
top_gain$label <- mk_lab(top_gain); top_loss$label <- mk_lab(top_loss)
top_gain$panel <- "Gained in CRSwNP · top 12"
top_loss$panel <- "Lost in CRSwNP · top 12"
both <- rbind(top_gain, top_loss)
both <- both[order(both$panel, both$delta), ]
both$label <- factor(both$label, levels = both$label)
both$panel <- factor(both$panel, levels = c("Gained in CRSwNP · top 12",
                                            "Lost in CRSwNP · top 12"))

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))
p <- ggplot(both, aes(x = delta, y = label, colour = panel)) +
  geom_vline(xintercept = 0, colour = COL$grey_ax, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = delta, y = label, yend = label),
               colour = "#CCCCCC", linewidth = 0.5) +
  geom_point(size = 2.4) +
  facet_wrap(~panel, ncol = 1, scales = "free") +
  scale_colour_manual(values = c("Gained in CRSwNP · top 12" = COL$disease1,
                                 "Lost in CRSwNP · top 12" = COL$disease2),
                      guide = "none") +
  scale_x_continuous(labels = function(v) format(v, scientific = TRUE, digits = 1)) +
  labs(x = "Delta probability (CRSwNP − control)", y = NULL) + th +
  theme(axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.3))
save_fig(p, "analysis/step8_expansion/lr/cellchat_rerun_norm/SuppFig1_norm", width_mm = 180, height_mm = 175)
cat("SuppFig1_norm done. 12 gained + 12 lost interactions, selection as submitted\n")
