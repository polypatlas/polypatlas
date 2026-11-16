# Fig 7 (v3 redraw): Consensus genes, top 15 up and top 15 down.
# Conclusion: the strongest disease-only consensus genes are direction-
# consistent in 80-100% of the 10 disease cohorts, with FDR support in 4-6.
# Chart type: diverging lollipop (ranking -> lollipop per v3 selection table;
# replaces the submitted horizontal bars). Gene selection and ordering
# replicate analysis/step6_validation/step6_figures.R exactly (top 15 per
# direction by n_sig_FDR005, then n_consistent, in consensus_genes.csv
# version "disease_cohorts_10"); gene sets are asserted against the
# submitted figure.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(dplyr); library(patchwork)})

cg <- read.csv("analysis/step6_validation/consensus_genes.csv", stringsAsFactors = FALSE)
cb <- cg %>% filter(version == "disease_cohorts_10") %>% group_by(direction) %>%
  arrange(desc(n_sig_FDR005), desc(n_consistent)) %>% slice_head(n = 15) %>%
  ungroup()
up_genes <- c("ATP8A2", "ARRB2", "SH2D7", "PNPLA6", "NTRK1", "JSRP1", "HMOX1",
              "SPAG4", "SERPINB4", "RGS19", "P2RY6", "CDH22", "A4GALT",
              "EFHD2", "CDH26")
dn_genes <- c("CD300LG", "CCBE1", "CA3", "PON3", "NELL2", "NDNF", "GPM6A",
              "BMP4", "BCKDHB", "PDP2", "LANCL1", "SEMA3E", "DLG2", "AR",
              "ALDH6A1")
stopifnot(setequal(cb$gene[cb$direction == "up"], up_genes),
          setequal(cb$gene[cb$direction == "down"], dn_genes),
          nrow(cb) == 30)
cb$x <- ifelse(cb$direction == "up", cb$frac_consistent, -cb$frac_consistent)
cb$lab <- paste0(cb$n_sig_FDR005, " sig · ", cb$n_consistent, "/",
                 cb$n_cohorts_measured)

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))
mk <- function(dd, col) {
  d <- cb[cb$direction == dd, ]
  d <- d[order(d$frac_consistent, d$n_sig_FDR005), ]
  d$gene <- factor(d$gene, levels = d$gene)
  ggplot(d, aes(x = x, y = gene)) +
    geom_vline(xintercept = 0, colour = COL$grey_ax, linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = x, y = gene, yend = gene),
                 colour = "#CCCCCC", linewidth = 0.5) +
    geom_point(size = 2.6, colour = col) +
    geom_text(aes(x = x + ifelse(dd == "up", 0.035, -0.035), label = lab),
              hjust = ifelse(dd == "up", 0, 1),
              size = 2.5, colour = COL$grey_txt) +
    scale_x_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                       labels = function(v) abs(v),
                       limits = if (dd == "up") c(0, 1.52) else c(-1.52, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0.02, 0.15))) +
    labs(x = NULL, y = NULL) + th +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.3))
}
pU <- mk("up", COL$up) +
  annotate("text", x = 0.02, y = 16.1, hjust = 0, size = 3.2, fontface = "bold",
           label = "Up · 825 genes")
pD <- mk("down", COL$down) +
  labs(x = "Fraction of cohorts direction-consistent") +
  annotate("text", x = -0.02, y = 16.1, hjust = 1, size = 3.2, fontface = "bold",
           label = "Down · 196 genes")
fig <- (pU / pD) + plot_layout(heights = c(1, 1))
save_fig(fig, "manuscript/figures/Fig7", width_mm = 180, height_mm = 175)
cat("Fig7 v3 done. 15 up + 15 down genes, sets match the submitted figure\n")
