# Fig 5 (v1.8 redraw): Single-cell composition, flagship dataset GSE235711.
# Release v1.8 restores 'secretory_club' as a label distinct from 'secretory'
# in the composition table (the two had been collapsed into one name, which
# double-counted donors in the Kruskal-Wallis test of the earlier figure).
# Conclusions: (a) donor-level composition of 19 populated classes across
# 10 CRSwNP / 5 CRSsNP / 5 control donors; (b) goblet, mast, secretory_club
# and ciliated compartments shift most across groups (Kruskal-Wallis with
# Benjamini-Hochberg FDR across the 19 classes; full table in
# analysis/step3_scrna/fig5_composition_stats_v1.8.csv).
# Chart types: 100% stacked composition bars (composition exception) and
# violin + per-donor points (distribution comparison, v3 selection table).
source("~/workspace/custom-skills/figure_style_v3/figstyle_v3.R")
suppressPackageStartupMessages({library(ggbeeswarm); library(patchwork)})

comp <- read.csv("analysis/step3_scrna/composition_all.csv", stringsAsFactors = FALSE)
d <- comp[comp$dataset == "GSE235711", ]
stopifnot(length(unique(d$donor)) == 20, length(unique(d$celltype)) == 19)
d$group <- factor(d$group, levels = c("control", "CRSsNP", "CRSwNP"))

## donor order within group (order of appearance)
ord <- do.call(rbind, lapply(levels(d$group), function(g) {
  dd <- unique(d[d$group == g, c("donor", "group")])
  data.frame(donor = dd$donor, group = g, stringsAsFactors = FALSE)
}))
stopifnot(table(ord$group) == c(control = 5, CRSsNP = 5, CRSwNP = 10))
d$donor <- factor(d$donor, levels = ord$donor)

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7),
                            legend.title = element_text(size = 8))
ct_lab <- c(CD4_T = "CD4 T", CD8_T = "CD8 T")   # display-only label cleaning (V8)
pA <- ggplot(d, aes(x = donor, y = fraction, fill = celltype)) +
  geom_col(colour = "grey40", linewidth = 0.15) +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_cat(19), name = "Cell type",
                    labels = function(x) ifelse(x %in% names(ct_lab), ct_lab[x], x),
                    guide = guide_legend(ncol = 2)) +
  labs(x = "Donor (n = 20)", y = "Fraction of cells") + th +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5))

## Kruskal-Wallis per class + BH across the 19 classes (matches the stats file)
st <- read.csv("analysis/step3_scrna/fig5_composition_stats_v1.8.csv",
               stringsAsFactors = FALSE)
top4 <- st$celltype[1:4]
stopifnot(setequal(top4, c("goblet", "mast", "secretory_club", "ciliated")))
pv <- setNames(st$kw_p[1:4], top4); qv <- setNames(st$bh_fdr[1:4], top4)
stopifnot(abs(pv[["goblet"]] - 0.0030) < 5e-4, abs(pv[["mast"]] - 0.0044) < 5e-4,
          abs(pv[["secretory_club"]] - 0.0082) < 5e-4, abs(pv[["ciliated"]] - 0.0136) < 5e-4)
d2 <- d[d$celltype %in% top4, ]
d2$celltype <- factor(d2$celltype, levels = top4)
fmt_p <- function(p) ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.3f", p))
strip_lab <- setNames(sprintf("%s\np = %s\nFDR = %.3f", top4, fmt_p(pv), qv), top4)
grp_cols <- c(control = COL$neutral, CRSsNP = COL$disease2, CRSwNP = COL$disease1)
pB <- ggplot(d2, aes(x = group, y = fraction, fill = group)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.75, scale = "width") +
  geom_quasirandom(width = 0.18, size = 1.1, colour = "#222222", alpha = 0.75) +
  stat_summary(fun = median, geom = "point", size = 2.2, colour = "white",
               shape = 21, fill = "black", stroke = 0.6) +
  facet_wrap(~celltype, ncol = 4, labeller = labeller(celltype = strip_lab)) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  coord_cartesian(ylim = c(0, 0.36)) +
  labs(x = NULL, y = "Fraction of cells") + th +
  theme(strip.text = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5))

fig <- (pA / pB) + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))
save_fig(fig, "manuscript/figures/Fig5", width_mm = 180, height_mm = 200)
cat("Fig5 v1.8 done. Top KW shifts:\n"); print(st[1:4, c("celltype", "kw_p", "bh_fdr")])
