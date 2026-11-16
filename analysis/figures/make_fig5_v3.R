# Fig 5 (v3 redraw): Single-cell composition, flagship dataset GSE235711.
# Conclusions: (a) donor-level composition of 18 unified classes across
# 10 CRSwNP / 5 CRSsNP / 5 control donors; (b) goblet, secretory, mast and
# ciliated compartments shift most across groups (Kruskal-Wallis).
# Chart types: 100% stacked composition bars (composition exception) and
# violin + per-donor points (distribution comparison, v3 selection table).
# Data and KW computation identical to the original script.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(ggbeeswarm); library(patchwork)})

comp <- read.csv("analysis/step3_scrna/composition_all.csv", stringsAsFactors = FALSE)
d <- comp[comp$dataset == "GSE235711", ]
stopifnot(length(unique(d$donor)) == 20, length(unique(d$celltype)) == 18)
d$group <- factor(d$group, levels = c("control", "CRSsNP", "CRSwNP"))

## donor order within group (order of appearance, as original)
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
  scale_fill_manual(values = pal_cat(18), name = "Cell type",
                    labels = function(x) ifelse(x %in% names(ct_lab), ct_lab[x], x),
                    guide = guide_legend(ncol = 2)) +
  labs(x = "Donor (n = 20)", y = "Fraction of cells") + th +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5))

## Kruskal-Wallis per cell type across the 3 groups (as original)
kw <- do.call(rbind, lapply(unique(d$celltype), function(ct) {
  x <- d[d$celltype == ct, ]
  data.frame(celltype = ct,
             p = tryCatch(kruskal.test(fraction ~ group, data = x)$p.value,
                          error = function(e) NA))
}))
kw <- kw[order(kw$p), ]
top4 <- kw$celltype[1:4]
stopifnot(setequal(top4, c("goblet", "secretory", "mast", "ciliated")))
pv <- setNames(kw$p[1:4], top4)
stopifnot(abs(pv[["goblet"]] - 0.0030) < 5e-4, abs(pv[["secretory"]] - 0.0030) < 5e-4,
          abs(pv[["mast"]] - 0.0044) < 5e-4, abs(pv[["ciliated"]] - 0.014) < 5e-4)
d2 <- d[d$celltype %in% top4, ]
d2$celltype <- factor(d2$celltype, levels = top4)
## strip labels formatted to match the caption exactly (0.0030/0.0030/0.0044/0.014)
fmt_p <- function(p) ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.3f", p))
stopifnot(identical(unname(fmt_p(pv[c("goblet", "secretory", "mast", "ciliated")])),
                    c("0.0030", "0.0030", "0.0044", "0.014")))
strip_lab <- setNames(sprintf("%s\np = %s", top4, fmt_p(kw$p[1:4])), top4)
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
cat("Fig5 v3 done. Top KW shifts:\n"); print(kw[1:4, ])
