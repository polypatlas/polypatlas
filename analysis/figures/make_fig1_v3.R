# Fig 1 (v3 redraw): PolypAtlas full-pipeline flowchart.
# Approved composition retained (user-approved layout; nodes and counts unchanged).
# v3 normalisation only: locked palette, Arial, in-figure text >= 7 pt, alignment.
# User ruling 2026-11-09 (workflow figures: text only in boxes, no decorative
# icons) takes precedence over v3 V5 for this figure -> V5 declared N/A.
# Sources: curation/curation_table.csv (66 series); counts asserted below.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages(library(grid))

cur <- read.csv("curation/curation_table.csv", stringsAsFactors = FALSE)
tab <- table(factor(cur$verdict, levels = c("core", "pending", "excluded")))
stopifnot(nrow(cur) == 66, sum(cur$n_samples) == 2046,
          tab[["core"]] == 32, tab[["pending"]] == 9, tab[["excluded"]] == 25,
          sum(cur$n_samples[cur$verdict == "core"]) == 802,
          sum(cur$n_samples[cur$verdict == "pending"]) == 227,
          sum(cur$n_samples[cur$verdict == "excluded"]) == 1017)

cols <- list(
  census = COL$disease2, cur = COL$disease2,
  bulk = "#00A087", scrna = "#5B4E77", spat = COL$accent,
  ana = COL$disease1, val = COL$disease2, web = "#404040",
  edge = "#7F7F7F", tx = "white", txd = "#222222"
)

plot_fig <- function() {
  grid.newpage()
  pushViewport(viewport(w = 1, h = 1))
  box <- function(x, y, wd, ht, fill, l1, l2 = NULL, l3 = NULL, fs1 = 8.5, fs2 = 7) {
    grid.roundrect(x = x, y = y, width = wd, height = ht, r = unit(0.02, "npc"),
                   gp = gpar(fill = fill, col = "#333333", lwd = 1.2))
    yy <- y + ht * 0.24
    grid.text(l1, x = x, y = yy, gp = gpar(col = cols$tx, fontsize = fs1, fontface = "bold", family = "Arial"))
    if (!is.null(l2)) grid.text(l2, x = x, y = yy - ht * 0.24, gp = gpar(col = cols$tx, fontsize = fs2, family = "Arial"))
    if (!is.null(l3)) grid.text(l3, x = x, y = yy - ht * 0.44, gp = gpar(col = cols$tx, fontsize = fs2, family = "Arial"))
  }
  arrow <- function(x1, y1, x2, y2, curve = 0) {
    grid.curve(x1, y1, x2, y2, curvature = curve,
               arrow = grid::arrow(type = "closed", length = unit(0.02, "npc")),
               gp = gpar(col = cols$edge, lwd = 1.4, fill = cols$edge))
  }
  lab <- function(x, y, s, fs = 7) grid.text(s, x = x, y = y, gp = gpar(col = cols$txd, fontsize = fs, fontface = "italic", family = "Arial"))

  box(0.15, 0.899, 0.26, 0.10, cols$census, "Census (GEO, 2026-10-01)", "66 human series", "2,046 samples total")
  box(0.50, 0.899, 0.30, 0.10, cols$cur, "Per-series curation", "32 core / 9 pending / 25 excluded", "verdicts + evidence notes")
  arrow(0.28, 0.899, 0.35, 0.899)

  box(0.17, 0.681, 0.26, 0.134, cols$bulk, "Bulk reanalysis", "16 cohorts · 30 contrasts", "19,266 genes · DESeq2/limma")
  box(0.50, 0.681, 0.30, 0.134, cols$scrna, "Single-cell unified annotation", "10 datasets · 476,254 cells", "105 donor units · 41 pseudobulk contrasts")
  box(0.83, 0.681, 0.26, 0.134, cols$spat, "Spatial layer", "GSE235714 GeoMx · 40 ROIs", "2 compartments · descriptive scope")
  lab(0.17, 0.795, "uniform pipelines + scale checks")
  lab(0.50, 0.795, "positive-control assertions")
  arrow(0.35, 0.849, 0.17, 0.748, curve = 0.08)
  arrow(0.50, 0.849, 0.50, 0.748)
  arrow(0.65, 0.849, 0.83, 0.748, curve = -0.08)

  box(0.50, 0.496, 0.50, 0.118, cols$ana, "Analysis layer (new)", "enrichment 57 ORA sets + 62 GSEA · 13 program scores", "MCP-style deconvolution · CellChat LR")
  arrow(0.17, 0.614, 0.42, 0.555, curve = 0.06)
  arrow(0.50, 0.614, 0.50, 0.555)
  arrow(0.83, 0.614, 0.58, 0.555, curve = -0.06)

  box(0.50, 0.328, 0.58, 0.109, cols$val, "Validation", "marker recovery 87.5% · dupilumab 7/8 reversed", "audit: 23/60 reproduced · 53/60 direction-consistent · 825/196 consensus")
  arrow(0.50, 0.437, 0.50, 0.383)

  box(0.50, 0.210, 0.44, 0.092, cols$web, "PolypAtlas website (7 pages, bilingual)", "the live resource - one box per page below")
  arrow(0.50, 0.273, 0.50, 0.256)

  pages <- list(
    list(x = 0.0825, l1 = "Home"),
    list(x = 0.2192, l1 = "Gene", l2 = "query"),
    list(x = 0.3558, l1 = "Pathway", l2 = "query"),
    list(x = 0.4925, l1 = "Cohorts &", l2 = "datasets"),
    list(x = 0.6292, l1 = "Treatment", l2 = "response"),
    list(x = 0.7658, l1 = "Downloads"),
    list(x = 0.9025, l1 = "Methods &", l2 = "limitations")
  )
  bus_y <- 0.140; page_top <- 0.125; page_y <- 0.075; page_h <- 0.10; page_w <- 0.125
  grid.lines(x = c(0.50, 0.50), y = c(0.164, bus_y), gp = gpar(col = cols$edge, lwd = 1.4))
  grid.lines(x = c(pages[[1]]$x, pages[[7]]$x), y = c(bus_y, bus_y), gp = gpar(col = cols$edge, lwd = 1.4))
  for (p in pages) {
    arrow(p$x, bus_y, p$x, page_top + 0.002)
    box(p$x, page_y, page_w, page_h, cols$web, p$l1, p$l2, fs1 = 7, fs2 = 7)
  }
}

w <- 180 / 25.4; h <- 140 / 25.4
grDevices::cairo_pdf("manuscript/figures/Fig1.pdf", width = w, height = h, family = "Arial")
plot_fig(); dev.off()
grDevices::png("manuscript/figures/Fig1.png", width = w, height = h, units = "in", res = 300)
plot_fig(); dev.off()
cat("Fig1 v3 done\n")
