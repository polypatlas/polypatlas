#!/usr/bin/env Rscript
# Fig. 5 composition statistics, release v1.8 (GSE235711, donor-level).
# For each of the 19 populated composition classes of the flagship dataset:
# Kruskal-Wallis test across control / CRSsNP / CRSwNP donor fractions
# (n = 5 / 5 / 10 donors), Benjamini-Hochberg FDR across the 19 classes,
# and Hodges-Lehmann location estimates with 95% CIs (wilcox.test) for the
# two disease-vs-control comparisons. Input: the restored composition table
# (composition_all.csv, release v1.8; secretory and secretory_club are
# distinct labels again). Output: fig5_composition_stats_v1.8.csv.
comp <- read.csv("analysis/step3_scrna/composition_all.csv", stringsAsFactors = FALSE)
d <- comp[comp$dataset == "GSE235711", ]
stopifnot(length(unique(d$donor)) == 20, length(unique(d$celltype)) == 19)
d$group <- factor(d$group, levels = c("control", "CRSsNP", "CRSwNP"))
stopifnot(all(table(unique(d[, c("donor", "group")])$group) ==
              c(control = 5, CRSsNP = 5, CRSwNP = 10)))
hl <- function(x, y) {
  w <- suppressWarnings(wilcox.test(y, x, conf.int = TRUE, exact = FALSE))
  c(est = unname(w$estimate), lo = w$conf.int[1], hi = w$conf.int[2])
}
rows <- lapply(sort(unique(d$celltype)), function(ct) {
  x <- d[d$celltype == ct, ]
  g <- split(x$fraction, x$group)
  a <- hl(g$control, g$CRSwNP); b <- hl(g$control, g$CRSsNP)
  data.frame(celltype = ct,
             median_control = median(g$control), median_CRSsNP = median(g$CRSsNP),
             median_CRSwNP = median(g$CRSwNP),
             kw_p = kruskal.test(fraction ~ group, data = x)$p.value,
             hl_CRSwNP_vs_control = a["est"], hl_lo = a["lo"], hl_hi = a["hi"],
             hl_CRSsNP_vs_control = b["est"], hl2_lo = b["lo"], hl2_hi = b["hi"])
})
st <- do.call(rbind, rows)
st$bh_fdr <- p.adjust(st$kw_p, method = "BH")
st <- st[order(st$kw_p), ]
rownames(st) <- NULL
write.csv(st, "analysis/step3_scrna/fig5_composition_stats_v1.8.csv", row.names = FALSE)
print(st, digits = 4)
