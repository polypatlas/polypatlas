# Fig 9 (v3 redraw): Pathway enrichment results.
# Conclusion: 20 up-regulated GO-BP terms recur (FDR < 0.05) in >= 3 of the 4
# major bulk contrasts. Chart type: enrichment dot plot (v3 selection table;
# size = gene count, colour = -log10 FDR on the locked gradient).
# Term selection replicates figE1_enrichment_fig9corrected.R exactly.
source("analysis/figures/figstyle_v3.R")

ora <- read.csv("analysis/step8_expansion/enrichment/enrichment_ora.csv", stringsAsFactors = FALSE)
key_contrasts <- c("GSE136825_polyp_vs_ctrl", "GSE179265_CRSwNP_vs_ctrl",
                   "GSE338426_CRSwNP_vs_ctrl", "GSE338426_eCRSwNP_vs_ctrl")
short <- c("GSE136825\npolyp vs ctrl", "GSE179265\nCRSwNP vs ctrl",
           "GSE338426\nCRSwNP vs ctrl", "GSE338426\neCRSwNP vs ctrl")
names(short) <- key_contrasts
d <- ora[ora$contrast %in% key_contrasts & ora$direction == "up" &
         ora$database == "GO_BP" & ora$FDR < 0.05, ]
d$neglog <- -log10(pmax(d$FDR, 1e-300))
dup <- unique(d[, c("contrast", "term_name")])
rec_tab <- table(factor(dup$term_name))
rec <- names(rec_tab[rec_tab >= 3])
mx <- tapply(d$neglog[d$term_name %in% rec], d$term_name[d$term_name %in% rec], max)
ord <- order(rec_tab[rec], mx[rec], decreasing = TRUE)
rec <- rec[ord[1:min(20, length(rec))]]
stopifnot(length(rec) == 20)
d2 <- d[d$term_name %in% rec, ]
d2$term_lab <- stringr::str_wrap(d2$term_name, width = 42)
d2$term_lab <- factor(d2$term_lab, levels = rev(unique(d2$term_lab[order(d2$neglog, decreasing = TRUE)])))
d2$cshort <- factor(short[d2$contrast], levels = short)

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))
p <- ggplot(d2, aes(x = cshort, y = term_lab, size = count, colour = neglog)) +
  geom_point(alpha = 0.9) +
  scale_colour_gradient(low = COL$down, high = COL$up, name = "-log10 FDR") +
  scale_size_continuous(name = "Genes", range = c(2, 6.5)) +
  labs(x = NULL, y = NULL) + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 7.5))
save_fig(p, "manuscript/figures/Fig9", width_mm = 180, height_mm = 140)
cat("Fig9 v3 done. recurrent terms:", length(rec), "| plotted cells:", nrow(d2), "\n")
print(as.data.frame(table(d2$cshort)))
