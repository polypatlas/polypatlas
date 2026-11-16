# Fig 6 (v3 redraw): Known-marker panel recovery heatmap.
# Conclusion: the 8 pre-defined type-2 markers are recovered across layers -
# bulk disease 49/56 = 87.5% direction-correct, dupilumab treatment reverses
# 7/8 (IL33 unchanged), single-cell layer sparse for 4 factors (not measured).
# Single-panel tile heatmap; cell text carries both percentages, so no colour
# legend is printed (v3 V10: numbers appear once; V4: direct annotation over
# legend). Not-measured cells are grey, not mid-scale coloured.
# Data prep replicates analysis/step6_validation/step6_figures.R exactly.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(dplyr)})

mr <- read.csv("analysis/step6_validation/marker_recovery.csv", stringsAsFactors = FALSE)
ha <- mr %>% filter(!grepl("spatial", layer), layer %in% c("bulk", "scrna")) %>%
  mutate(cls = ifelse(layer == "bulk", paste0("Bulk: ", contrast_class),
                      "scRNA: GSE235711 CRSwNP"),
         frac_dir_correct = as.numeric(frac_dir_correct),
         frac_dir_and_sig = as.numeric(frac_dir_and_sig))
## headline-number assertions against the same source table
bd <- ha[ha$cls == "Bulk: disease", ]
stopifnot(nrow(bd) == 8, all(bd$n_measured == 7),
          round(sum(bd$frac_dir_correct * bd$n_measured)) == 49,  # 49/56 = 87.5%
          round(sum(bd$frac_dir_and_sig * bd$n_measured)) == 20)   # 20/56 = 35.7%
bt <- ha[ha$cls == "Bulk: treatment", ]
stopifnot(sum(bt$frac_dir_correct == 1) == 7,                      # 7/8 reversed
          bt$frac_dir_correct[bt$gene == "IL33"] == 0)             # IL33 unchanged
sc <- ha[ha$cls == "scRNA: GSE235711 CRSwNP", ]
stopifnot(sum(is.na(sc$frac_dir_correct)) == 4)                    # 4 sparse factors

ha <- ha %>% mutate(
  lbl = ifelse(is.na(frac_dir_correct), "not\nmeasured",
               sprintf("%.0f%%\n(%.0f%%)", 100 * frac_dir_correct, 100 * frac_dir_and_sig)),
  gene = factor(gene, levels = rev(c("CLC", "CCL26", "CCL13", "IL33", "TSLP",
                                     "POSTN", "SERPINB2", "ALOX15"))),
  cls = factor(cls, levels = c("Bulk: disease", "Bulk: treatment",
                               "scRNA: GSE235711 CRSwNP")))
cls_lab <- c("Bulk: disease" = "Bulk · disease",
             "Bulk: treatment" = "Bulk · dupilumab treatment",
             "scRNA: GSE235711 CRSwNP" = "Single-cell · GSE235711")

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))
p <- ggplot(ha, aes(x = cls, y = gene, fill = frac_dir_correct)) +
  geom_tile(colour = "grey40", linewidth = 0.4) +
  geom_text(aes(label = lbl), size = 3.2, colour = "black") +
  scale_fill_gradient2(low = COL$up, mid = "#F5EFE0", high = "#00A087",
                       midpoint = 0.5, limits = c(0, 1), na.value = "grey85",
                       guide = "none") +
  scale_x_discrete(labels = cls_lab) +
  labs(x = "Direction-correct % (parentheses: also FDR < 0.05)", y = NULL) + th +
  theme(axis.text.x = element_text(size = 8.5), panel.grid.major = element_blank())
save_fig(p, "manuscript/figures/Fig6", width_mm = 180, height_mm = 115)
cat("Fig6 v3 done. bulk 49/56 = 87.5%; treatment 7/8; scRNA not-measured genes: 4\n")
