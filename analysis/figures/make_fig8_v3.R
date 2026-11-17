# Fig 8 (v3 redraw; v1.7 data, 2026-11-17): Internal reproduction audit
# across the 12 cohorts with testable claims. Conclusion: of 60 testable
# claims, 23 were reproduced (direction + FDR) and 53/60 = 88.3% are
# direction-consistent; 3 were not reproduced (all in the proxy stratum).
# Chart type: stacked count bars per cohort (composition of outcomes; the
# v3-sanctioned bar use). Data prep replicates
# analysis/step6_validation/step6_figures.R; totals asserted.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(dplyr)})

au <- read.csv("analysis/step6_validation/reproduction_audit_v1.7.csv", stringsAsFactors = FALSE)
testable <- au %>% filter(verdict %in% c("reproduced", "direction_only",
                                         "not_reproduced", "not_measured")) %>%
  mutate(cohort = paste0(gse, "\nPMID ", pmid))
cnt <- testable %>% count(cohort, verdict)
tot <- tapply(cnt$n, cnt$verdict, sum)
stopifnot(length(unique(testable$cohort)) == 12, sum(cnt$n) == 60,
          tot[["reproduced"]] == 23, tot[["direction_only"]] == 30,
          tot[["not_measured"]] == 4, tot[["not_reproduced"]] == 3)
## ggplot2 stacks the first factor level at the top; this level order
## reproduces the submitted stacking (bottom -> top: not reproduced,
## not measured, direction only, reproduced)
cnt$verdict <- factor(cnt$verdict,
                      levels = c("reproduced", "direction_only",
                                 "not_measured", "not_reproduced"))
cnt$cohort <- factor(cnt$cohort, levels = sort(unique(cnt$cohort)))

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))
p <- ggplot(cnt, aes(x = cohort, y = n, fill = verdict)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 2.6, colour = "white", fontface = "bold") +
  scale_fill_manual(name = "Outcome",
                    values = c(reproduced = "#00A087", direction_only = "#F39B7F",
                               not_measured = COL$neutral, not_reproduced = COL$up),
                    breaks = c("reproduced", "direction_only", "not_measured",
                               "not_reproduced"),
                    labels = c("Reproduced (direction + FDR < 0.05)",
                               "Direction only", "Not measured",
                               "Not reproduced")) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20)) +
  labs(x = NULL, y = "Claims") + th +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5))
save_fig(p, "manuscript/figures/Fig8", width_mm = 180, height_mm = 120)
cat("Fig8 v3 (v1.7) done. 23 reproduced / 30 direction-only / 4 not measured / 3 not reproduced = 60\n")
