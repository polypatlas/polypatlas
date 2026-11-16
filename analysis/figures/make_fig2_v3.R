# Fig 2 (v3 redraw): The state of the source data.
# Conclusion: of 66 curated series only 32 are core; exclusions are dominated by
# in-vitro material; sample sizes span 2-266; deposited metadata are incomplete.
# Chart types (v3 selection table): (a) 100% stacked composition bar (verdict
# shares), (b) lollipop (ranked exclusion reasons), (c) histogram + density
# (sample-size distribution), (d) lollipop (metadata completeness).
# All numbers computed from curation/curation_table.csv (v1.4) and asserted;
# exclusion-reason categories are the verified curation_report.md S5 breakdown.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(patchwork)})

cur <- read.csv("curation/curation_table.csv", stringsAsFactors = FALSE)
stopifnot(nrow(cur) == 66, sum(cur$n_samples) == 2046)
vd <- as.data.frame(table(factor(cur$verdict, levels = c("core", "pending", "excluded"))))
names(vd) <- c("verdict", "series")
vd$samples <- sapply(vd$verdict, function(v) sum(cur$n_samples[cur$verdict == v]))
stopifnot(vd$series == c(32, 9, 25), vd$samples == c(802, 227, 1017))

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))

## (a) verdict composition: single 100% stacked bar, segment = series count,
## direct labels carry the sample counts (each number printed once).
vd$verdict <- factor(vd$verdict, levels = c("core", "pending", "excluded"))
vd$lab <- c("Core\n802 samples", "Pending\n227 samples", "Excluded\n1,017 samples")
pA <- ggplot(vd, aes(x = "66 series", y = series, fill = verdict)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), position = position_stack(vjust = 0.5),
            size = 2.6, colour = "white", fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = c(core = COL$disease2, pending = COL$accent,
                               excluded = COL$neutral), guide = "none") +
  scale_y_continuous(breaks = c(0, 20, 40, 66), limits = c(0, 66)) +
  labs(x = NULL, y = "Series") + th

## (b) exclusion reasons (verified breakdown, curation_report S5; sums to 25)
reasons <- data.frame(
  reason = c("In-vitro perturbation", "Airway epithelial culture models",
             "Neighbouring-disease cohorts", "Non-expression omics",
             "Superseries", "In-vitro + mouse", "Non-nasal tissue"),
  n = c(7, 6, 4, 3, 3, 1, 1),
  eg = c("7 in-vitro series", "GSE22141 etc.", "incl. GSE250580, n = 266",
         "3 omics types", "GSE235715 etc.", "GSE140900, n = 113", "GSE3874")
)
stopifnot(sum(reasons$n) == 25)
reasons$reason <- factor(reasons$reason, levels = rev(reasons$reason))
pB <- ggplot(reasons, aes(x = n, y = reason)) +
  geom_segment(aes(x = 0, xend = n, y = reason, yend = reason),
               colour = "#CCCCCC", linewidth = 0.5) +
  geom_point(size = 2.6, colour = COL$disease1) +
  geom_text(aes(x = n + 0.25, label = paste0(n, "  (", eg, ")")),
            hjust = 0, size = 2.5, colour = COL$grey_txt) +
  scale_x_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 11.5)) +
  labs(x = "Series excluded", y = NULL) + th +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.3))

## (c) sample-size distribution (log scale), histogram + per-series rug
pC <- ggplot(cur, aes(x = n_samples)) +
  geom_histogram(bins = 18, fill = COL$disease2, colour = "white", linewidth = 0.3) +
  geom_rug(colour = COL$neutral, linewidth = 0.4, length = unit(0.02, "npc")) +
  scale_x_log10(breaks = c(2, 10, 50, 266)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9)) +
  labs(x = "Samples per series (log scale)", y = "Series") + th

## (d) metadata completeness across the 32 core series (computed in-script)
core <- cur[cur$verdict == "core", ]
n_grp   <- sum(core$groups_extractable == "yes")
n_part  <- sum(grepl("^partial", core$groups_extractable))
n_file  <- sum(grepl("^yes \\(resolved", core$groups_extractable))
n_pmid  <- sum(!is.na(core$pmid) & nzchar(trimws(core$pmid)))
n_endot <- sum(grepl("yes", core$endotype_labels, ignore.case = TRUE))
n_treat <- sum(grepl("yes", core$treatment_longitudinal, ignore.case = TRUE))
stopifnot(nrow(core) == 32, n_grp == 29, n_part == 1, n_file == 2,
          n_pmid == 27, n_endot == 6, n_treat == 1)
md <- data.frame(
  field = c("Group fields extractable", "PMID link in record",
            "Endotype labels in metadata", "Longitudinal treatment design"),
  n = c(n_grp, n_pmid, n_endot, n_treat),
  lab = c("29 (+1 partial,\n+2 resolved)", "27", "6", "1")
)
md$field <- factor(md$field, levels = rev(md$field))
pD <- ggplot(md, aes(x = n, y = field)) +
  geom_segment(aes(x = 0, xend = 32, y = field, yend = field),
               colour = "#EEEEEE", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = n, y = field, yend = field),
               colour = "#CCCCCC", linewidth = 0.5) +
  geom_point(size = 2.6, colour = COL$disease2) +
  geom_text(aes(x = n + 1, label = lab), hjust = 0, size = 2.6,
            colour = COL$grey_txt) +
  scale_x_continuous(breaks = c(0, 10, 20, 32), limits = c(0, 50)) +
  labs(x = "Core series (of 32)", y = NULL) + th +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.3))

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))
save_fig(fig, "manuscript/figures/Fig2", width_mm = 180, height_mm = 150)
cat("Fig2 v3 done. core 32/802, pending 9/227, excluded 25/1017\n")
