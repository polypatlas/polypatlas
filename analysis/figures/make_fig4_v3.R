# Fig 4 (v3 redraw): Evidence and consistency.
# Conclusions: (a) the top-30 consensus-up genes are directionally coherent
# across the 16 bulk cohorts; (b) CLC is up in nearly every disease cohort and
# reversed by dupilumab, RE meta +1.865 [1.317, 2.413]; (c) consensus genes are
# direction-consistent in a median 0.90 (up) / 0.80 (down) of measured cohorts.
# Chart types: heatmap (anchor; row order = pheatmap clustering of the same
# matrix), forest plot with aligned n column, faceted distribution histograms.
# Data prep identical to the original script; presentation follows v3.
# The cohort-design annotation bar tried in drafting was dropped: it would
# have been this figure's third legend (v3 V4 limit is 2); the matrix itself
# is the single information layer, as in the submitted version.
source("analysis/figures/figstyle_v3.R")
suppressPackageStartupMessages({library(pheatmap); library(jsonlite); library(patchwork)})

ev  <- read.csv("analysis/step2_bulk/evidence_matrix_logFC.csv", row.names = 1, check.names = FALSE)
stopifnot(ncol(ev) == 16)
dj <- fromJSON("analysis/step6_validation/consensus_genes_disease_only.json")
stopifnot(nrow(dj$up) == 825, nrow(dj$down) == 196)
updf <- dj$up[order(-dj$up$median_logFC), ]
up30 <- updf$gene[1:30]
m <- as.matrix(ev[up30, ])
m <- pmin(pmax(m, -4), 4)   # display cap as in the original; NA = not measured (grey)

## row order from the same pheatmap clustering used by the submitted figure
ph <- pheatmap(m, cluster_rows = TRUE, cluster_cols = FALSE, silent = TRUE,
               na_col = "grey85",
               color = colorRampPalette(c("#2e75b6", "white", "#c00000"))(60),
               breaks = seq(-4, 4, length.out = 61))
gene_ord <- rownames(m)[ph$tree_row$order]
hm <- data.frame(gene = rep(rownames(m), times = ncol(m)),
                 cohort = rep(colnames(m), each = nrow(m)),
                 value = as.vector(m))
hm$gene <- factor(hm$gene, levels = rev(gene_ord))
hm$cohort <- factor(hm$cohort, levels = colnames(m))

th <- theme_fig(10) + theme(axis.text = element_text(size = 7.5),
                            legend.text = element_text(size = 7.5),
                            legend.title = element_text(size = 8))
pA <- ggplot(hm, aes(x = cohort, y = gene, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = COL$down, mid = "white", high = COL$up,
                       limits = c(-4, 4), na.value = "grey85", name = "logFC") +
  labs(x = NULL, y = NULL) + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid.major = element_blank(),
        legend.position = "right")

## (b) CLC forest (data prep as original)
runlog <- read.csv("analysis/step2_bulk/runlog.csv", stringsAsFactors = FALSE)
prim <- list(GSE136825 = "polyp_vs_ctrl", GSE150846 = "NP_ILC2_vs_blood_ILC2",
             GSE174681 = "NP_ASC_vs_PBLnaive", GSE179265 = "CRSwNP_vs_ctrl",
             GSE189690 = "NERD_vs_nonNERD", GSE194282 = "NP_vs_ctrl",
             GSE198950 = "CRSwNP_vs_NonCRS", GSE213809 = "NP_vs_PB_stimulated",
             GSE267331 = "dupilumab_post_vs_pre", GSE313186 = "exposed_vs_ctrl",
             GSE330761 = "CRSwNP_MCTC_vs_CRSsNP_MCTC", GSE338426 = "CRSwNP_vs_ctrl",
             GSE36830 = "CRSwNP_vs_UTctrl", GSE69093 = "case_vs_ctrl_celltypeAdj",
             GSE70898_70900 = "polypT_vs_nasalT", GSE72713 = "CRSwNP_vs_ctrl")
de_files <- list.files("analysis/step2_bulk/de", pattern = "[.]csv$", full.names = TRUE)
file_for <- function(gse, contrast) {
  stems <- sub("[.]csv$", "", sub(paste0("^", gse, "_"), "", basename(de_files)))
  hit <- de_files[tolower(stems) == tolower(contrast)]
  if (length(hit) != 1) stop(paste("forest file problem:", gse, contrast))
  hit
}
forest <- do.call(rbind, lapply(names(prim), function(g) {
  d <- read.csv(file_for(g, prim[[g]]), stringsAsFactors = FALSE)
  r <- d[d$symbol == "CLC", ]
  rr <- runlog[runlog$gse == g & runlog$contrast == prim[[g]], ]
  if (nrow(r) == 0) {
    return(data.frame(cohort = g, logFC = NA_real_, se = NA_real_,
                      n = paste0(rr$n_case[1], "v", rr$n_control[1]),
                      treat = FALSE, measured = FALSE))
  }
  stopifnot(nrow(r) == 1)
  data.frame(cohort = g, logFC = r$logFC, se = abs(r$logFC / r$t),
             n = paste0(rr$n_case[1], "v", rr$n_control[1]),
             treat = (g == "GSE267331"), measured = TRUE)
}))
forest$lo <- forest$logFC - 1.96 * forest$se
forest$hi <- forest$logFC + 1.96 * forest$se
stopifnot(sum(forest$measured) == 15, !forest$measured[forest$cohort == "GSE213809"])
f1 <- read.csv("analysis/step9_features/f1_meta.csv", stringsAsFactors = FALSE)
f1c <- f1[f1$gene == "CLC", ]
stopifnot(nrow(f1c) == 1, round(f1c$meta_logFC, 3) == 1.865,
          round(f1c$ci_lo, 3) == 1.317, round(f1c$ci_hi, 3) == 2.413)
forest <- forest[order(forest$logFC, na.last = TRUE), ]
meta_df <- data.frame(cohort = "Meta (RE, k = 10)", logFC = f1c$meta_logFC,
                      lo = f1c$ci_lo, hi = f1c$ci_hi)
forest$cohort <- factor(forest$cohort, levels = c("Meta (RE, k = 10)", forest$cohort))
meta_df$cohort <- factor(meta_df$cohort, levels = levels(forest$cohort))

pB <- ggplot(subset(forest, measured),
             aes(x = logFC, y = cohort, colour = treat, shape = treat)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL$neutral, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, linewidth = 0.5) +
  geom_point(size = 2.2) +
  geom_errorbarh(data = meta_df, aes(xmin = lo, xmax = hi, y = cohort),
                 inherit.aes = FALSE, colour = COL$up, linewidth = 0.9, height = 0.3) +
  geom_point(data = meta_df, aes(x = logFC, y = cohort), inherit.aes = FALSE,
             colour = COL$up, fill = COL$up, shape = 23, size = 3.4) +
  geom_text(data = subset(forest, measured), aes(x = 13.4, label = n),
            hjust = 0, size = 2.5, colour = COL$grey_txt) +
  scale_colour_manual(values = c("FALSE" = COL$disease2, "TRUE" = COL$disease1),
                      labels = c("FALSE" = "Disease cohort",
                                 "TRUE" = "Dupilumab\n(reversal expected)"),
                      name = NULL) +
  scale_shape_manual(values = c("FALSE" = 19, "TRUE" = 17), guide = "none") +
  scale_x_continuous(breaks = c(0, 4, 8, 12)) +
  coord_cartesian(xlim = c(-2.2, 14.8)) +   # clip wide CIs at the frame, do not drop them
  labs(x = "CLC logFC (95% CI)", y = NULL) + th +
  theme(legend.position = "top", axis.text.y = element_text(size = 7.5),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.3))

## (c) direction-consistency distribution of the 1,021 consensus genes
cf <- do.call(rbind, lapply(c("up", "down"), function(dd) {
  data.frame(dir = dd, frac = dj[[dd]]$frac_consistent)
}))
med_u <- median(cf$frac[cf$dir == "up"]); med_d <- median(cf$frac[cf$dir == "down"])
stopifnot(round(med_u, 2) == 0.90, round(med_d, 2) == 0.80, nrow(cf) == 1021)
cf$dir <- factor(cf$dir, levels = c("up", "down"))
strip_lab <- c(up = sprintf("Up · n = 825 · median %.2f", med_u),
               down = sprintf("Down · n = 196 · median %.2f", med_d))
pC <- ggplot(cf, aes(x = frac, fill = dir)) +
  geom_histogram(bins = 12, colour = "grey25", linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~dir, ncol = 1, labeller = labeller(dir = strip_lab), scales = "free_x") +
  geom_vline(data = data.frame(dir = factor(c("up", "down"), levels = c("up", "down")),
                               med = c(med_u, med_d)),
             aes(xintercept = med), linetype = "dashed", colour = COL$grey_txt) +
  scale_fill_manual(values = c(up = COL$up, down = COL$down), guide = "none") +
  labs(x = "Fraction of cohorts with consistent direction",
       y = "Consensus genes") + th

fig <- pA + pB + pC +
  plot_layout(design = "AB\nAC", widths = c(0.48, 0.52),
              heights = c(0.44, 0.56)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))
save_fig(fig, "manuscript/figures/Fig4", width_mm = 180, height_mm = 190)
cat("Fig4 v3 done. medians:", med_u, med_d, "\n")
