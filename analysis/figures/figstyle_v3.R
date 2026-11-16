# figstyle_v3.R — 绘图规范 v3 配套函数库（v2 函数沿用，2026-10-30 升 v3）
# 用法：每张图的脚本开头 source("~/workspace/custom-skills/figure_style_v3/figstyle_v3.R")
# 规范全文见同目录 SKILL.md。色板为 AGENTS.md 锁定值，不得在此另配。

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

## ---- 色板（锁定） -----------------------------------------------------------
COL <- list(
  disease1 = "#E64B35",  # 主病
  disease2 = "#3C5488",  # 对照病
  up       = "#DC0000",  # 上调
  down     = "#4DBBD5",  # 下调
  neutral  = "#8491B4",  # 中性/邻居
  accent   = "#7E6148",  # 点缀棕
  grey_txt = "#333333",
  grey_ax  = "#666666"
)

# 分类色扩展（>2 组时）：以锁定色相为锚的 12 色，禁止 ggplot 默认色
pal_cat <- function(n) {
  base <- c("#E64B35", "#3C5488", "#4DBBD5", "#7E6148", "#8491B4",
            "#00A087", "#F39B7F", "#91D1C2", "#DC0000", "#7E7E7E",
            "#B09C00", "#5B4E77")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

## ---- 统一 theme --------------------------------------------------------------
# base_size 按最终出图尺寸设定：主图 10，补充图 9。字号梯度见 SKILL.md §三。
theme_fig <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      axis.line.x        = element_line(colour = COL$grey_ax, linewidth = 0.4),
      axis.line.y        = element_line(colour = COL$grey_ax, linewidth = 0.4),
      axis.ticks         = element_line(colour = COL$grey_ax, linewidth = 0.4),
      axis.text          = element_text(colour = COL$grey_txt, size = base_size - 2),
      axis.title         = element_text(colour = COL$grey_txt, size = base_size - 1),
      legend.position    = "right",
      legend.title       = element_text(size = base_size - 2),
      legend.text        = element_text(size = base_size - 2),
      legend.key         = element_blank(),
      legend.background  = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(size = base_size - 1, face = "bold",
                                        colour = COL$grey_txt),
      plot.tag           = element_text(size = base_size + 2, face = "bold"),
      text               = element_text(family = "Arial")
    )
}

## ---- 标准保存：按 mm 出图，PDF + PNG 300 dpi --------------------------------
# width_mm: 单栏 89 / 1.5 栏 130 / 双栏 180；height_mm ≤ 240
save_fig <- function(plot, path_base, width_mm = 180, height_mm = 120) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(path_base, ".pdf"), width = w, height = h,
                       family = "Arial")
  print(plot); grDevices::dev.off()
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(path_base, ".png"), width = w, height = h,
                  units = "in", res = 300)
  } else {
    grDevices::png(paste0(path_base, ".png"), width = w, height = h,
                   units = "in", res = 300)
  }
  print(plot); grDevices::dev.off()
  invisible(path_base)
}

## ---- 图型封装 -----------------------------------------------------------------

# 分布比较：小提琴 + 逐样本点 + 中位数点（替代均值条图）
# df 需含列 x（分组）, y（数值）；可选列 fill（默认同 x）
plot_dist <- function(df, x = "group", y = "value", ylab = y, cols = NULL) {
  d <- df; d$.x <- d[[x]]; d$.y <- d[[y]]
  if (is.null(cols)) cols <- pal_cat(length(unique(d$.x)))
  ggplot(d, aes(.x, .y, fill = .x)) +
    geom_violin(trim = FALSE, colour = NA, alpha = 0.75, scale = "width") +
    ggbeeswarm::geom_quasirandom(width = 0.18, size = 1.0, colour = "#222222",
                                 alpha = 0.75) +
    stat_summary(fun = median, geom = "point", size = 2.4,
                 colour = "white", shape = 21, fill = "black", stroke = 0.6) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = ylab) +
    theme_fig() + theme(legend.position = "none")
}

# 雨云图：半提琴（ggdist）+ 箱线 + 逐点 —— 分布形态本身是结论时用
plot_raincloud <- function(df, x = "group", y = "value", ylab = y, cols = NULL) {
  stopifnot(requireNamespace("ggdist", quietly = TRUE))
  d <- df; d$.x <- d[[x]]; d$.y <- d[[y]]
  if (is.null(cols)) cols <- pal_cat(length(unique(d$.x)))
  ggplot(d, aes(.x, .y, fill = .x)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = -0.12, .width = 0,
                         point_colour = NA, alpha = 0.75) +
    geom_boxplot(width = 0.12, outlier.shape = NA, colour = COL$grey_ax,
                 fill = "white", linewidth = 0.4) +
    ggbeeswarm::geom_quasirandom(width = 0.08, size = 1.0, colour = "#222222",
                                 alpha = 0.75, dodge.width = 0.12) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = ylab) +
    theme_fig() + theme(legend.position = "none")
}

# 棒棒糖：排名/Top-N（替代横条图）；df 需含 y（类别，已按 x 排序）, x（数值）
plot_lollipop <- function(df, x = "value", y = "label", xlab = x,
                          col = COL$disease1) {
  d <- df; d$.x <- d[[x]]; d$.y <- d[[y]]
  ggplot(d, aes(.x, .y)) +
    geom_segment(aes(x = 0, xend = .x, y = .y, yend = .y),
                 colour = "#CCCCCC", linewidth = 0.5) +
    geom_point(size = 2.6, colour = col) +
    labs(x = xlab, y = NULL) +
    theme_fig() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.3))
}

# 富集点图：x=富集强度, y=通路, 点大小=基因数, 色=-log10 FDR（替代通路条图）
plot_enrich_dot <- function(df, x = "rich_factor", y = "pathway",
                            size = "count", color = "neglog10FDR",
                            xlab = "Rich factor") {
  d <- df; d$.x <- d[[x]]; d$.y <- d[[y]]; d$.s <- d[[size]]; d$.c <- d[[color]]
  ggplot(d, aes(.x, .y, size = .s, colour = .c)) +
    geom_point(alpha = 0.9) +
    scale_colour_gradient(low = COL$down, high = COL$up,
                          name = expression(-log[10]~FDR)) +
    scale_size_continuous(name = "Genes", range = c(2, 6)) +
    labs(x = xlab, y = NULL) +
    theme_fig()
}

# 配对斜率图：个体灰线 + 组均值粗线（替代分组条图）
# df 需含列 id, time（两水平）, value；可选列 group 决定粗线颜色
plot_paired_slope <- function(df, id = "id", time = "time", y = "value",
                              group = NULL, ylab = y) {
  d <- df; d$.id <- d[[id]]; d$.t <- d[[time]]; d$.y <- d[[y]]
  p <- ggplot(d, aes(.t, .y, group = .id)) +
    geom_line(colour = "#BBBBBB", alpha = 0.6, linewidth = 0.4) +
    geom_point(size = 1.0, colour = "#888888")
  if (!is.null(group)) {
    d$.g <- d[[group]]
    m <- d %>% group_by(.g, .t) %>% summarise(.y = mean(.y), .groups = "drop")
    p <- ggplot() +
      geom_line(data = d, aes(.t, .y, group = .id),
                colour = "#BBBBBB", alpha = 0.6, linewidth = 0.4) +
      geom_line(data = m, aes(.t, .y, group = .g, colour = .g),
                linewidth = 1.4) +
      geom_point(data = m, aes(.t, .y, colour = .g), size = 2.6) +
      scale_colour_manual(values = pal_cat(length(unique(m$.g))), name = NULL)
  }
  p + labs(x = NULL, y = ylab) + theme_fig()
}

# 森林图：效应+CI，按效应排序（df 需含 label, est, lo, hi；可选 p 文本列注）
plot_forest <- function(df, label = "label", est = "est", lo = "lo",
                        hi = "hi", xlab = "Effect", null = 0) {
  d <- df; d$.l <- d[[label]]; d$.e <- d[[est]]; d$.lo <- d[[lo]]; d$.hi <- d[[hi]]
  ggplot(d, aes(.e, .l)) +
    geom_vline(xintercept = null, colour = "#999999", linetype = "dashed",
               linewidth = 0.4) +
    geom_errorbarh(aes(xmin = .lo, xmax = .hi), height = 0.22,
                   colour = COL$grey_ax, linewidth = 0.5) +
    geom_point(size = 2.4, colour = COL$disease1) +
    labs(x = xlab, y = NULL) +
    theme_fig() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.3))
}
