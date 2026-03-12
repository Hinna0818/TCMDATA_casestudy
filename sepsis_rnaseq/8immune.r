## 8immune.R - Bindea immune ssGSEA analysis using IOBR built-in signatures
rm(list = ls())
dir.create("output", showWarnings = FALSE)

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "pheatmap", "tibble")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("IOBR", quietly = TRUE)) {
  stop(
    "IOBR is not available in the current R library. ",
    "Please run this script with your R 4.4 environment where IOBR is already installed."
  )
}

library(IOBR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)
library(patchwork)

theme_sci <- function(base_size = 9) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(color = "black", size = base_size),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size)
    )
}

pal_group <- c("Hlty" = "#4DBBD5", "Seps_P" = "#E64B35")
pal_heat <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

pretty_sig_names <- function(x) {
  x <- gsub("_ssGSEA$", "", x)
  x <- gsub("_Bindea_et_al$", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("TFH", "Tfh", x)
  x <- gsub("Tcm", "Tcm", x)
  x <- gsub("Tem", "Tem", x)
  x <- gsub("Tgd", "Tgd", x)
  x <- gsub("iDC", "iDC", x)
  x <- gsub("aDC", "aDC", x)
  x <- gsub("DC", "DC", x)
  x <- gsub("CD8 T", "CD8+ T", x)
  x <- gsub("Th1", "Th1", x)
  x <- gsub("Th2", "Th2", x)
  x <- gsub("Th17", "Th17", x)
  x
}

make_scatter_panel <- function(df, x_var, y_var, x_lab, y_lab) {
  sub_df <- df[, c(x_var, y_var, "condition"), drop = FALSE]
  colnames(sub_df) <- c("x", "y", "condition")

  ggplot(sub_df, aes(x = .data$x, y = .data$y)) +
    geom_point(aes(color = .data$condition), size = 2.0, alpha = 0.8) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "#4D4D4D",
      linewidth = 0.6,
      linetype = "dashed"
    ) +
    scale_color_manual(values = pal_group, name = NULL) +
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title.x = element_text(face = "italic", size = 9),
      axis.title.y = element_text(face = "plain", size = 9)
    )
}

load_iobr_bindea <- function() {
  data("signature_collection", package = "IOBR", envir = environment())
  data("sig_group", package = "IOBR", envir = environment())
  sc <- get("signature_collection", envir = environment())
  sg <- get("sig_group", envir = environment())

  if (is.null(sc) || is.null(sg)) {
    stop("Failed to load IOBR signature_collection or sig_group.")
  }
  if (!"Bindea_et_al" %in% names(sg)) {
    stop("Bindea_et_al group was not found in IOBR sig_group.")
  }

  bindea_names <- sg[["Bindea_et_al"]]
  bindea_names <- setdiff(
    bindea_names,
    c("SW480_cancer_cells_Bindea_et_al", "Normal_mucosa_Bindea_et_al", "Lymph_vessels_Bindea_et_al")
  )
  bindea_sets <- sc[intersect(bindea_names, names(sc))]
  bindea_sets <- lapply(bindea_sets, function(x) unique(as.character(stats::na.omit(x))))
  bindea_sets
}

vsd_mat <- readRDS("output/vsd.rds")
meta <- readRDS("output/sample_meta.rds")
key_targets <- readRDS("output/key_targets.rds")

manual_remove <- c("GSM4683591", "GSM4683603")
keep <- setdiff(colnames(vsd_mat), manual_remove)
vsd_mat <- vsd_mat[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
meta$condition <- factor(meta$condition, levels = c("Hlty", "Seps_P"))

cat("Key targets:", paste(key_targets, collapse = ", "), "\n")

bindea_sets <- load_iobr_bindea()
bindea_sets <- lapply(bindea_sets, function(gs) intersect(gs, rownames(vsd_mat)))
bindea_sets <- bindea_sets[vapply(bindea_sets, length, integer(1)) >= 5]

if (length(bindea_sets) < 4) {
  stop("Too few Bindea signatures matched the expression matrix.")
}

cat("Using Bindea immune signatures:\n")
cat(paste(names(bindea_sets), collapse = "\n"), "\n")

ssgsea_res <- IOBR::calculate_sig_score(
  pdata = NULL,
  eset = vsd_mat,
  signature = bindea_sets,
  method = "ssgsea",
  mini_gene_count = 5,
  adjust_eset = FALSE,
  parallel.size = 1L
)

ssgsea_scores <- ssgsea_res |>
  as.data.frame() |>
  tibble::column_to_rownames("ID")

colnames(ssgsea_scores) <- pretty_sig_names(colnames(ssgsea_scores))
ssgsea_scores <- t(ssgsea_scores)
rownames(ssgsea_scores) <- make.unique(rownames(ssgsea_scores))
ssgsea_scores <- as.data.frame(ssgsea_scores)

sample_order <- rownames(meta)[order(meta$condition)]
score_mat <- as.matrix(ssgsea_scores[, sample_order, drop = FALSE])
ann_col <- data.frame(Group = meta[sample_order, "condition", drop = TRUE])
rownames(ann_col) <- sample_order

pdf("output/fig8a_ssgsea_heatmap.pdf", width = 10.5, height = 7.2)
pheatmap::pheatmap(
  score_mat,
  scale = "row",
  color = pal_heat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = list(Group = pal_group),
  show_colnames = FALSE,
  border_color = NA,
  fontsize_row = 9,
  fontsize = 9,
  main = "Bindea immune ssGSEA"
)
dev.off()

score_long <- t(score_mat) |>
  as.data.frame(check.names = FALSE) |>
  tibble::rownames_to_column("Sample") |>
  dplyr::left_join(
    meta |>
      tibble::rownames_to_column("Sample") |>
      dplyr::select("Sample", "condition"),
    by = "Sample"
  ) |>
  tidyr::pivot_longer(
    cols = all_of(rownames(score_mat)),
    names_to = "Immune_cell",
    values_to = "Abundance"
  )

immune_test <- score_long |>
  dplyr::group_by(.data$Immune_cell) |>
  dplyr::summarise(
    p_value = tryCatch(wilcox.test(Abundance ~ condition)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    FDR = p.adjust(.data$p_value, method = "BH"),
    sig_label = dplyr::case_when(
      .data$FDR < 0.001 ~ "***",
      .data$FDR < 0.01 ~ "**",
      .data$FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

bar_df <- score_long |>
  dplyr::group_by(.data$Immune_cell, .data$condition) |>
  dplyr::summarise(
    mean_abundance = mean(.data$Abundance, na.rm = TRUE),
    se = sd(.data$Abundance, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

cell_order <- bar_df |>
  tidyr::pivot_wider(names_from = "condition", values_from = "mean_abundance") |>
  dplyr::mutate(diff = .data$Seps_P - .data$Hlty) |>
  dplyr::arrange(dplyr::desc(abs(.data$diff))) |>
  dplyr::pull("Immune_cell") |>
  unique()

bar_df$Immune_cell <- factor(bar_df$Immune_cell, levels = cell_order)
p8b_df <- score_long
p8b_df$Immune_cell <- factor(p8b_df$Immune_cell, levels = cell_order)

label_df <- p8b_df |>
  dplyr::group_by(.data$Immune_cell) |>
  dplyr::summarise(y_pos = max(.data$Abundance, na.rm = TRUE) + 0.06 * diff(range(.data$Abundance, na.rm = TRUE)), .groups = "drop") |>
  dplyr::left_join(immune_test, by = "Immune_cell")

p8b <- ggplot(p8b_df, aes(x = .data$Immune_cell, y = .data$Abundance, fill = .data$condition)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.35,
    color = "grey20"
  ) +
  geom_jitter(
    aes(color = .data$condition),
    position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.75),
    size = 1.1,
    alpha = 0.75,
    stroke = 0
  ) +
  geom_text(
    data = label_df,
    aes(x = .data$Immune_cell, y = .data$y_pos, label = .data$sig_label),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  scale_fill_manual(values = pal_group, name = "Group") +
  scale_color_manual(values = pal_group, guide = "none") +
  labs(x = "Immune cells", y = "Immune cell abundance") +
  theme_sci(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

ggsave(
  "output/fig8b_immune_boxplot.pdf",
  p8b,
  width = max(180, 7 * length(unique(bar_df$Immune_cell))),
  height = 140,
  units = "mm",
  dpi = 600,
  device = cairo_pdf
)

expr_key <- t(vsd_mat[key_targets, , drop = FALSE]) |>
  as.data.frame(check.names = FALSE) |>
  tibble::rownames_to_column("Sample")

score_df <- t(score_mat) |>
  as.data.frame(check.names = FALSE) |>
  tibble::rownames_to_column("Sample")

corr_pairs <- expand.grid(
  Gene = key_targets,
  Immune_cell = rownames(score_mat),
  stringsAsFactors = FALSE
)

corr_df <- corr_pairs |>
  dplyr::rowwise() |>
  dplyr::mutate(
    cor = suppressWarnings(cor(expr_key[[Gene]], score_df[[Immune_cell]], method = "spearman")),
    p = suppressWarnings(cor.test(expr_key[[Gene]], score_df[[Immune_cell]], method = "spearman")$p.value)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR = p.adjust(.data$p, method = "BH"),
    R2 = .data$cor^2
  ) |>
  dplyr::mutate(
    label = ifelse(
      .data$FDR < 0.001, paste0(sprintf("%.2f", .data$cor), "\n***"),
      ifelse(
        .data$FDR < 0.01, paste0(sprintf("%.2f", .data$cor), "\n**"),
        ifelse(.data$FDR < 0.05, paste0(sprintf("%.2f", .data$cor), "\n*"), sprintf("%.2f", .data$cor))
      )
    ),
    Gene = factor(.data$Gene, levels = key_targets),
    Immune_cell = factor(.data$Immune_cell, levels = rev(rownames(score_mat)))
  )

p8c <- ggplot(corr_df, aes(x = .data$Gene, y = .data$Immune_cell, fill = .data$cor)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = .data$label), size = 3.0, lineheight = 0.9) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
  labs(x = NULL, y = NULL) +
  theme_sci(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "bold.italic", angle = 0, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  "output/fig8c_target_immune_correlation.pdf",
  p8c,
  width = 155,
  height = max(130, 5.5 * length(unique(corr_df$Immune_cell))),
  units = "mm",
  dpi = 600,
  device = cairo_pdf
)

strong_pairs <- corr_df |>
  dplyr::filter(.data$FDR < 0.05) |>
  dplyr::mutate(direction = ifelse(.data$cor >= 0, "Positive", "Negative")) |>
  dplyr::group_by(.data$Gene, .data$direction) |>
  dplyr::arrange(.data$FDR, dplyr::desc(abs(.data$cor)), .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::arrange(match(.data$Gene, key_targets), dplyr::desc(.data$direction))

if (nrow(strong_pairs) > 0) {
  plot_input_df <- dplyr::left_join(expr_key, score_df, by = "Sample")
  plot_input_df <- dplyr::left_join(
    plot_input_df,
    meta |>
      tibble::rownames_to_column("Sample") |>
      dplyr::select("Sample", "condition"),
    by = "Sample"
  )

  scatter_plots <- lapply(seq_len(nrow(strong_pairs)), function(i) {
    pair_i <- strong_pairs[i, ]
    x_var <- as.character(pair_i$Gene)
    y_var <- as.character(pair_i$Immune_cell)
    stat_label <- paste0(
      "FDR < 0.01",
      "\nR2 = ", sprintf("%.2f", pair_i$R2)
    )

    make_scatter_panel(
      df = plot_input_df,
      x_var = x_var,
      y_var = y_var,
      x_lab = paste0(x_var, " expression"),
      y_lab = paste0(y_var, " abundance")
    ) +
      annotate(
        "label",
        x = -Inf,
        y = Inf,
        label = stat_label,
        hjust = -0.01,
        vjust = 1.01,
        size = 2.3,
        label.size = 0.18,
        fill = "white",
        color = "black",
        label.padding = unit(0.08, "lines")
      )
  })

  p8d <- patchwork::wrap_plots(scatter_plots, ncol = 4, nrow = 2) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "top")

  ggsave(
    "output/fig8d_keytarget_immune_scatter.pdf",
    p8d,
    width = 240,
    height = 120,
    units = "mm",
    dpi = 600,
    device = cairo_pdf
  )
}

saveRDS(ssgsea_scores, "output/ssgsea_scores.rds")
write.csv(bar_df, "output/immune_group_compare.csv", row.names = FALSE)
write.csv(corr_df, "output/keytarget_immune_correlation.csv", row.names = FALSE)
write.csv(strong_pairs, "output/keytarget_immune_scatter_pairs.csv", row.names = FALSE)

