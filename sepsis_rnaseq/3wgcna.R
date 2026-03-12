## ============================================================
## 3wgcna.R — WGCNA 加权基因共表达网络分析
rm(list = ls())

## ---- 0. 路径 & 包 ----
dir.create("output", showWarnings = FALSE)

library(WGCNA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ggfittext)
library(ggtree)
library(ggfun)
library(aplot)

allowWGCNAThreads()

insert_right <- function(p, p2, width = .05) {
  p + p2 + patchwork::plot_layout(widths = c(1, width))
}

## ---- 1. 加载数据 ----
vsd_mat  <- readRDS("output/vsd.rds")
meta     <- readRDS("output/sample_meta.rds")
sig_degs <- readRDS("output/sig_degs.rds")

cat("VST 矩阵:", nrow(vsd_mat), "genes ×", ncol(vsd_mat), "samples\n")
cat("显著 DEGs:", nrow(sig_degs), "\n")

## ---- 2. 选取高变异基因（top 5000）& 构建表达矩阵 ----
var_genes <- apply(vsd_mat, 1, var)
top5k     <- names(sort(var_genes, decreasing = TRUE))[1:5000]
expr      <- t(vsd_mat[top5k, ])   # samples × genes

cat("WGCNA 输入:", nrow(expr), "samples ×", ncol(expr), "genes\n")

## 样本/基因质控
gsg <- goodSamplesGenes(expr, verbose = 0)
if (!gsg$allOK) {
  expr <- expr[gsg$goodSamples, gsg$goodGenes]
  cat("移除不合格样本/基因后:", nrow(expr), "×", ncol(expr), "\n")
}


## 删除样本后重新进行 soft-threshold 选择
manual_remove <- c("GSM4683591", "GSM4683603")
manual_remove <- intersect(manual_remove, rownames(expr))
if (length(manual_remove) > 0) {
  cat("手动移除样本:", paste(manual_remove, collapse = ", "), "\n")
  expr <- expr[!rownames(expr) %in% manual_remove, , drop = FALSE]
  cat("手动移除后:", nrow(expr), "samples ×", ncol(expr), "genes\n")
}

## ---- 3. 构建性状数据 ----
trait <- meta[rownames(expr), , drop = FALSE]
trait$condition <- factor(trait$condition, levels = c("Hlty", "Seps_P"))
trait$sex <- factor(trait$sex, levels = c("F", "M"))
nGenes   <- ncol(expr)
nSamples <- nrow(expr)

## ---- 4. Soft-threshold 选择 ----
powers <- c(seq(10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(expr, powerVector = powers, verbose = 5)

cat("\n--- Soft Threshold 拟合 ---\n")
print(sft$fitIndices[, c(1, 2, 5)])

## 可视化
sft_df <- sft$fitIndices
sft_df$signed_R2 <- -sign(sft_df$slope) * sft_df$SFT.R.sq

p1 <- ggplot(sft_df, aes(Power, signed_R2)) +
  geom_text(aes(label = Power), color = "red", size = 3) +
  geom_hline(yintercept = 0.8, color = "blue") +
  labs(x = "Soft Threshold (power)",
       y = expression(Scale~Free~Topology~Model~Fit~signed~R^2),
       title = "Scale independence") +
  theme_bw()

p2 <- ggplot(sft_df, aes(Power, mean.k.)) +
  geom_text(aes(label = Power), color = "red", size = 3) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean connectivity") +
  theme_bw()

p1+p2
ggsave("output/sft_threshold.pdf", p1 + p2,
       width = 10, height = 5, device = cairo_pdf)

soft_power <- 12

## ---- 5. 构建共表达网络 ----
set.seed(42)
net <- blockwiseModules(
  expr,
  power             = soft_power,
  maxBlockSize       = 6000,
  TOMType            = "unsigned",
  minModuleSize      = 30,
  reassignThreshold  = 0,
  mergeCutHeight     = 0.25,
  numericLabels      = TRUE,
  pamRespectsDendro  = FALSE,
  saveTOMs           = TRUE,
  saveTOMFileBase    = "output/TOM",
  verbose            = 3
)

## 将数字标签转为颜色标签
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(expr)
cat("\n--- 模块统计 ---\n")
print(table(moduleColors))
cat("模块数量:", length(unique(moduleColors)) - 1, "（不含 grey）\n")

## ---- 6. 模块-表型相关性 ----
## 使用二值性状矩阵做相关性：Condition = 0/1, Sex = 0/1
trait_num <- data.frame(
  Condition = as.numeric(trait$condition == "Seps_P"),
  Sex = as.numeric(trait$sex == "M"),
  row.names = rownames(trait),
  check.names = FALSE
)

## 计算模块特征向量
MEs0 <- moduleEigengenes(expr, moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)
MEs  <- as.data.frame(lapply(MEs, as.numeric), row.names = rownames(MEs))

## 对齐样本并检查异常列
common_samples <- intersect(rownames(MEs), rownames(trait_num))
MEs <- MEs[common_samples, , drop = FALSE]
trait_num <- trait_num[common_samples, , drop = FALSE]

zero_var_me <- vapply(MEs, function(z) sd(z, na.rm = TRUE) == 0, logical(1))
if (any(zero_var_me)) {
  cat("移除零方差模块特征向量:", paste(names(MEs)[zero_var_me], collapse = ", "), "\n")
  MEs <- MEs[, !zero_var_me, drop = FALSE]
}

## 模块-性状相关性 & p 值
moduleTraitCor <- cor(
  as.matrix(MEs),
  as.matrix(trait_num),
  use = "pairwise.complete.obs",
  method = "pearson"
)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

cat("\n--- 模块-表型相关性（前10行）---\n")
print(head(round(moduleTraitCor, 3), 10))

## ---- 7. Fig 5B: 模块-性状热图（demo 风格）----
plot_module_trait <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)

  valid_rows <- rowSums(is.na(x)) == 0
  valid_cols <- colSums(is.na(x)) == 0
  x <- x[valid_rows, valid_cols, drop = FALSE]
  y <- y[valid_rows, valid_cols, drop = FALSE]

  if (nrow(x) == 0 || ncol(x) == 0) {
    stop("moduleTraitCor 全为 NA；请检查 trait 编码或模块特征向量是否含 NA。")
  }

  df <- 1 - x
  module_levels <- rownames(x)
  trait_levels <- colnames(x)

  if (nrow(df) >= 2) {
    module.dendro <- df |>
      dist() |>
      hclust(method = "average")
    module_levels <- rev(module.dendro$labels)
  }
  if (ncol(df) >= 2) {
    trait.dendro <- df |>
      t() |>
      dist() |>
      hclust(method = "average")
    trait_levels <- trait.dendro$labels
  }

  x_long <- x |>
    as.data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = "module") |>
    tidyr::pivot_longer(cols = -all_of("module"), values_to = "cor")

  y_long <- y |>
    as.data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = "module") |>
    tidyr::pivot_longer(cols = -all_of("module"), values_to = "pval")

  dat <- x_long |>
    dplyr::left_join(y_long, by = c("module", "name")) |>
    dplyr::mutate(
      label = paste0(signif(.data$cor, 2), "\n(", signif(.data$pval, 1), ")"),
      module = factor(.data$module, levels = module_levels),
      name = factor(.data$name, levels = trait_levels),
      name = dplyr::recode(as.character(.data$name),
        Condition = "Condition",
        Sex = "Sex"
      )
    )

  p <- dat |>
    ggplot(aes(x = .data$name, y = .data$module, fill = .data$cor)) +
    geom_tile() +
    ggfittext::geom_fit_text(mapping = aes(label = .data$label)) +
    scale_fill_gradient2(
      low = "deepskyblue",
      mid = "white",
      high = "orangered",
      midpoint = 0,
      limits = c(-1, 1),
      name = "correlation"
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  p2 <- dat |>
    dplyr::select(.data$module) |>
    dplyr::distinct() |>
    ggplot(aes(x = "module", y = .data$module, fill = I(gsub("ME", "", .data$module)))) +
    geom_tile() +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none"
    )

  p <- p |> insert_right(p2, width = .05)
  return(p)
}

p5b <- plot_module_trait(moduleTraitCor, moduleTraitPvalue)
p5b
heatmap_height_mm <- max(110, min(180, 5 * nrow(moduleTraitCor) + 24))
ggsave("output/fig5b_module_trait.pdf", p5b,
       width = 100, height = heatmap_height_mm, units = "mm",
       dpi = 600, device = cairo_pdf)

## ---- 7.1 Fig 5B-补充: WGCNA 树状图 + 模块颜色条 ----
plot_wgcna <- function(x, linewidth = .05) {
  colors <- x$colors
  if (is.numeric(colors)) {
    colors <- yulab.utils::get_fun_from_pkg(fun = "labels2colors", pkg = "WGCNA")(colors)
  }

  block_genes <- x$blockGenes[[1]]
  module_vec <- as.character(colors)[block_genes]
  d <- data.frame(
    label = seq_along(block_genes),
    Module = module_vec,
    stringsAsFactors = FALSE
  )

  p <- ggtree::ggtree(x$dendrograms[[1]], layout = "dendrogram",
                      ladderize = FALSE, linewidth = linewidth) +
    ggfun::theme_noxaxis() +
    labs(x = "Height")

  p2 <- ggplot(data = d, aes(x = .data$label, y = "Module", fill = I(.data$Module))) +
    geom_tile() +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )

  p2 |> aplot::insert_top(p, height = 10)
}

p_wgcna <- plot_wgcna(net)
p_wgcna
ggsave("output/fig5b_wgcna_dendrogram.pdf", p_wgcna,
       width = 400, height = 110, units = "mm",
       dpi = 600, device = cairo_pdf)

## ---- 8. 提取关键模块 Hub 基因 ----
## 筛选标准：与任一 trait 列 |r| > 0.5 且 p < 0.05 的模块
sig_mask <- abs(moduleTraitCor) > 0.5 & moduleTraitPvalue < 0.05
sig_modules <- unique(gsub("^ME", "", rownames(moduleTraitCor)[rowSums(sig_mask) > 0]))
sig_modules <- sig_modules[sig_modules != "grey"]

cat("\n显著模块（|r|>0.5, p<0.05）:", paste(sig_modules, collapse = ", "), "\n")

wgcna_hub <- c()
wgcna_hub_tbl <- list()
for (mod in sig_modules) {
  mod_genes <- colnames(expr)[moduleColors == mod]
  mod_row <- paste0("ME", mod)
  me_col <- paste0("ME", mod)

  if (length(mod_genes) == 0) {
    cat("  模块", mod, ": 未找到对应基因，跳过\n")
    next
  }
  if (!me_col %in% colnames(MEs)) {
    cat("  模块", mod, ": 未找到对应 eigengene，跳过\n")
    next
  }

  best_trait <- colnames(moduleTraitCor)[which.max(abs(moduleTraitCor[mod_row, ]))]
  trait_vector <- trait_num[[best_trait]]

  ## GS (Gene Significance): 基因与当前最相关 trait 的相关性
  GS <- abs(cor(expr[, mod_genes, drop = FALSE],
                trait_vector, use = "pairwise.complete.obs"))
  ## MM (Module Membership): 基因与模块特征向量的相关性
  MM <- abs(cor(expr[, mod_genes, drop = FALSE],
                MEs[, me_col, drop = FALSE],
                use = "pairwise.complete.obs"))
  ## Hub 基因：GS > 0.2 且 MM > 0.8
  hub <- mod_genes[GS > 0.2 & MM > 0.8]
  wgcna_hub <- c(wgcna_hub, hub)
  wgcna_hub_tbl[[mod]] <- data.frame(
    gene = mod_genes,
    module = mod,
    trait = best_trait,
    GS = as.numeric(GS),
    MM = as.numeric(MM),
    is_hub = mod_genes %in% hub,
    stringsAsFactors = FALSE
  )
  cat("  模块", mod, "(trait:", best_trait, "):", length(mod_genes), "基因,",
      length(hub), "hub 基因\n")
}
wgcna_hub <- unique(wgcna_hub)
cat("\nWGCNA Hub 基因总数:", length(wgcna_hub), "\n")

## ---- 8.1 提取 turquoise 模块基因 ----
turquoise_genes <- colnames(expr)[moduleColors == "turquoise"]
turquoise_tbl <- data.frame(
  gene = turquoise_genes,
  module = "turquoise",
  stringsAsFactors = FALSE
)

write.csv(turquoise_tbl,
          "output/turquoise_module_genes.csv",
          row.names = FALSE)
saveRDS(turquoise_genes,
        "output/turquoise_module_genes.rds")
cat("turquoise 模块基因数:", length(turquoise_genes), "\n")

## ---- 9. 保存结果 ----
saveRDS(net,            "output/wgcna_net.rds")
saveRDS(moduleColors,   "output/wgcna_moduleColors.rds")
saveRDS(wgcna_hub,      "output/wgcna_hub.rds")
saveRDS(
  list(
    cor = moduleTraitCor,
    pvalue = moduleTraitPvalue,
    trait_matrix = trait_num,
    hub_table = dplyr::bind_rows(wgcna_hub_tbl)
  ),
  "output/module_trait.rds"
)

