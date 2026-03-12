## ============================================================
## 2degs.R — DESeq2 差异表达分析
rm(list = ls())

## ---- 0. 路径 & 包 ----
dir.create("output", showWarnings = FALSE)

library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ivolcano)

## ---- 1. 加载数据 ----
dds    <- readRDS("output/dds.rds")
meta   <- readRDS("output/sample_meta.rds")
counts <- readRDS("output/counts_filtered.rds")

cat("加载数据:", nrow(dds), "genes ×", ncol(dds), "samples\n")
cat("分组:", table(meta$condition), "\n")

## ---- 2. 运行 DESeq2 ----
dds <- DESeq(dds)

## 查看 resultsNames 确认 coef 名称
cat("Results names:", resultsNames(dds), "\n")

## LFC shrinkage（apeglm 方法，减少低表达基因噪声）
res <- lfcShrink(dds, coef = "condition_Seps_P_vs_Hlty", type = "apeglm")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

cat("\n--- DESeq2 结果摘要 ---\n")
summary(res, alpha = 0.05)

## ---- 3. 筛选显著 DEGs ----
sig_degs <- res_df %>%
  filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) >= 1) %>%
  arrange(padj)

cat("\n显著 DEGs（padj<0.05, |log2FC|>1）:", nrow(sig_degs), "\n")
cat("  上调:", sum(sig_degs$log2FoldChange > 0), "\n")
cat("  下调:", sum(sig_degs$log2FoldChange < 0), "\n")

## ---- 4. VST 变换（供下游 WGCNA / ML 使用）----
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

## ---- 5. Fig 5A: 火山图 ----
res_plot <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj <= 0.05 & log2FoldChange >=  1 ~ "Up",
      padj <= 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Up", "NS", "Down")),
    neg_log10p = -log10(pmax(padj, 1e-300))
  )

## 统计单阈值模式下各分组数量，供图例显示
legend_counts <- res_plot %>%
  mutate(
    sig3 = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Not_Significant"
    )
  ) %>%
  count(sig3)

legend_labels <- c(
  Up = "Up",
  Not_Significant = "NS",
  Down = "Down"
)

legend_levels <- names(legend_labels)
for (lvl in legend_levels) {
  n_val <- legend_counts$n[match(lvl, legend_counts$sig3)]
  if (is.na(n_val)) n_val <- 0L
  legend_labels[lvl] <- paste0(legend_labels[lvl], " (", n_val, ")")
}

p5a <- ivolcano(res_plot, 
               logFC_col = "log2FoldChange",
               pval_col = "padj",
               pval_cutoff = 0.05,
               logFC_cutoff = 1,
               gene_col = "gene",
               title = NULL,
               size_by = "manual",
               point_size = list(base = 1, medium = 1.2, large = 1.5),
               top_n = 10,
               interactive = FALSE) +
  scale_color_manual(
    values = c(
      Up = "#D55E5E",
      Not_Significant = "#CFCFCF",
      Down = "#4C88C7"
    ),
    breaks = c("Up", "Not_Significant", "Down"),
    labels = legend_labels,
    name = NULL
  ) +
  labs(
    x = expression(Log[2]~"(Fold Change)"),
    y = expression(-Log[10]~"("*italic("Adjusted P-value")*")")
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold", size = 11, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = grid::unit(3.5, "mm"),
    legend.text = element_text(size = 8, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(6, 8, 6, 6)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
print(p5a)
  
ggsave("output/fig5a_volcano.pdf", p5a,units = "mm", dpi = 600, device = cairo_pdf)

## ---- 6. 保存结果 ----
saveRDS(dds,     "output/dds_deseq.rds")
saveRDS(res_df,  "output/res_df.rds")
saveRDS(sig_degs,"output/sig_degs.rds")
saveRDS(vsd_mat, "output/vsd.rds")


