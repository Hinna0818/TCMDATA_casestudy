## 1dataloader.R — GSE154918 数据下载、导入与预处理
rm(list = ls())

## ---- 0. 路径设置 ----
dir.create("output", showWarnings = FALSE, recursive = TRUE)

## ---- 1. 加载必要的包 ----
library(GEOquery)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

## ---- 2. 数据文件下载 ----
## 2.1 Raw counts 文件
raw_counts_gz <- "GSE154918_raw_counts_GRCh38.p13_NCBI.tsv.gz"
if (!file.exists(raw_counts_gz)) {
  download.file(
    url      = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154918/suppl/GSE154918_raw_counts_GRCh38.p13_NCBI.tsv.gz",
    destfile = raw_counts_gz,
    mode     = "wb"
  )
}

## 2.2 基因注释文件
annot_gz <- "Human.GRCh38.p13.annot.tsv.gz"
if (!file.exists(annot_gz)) {
  download.file(
    url      = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154918/suppl/Human.GRCh38.p13.annot.tsv.gz",
    destfile = annot_gz,
    mode     = "wb"
  )
}

## ---- 3. 读取 Raw Counts ----
## 文件格式：第 1 列 = GeneID (NCBI Entrez ID)，后续列 = GSM 样本号
## 值为整数 raw counts
counts_raw <- read.delim(
  raw_counts_gz,
  header         = TRUE,
  row.names      = 1,     # GeneID 作为行名
  check.names    = FALSE,
  stringsAsFactors = FALSE
)
cat("Raw counts 维度:", nrow(counts_raw), "genes ×", ncol(counts_raw), "samples\n")

## ---- 4. 读取基因注释 ----
## 列：GeneID, Symbol, Description, Synonyms, GeneType, EnsemblGeneID,
##     Status, ChrAcc, ChrStart, ChrStop, Orientation, Length,
##     GOFunctionID, GOProcessID, GOComponentID, GOFunction, GOProcess, GOComponent
annot <- read.delim(
  annot_gz,
  header         = TRUE,
  stringsAsFactors = FALSE
)
cat("注释文件维度:", nrow(annot), "genes\n")

## 提取映射表：GeneID → Symbol, GeneType
gene_map <- annot %>%
  select(GeneID, Symbol, GeneType, EnsemblGeneID, Description) %>%
  distinct(GeneID, .keep_all = TRUE)

cat("GeneID → Symbol 映射:", nrow(gene_map), "entries\n")

## ---- 5. 通过 GEOquery 获取样本元数据 ----
## GEOquery 下载 Series Matrix 解析分组信息
## 样本特征字段：status (Hlty/Seps_P/...) 和 Sex (M/F)
gse <- getGEO("GSE154918", GSEMatrix = TRUE, getGPL = FALSE)
pdata <- pData(gse[[1]])

## 提取关键元数据列
meta <- pdata %>%
  transmute(
    geo_accession = geo_accession,
    title         = title,
    ## 解析 characteristics 字段
    ## characteristics_ch1 格式：  "status: Hlty"
    ## characteristics_ch1.1 格式："Sex: M"
    status = gsub("^status:\\s*", "", characteristics_ch1),
    sex    = gsub("^Sex:\\s*",    "", characteristics_ch1.1)
  ) %>%
  as.data.frame()

rownames(meta) <- meta$geo_accession

cat("\n--- 全部样本分组统计 ---\n")
print(table(meta$status))
cat("--- 性别分布 ---\n")
print(table(meta$status, meta$sex))

## ---- 6. 二分组筛选：仅保留 Hlty 和 Seps_P ----
meta_sub <- meta %>%
  filter(status %in% c("Hlty", "Seps_P"))

meta_sub$condition <- factor(meta_sub$status, levels = c("Hlty", "Seps_P"))

cat("\n筛选后样本数:", nrow(meta_sub), "\n")
cat("  Hlty  :", sum(meta_sub$condition == "Hlty"), "\n")
cat("  Seps_P:", sum(meta_sub$condition == "Seps_P"), "\n")

## 对应筛选 counts 矩阵（列名 = GSM 编号，与 meta_sub 行名匹配）
shared_samples <- intersect(colnames(counts_raw), rownames(meta_sub))
counts_sub <- counts_raw[, shared_samples]
meta_sub   <- meta_sub[shared_samples, ]

cat("匹配后 counts 维度:", nrow(counts_sub), "×", ncol(counts_sub), "\n")

## ---- 7. GeneID → Gene Symbol 映射 ----
## 将行名从 Entrez GeneID 转换为 Gene Symbol
gene_ids <- as.integer(rownames(counts_sub))

## 匹配注释
idx <- match(gene_ids, gene_map$GeneID)
matched <- !is.na(idx)

cat("GeneID 匹配率:", sum(matched), "/", length(gene_ids),
    "(", round(sum(matched)/length(gene_ids)*100, 1), "%)\n")

## 仅保留匹配到 Symbol 的基因
counts_mapped <- counts_sub[matched, ]
symbols       <- gene_map$Symbol[idx[matched]]
gene_types    <- gene_map$GeneType[idx[matched]]

## 去除 Symbol 为空或 NA 的行
valid <- !is.na(symbols) & nchar(symbols) > 0
counts_mapped <- counts_mapped[valid, ]
symbols       <- symbols[valid]
gene_types    <- gene_types[valid]

## 仅保留 protein-coding 基因（可选，推荐用于 WGCNA / ML）
counts_pc <- counts_mapped[gene_types == "protein-coding", ]
symbols_pc <- symbols[gene_types == "protein-coding"]

cat("Protein-coding 基因:", nrow(counts_pc), "\n")

## 处理 Symbol 重复：保留行均值最大的探针
dup_idx <- duplicated(symbols_pc)
if (any(dup_idx)) {
  cat("重复 Symbol:", sum(dup_idx), "个（取最高表达保留）\n")
  df_temp <- data.frame(symbol = symbols_pc,
                        mean_expr = rowMeans(counts_pc),
                        row_idx = seq_along(symbols_pc))
  df_keep <- df_temp %>%
    group_by(symbol) %>%
    slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
    ungroup()
  counts_pc  <- counts_pc[df_keep$row_idx, ]
  symbols_pc <- df_keep$symbol
}

rownames(counts_pc) <- symbols_pc

## ---- 8. 低表达基因过滤 ----
## 保留至少在最小分组样本数中 >= 10 counts 的基因
min_n <- min(table(meta_sub$condition))
keep  <- rowSums(counts_pc >= 10) >= min_n
counts_filtered <- counts_pc[keep, ]

cat("低表达过滤后基因数:", nrow(counts_filtered), "\n")
cat("（要求 >=10 counts 在至少", min_n, "个样本中）\n")

## ---- 9. 构建 DESeqDataSet ----
## 注意：此处仅构建对象，不运行 DESeq()
## DESeq() 应在 fig5.R 的 Step 1 中运行（以保持灵活性）
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = meta_sub,
  design    = ~ condition
)

cat("\nDESeqDataSet 构建完成:",
    nrow(dds), "genes ×", ncol(dds), "samples\n")

## ---- 10. QC 图 ----
pdf("output/qc_plots.pdf", width = 10, height = 8)

## 10.1 Library size 箱线图
lib_size <- data.frame(
  sample    = colnames(counts_filtered),
  total     = colSums(counts_filtered) / 1e6,     # million reads
  condition = meta_sub$condition
)

p_lib <- ggplot(lib_size, aes(x = condition, y = total, fill = condition)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, width = 0.5) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  scale_fill_manual(values = c("Hlty" = "#4DBBD5", "Seps_P" = "#E64B35")) +
  labs(x = NULL, y = "Library size (million filtered reads)",
       title = "Library size distribution") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
print(p_lib)

## 10.2 PCA（VST 变换后）
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Hlty" = "#4DBBD5", "Seps_P" = "#E64B35")) +
  labs(x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)"),
       title = "PCA — VST-transformed expression",
       color = "Condition") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")
print(p_pca)

## 10.3 样本距离热图
sample_dists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sample_dists)

annot_col <- data.frame(
  Condition = meta_sub$condition,
  Sex       = meta_sub$sex,
  row.names = rownames(meta_sub)
)
annot_colors <- list(
  Condition = c(Hlty = "#4DBBD5", Seps_P = "#E64B35"),
  Sex       = c(M = "#3C5488", F = "#F39B7F")
)

pheatmap(
  dist_mat,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(50),
  annotation_col  = annot_col,
  annotation_colors = annot_colors,
  show_rownames   = FALSE,
  show_colnames   = FALSE,
  main = "Sample-to-sample distance (VST)"
)

dev.off()
cat("\nQC 图已保存: output/qc_plots.pdf\n")

## ---- 11. 保存结果对象 ----
saveRDS(counts_filtered, "output/counts_filtered.rds")
saveRDS(meta_sub,        "output/sample_meta.rds")
saveRDS(dds,             "output/dds.rds")
saveRDS(gene_map,        "output/gene_annot.rds")
