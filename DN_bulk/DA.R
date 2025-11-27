#' Differential analysis for GSE142025
library(DESeq2)
library(dplyr)
library(ivolcano)
library(enrichplot)
library(TCMDATA)  # remotes::install_github("Hinna0818/TCMDATA")

load("../data/GSE142025/step1-output.Rdata")
cts <- symbol_matrix
cts[1:4,1:4]
dim(cts)
colnames(cts)

# add group info
coldata <- data.frame(condition = group_list, row.names = colnames(cts))
head(coldata)
table(coldata$condition)
coldata$condition <- factor(coldata$condition, levels = c("Control", "Early_DN","Advanced_DN"))
str(coldata)

# DA
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Pre-filtering
min_group_size <- min(table(coldata$condition))
keep <- rowSums(counts(dds) >= 10) >= min_group_size
table(keep)  
dds <- dds[keep, ]

## diff analysis
dds <- DESeq(dds)
res <- results(dds)
res

# "Early_DN" vs "Control"
res1 <- results(dds, contrast = c("condition","Early_DN","Control"))
res1$names <- rownames(res1)
res1$g<-"normal"
res1$g[res1$log2FoldChange > log2(1.5) & res1$padj< 0.05] <- "up"
res1$g[res1$log2FoldChange < -log2(1.5) & res1$padj< 0.05] <- "down"
head(res1)
table(res1$g)
df1 <- as.data.frame(res1)
df1 <- df1[!is.na(df1$names), ]

p1 <- ivolcano(data = df1,
         logFC_col = "log2FoldChange",
         pval_col = "padj",
         gene_col = "names",
         top_n = 10,
         onclick_fun=onclick_genecards,
         title = "Early_DN vs Control")
p1


# "Advanced_DN" vs "Control"
res2 <- results(dds, contrast = c("condition","Advanced_DN", "Control"))
res2$names <- rownames(res2)
res2$g <-"normal"
res2$g[res2$log2FoldChange > log2(1.5) & res2$padj< 0.05] <-"up"
res2$g[res2$log2FoldChange < -log2(1.5) & res2$padj< 0.05] <-"down"
head(res2)
table(res2$g)
df2 <- as.data.frame(res2)
df2 <- df2[!is.na(df2$names), ]

p2 <- ivolcano(data = df2,
         logFC_col = "log2FoldChange",
         pval_col = "padj",
         gene_col = "names",
         top_n = 10,
         onclick_fun=onclick_genecards,
         title = "Advanced_DN vs Control")
p2


# "Advanced_DN" vs "Early_DN"
res3 <- results(dds, contrast = c("condition","Advanced_DN", "Early_DN"))
res3$names <- rownames(res3)
res3$g <-"normal"
res3$g[res3$log2FoldChange > log2(1.5) & res3$padj< 0.05] <-"up"
res3$g[res3$log2FoldChange < -log2(1.5) & res3$padj< 0.05] <-"down"
head(res3)
table(res3$g)
df3 <- as.data.frame(res3)
df3 <- df3[!is.na(df3$names), ]

p3 <- ivolcano(data = df3,
         logFC_col = "log2FoldChange",
         pval_col = "padj",
         gene_col = "names",
         top_n = 10,
         onclick_fun=onclick_genecards,
         title = "Advanced_DN vs Early_DN")
p3


## save DEGs for three contrast
df1_sig <- df1 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange)))

df2_sig <- df2 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange)))

df3_sig <- df3 %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange)))

deg_all_sig <- list(df1_sig = df1_sig, df2_sig = df2_sig, df3_sig = df3_sig)
deg_all <- list(df1 = df1, df2 = df2, df3 = df3)

dir.create("../data/GSE142025_DEGs")
save(deg_all_sig, file = "../data/GSE142025_DEGs/deg_all_sig.rda")
save(deg_all, file = "../data/GSE142025_DEGs/deg_all.rda")

## vennplot for DEGs
# 3 types of DEGs overlap
vd <- getvenndata(df1_sig$names, df2_sig$names, 
                  df3_sig$names, set_names = c("Early vs Control", "Advanced vs Control", "Advanced vs Early"))
p3 <- ggvenn_plot(vd)
p3

# early-control up ∩ advanced-early down (protection)
ec_up <- subset(df1_sig, subset = df1_sig$g == "up")$names
ae_down <- subset(df3_sig, subset = df3_sig$g == "down")$names
vd_protect <- getvenndata(ec_up, ae_down, set_names = c("Early_Control up", "Advanced_Early down"))
p4 <- ggvenn_plot(vd_protect)
p4


# early-control down ∩ advanced-early up (progress)
ec_down <- subset(df1_sig, subset = df1_sig$g == "down")$names
ae_up <- subset(df3_sig, subset = df3_sig$g == "up")$names
vd_progress <- getvenndata(ec_down, ae_up, set_names = c("Early_Control down", "Advanced_Early up"))
p5 <- ggvenn_plot(vd_progress)
p5

## find target herb for the DN progress
progress_g <- intersect(ec_down, ae_up)
progress_res <- herb_enricher(genes = progress_g, pvalueCutoff = 0.05, qvalueCutoff = 0.2, type = "Herb_cn_name")
p6 <- dotplot(progress_res)
p6  ## dotplot shows quanxie can be one of the herbs related to DN


# early-control up ∩ advanced-control up (linear mounting) 
ec_up <- subset(df1_sig, subset = df1_sig$g == "up")$names
ac_up <- subset(df2_sig, subset = df2_sig$g == "up")$names
vd_mounting <- getvenndata(ec_up, ac_up, set_names = c("Early_Control up", "Advanced_Control up"))
p7 <- ggvenn_plot(vd_mounting)
p7

