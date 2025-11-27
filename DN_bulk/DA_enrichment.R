#' DEGs enrichment analysis for the progression module
#' 思路：
#' 1. 首先对整体的advanced vs early做差异分析，然后富集到免疫通路；
#' 2. 再去验证进程基因（对early down基因做富集，确认没有免疫相关通路；再对advanced up基因
#' 富集，验证是否免疫通路恢复）；
#' 3. 最后去early down且advanced high的基因做富集。

library(dplyr)
library(TCMDATA)
library(clusterProfiler)

## DEGs of advance_DN and early_DN
load("../data/GSE142025_DEGs/deg_all_sig.rda")
AE_degs <- deg_all_sig[[3]] %>% 
  filter(g != "normal")

## KEGG and GO enrichment analysis of DEGs(advance_DN and early_DN)
g <- AE_degs$names
x <- enrichGO(g, OrgDb='org.Hs.eg.db', keyType="SYMBOL", ont = "BP")
lp <- gglollipop(x)
lp

## select early down and advance up DEGs
df1_sig <- deg_all_sig[[1]]
df3_sig <- deg_all_sig[[3]]
ec_down <- subset(df1_sig, subset = df1_sig$g == "down")$names
ae_up <- subset(df3_sig, subset = df3_sig$g == "up")$names
progress_g <- intersect(ec_down, ae_up)

# early down GO
gg <- enrichGO(ec_down, OrgDb='org.Hs.eg.db', keyType="SYMBOL", ont = "BP")
lpp <- gglollipop(gg)
lpp

# advance up GO
ggg <- enrichGO(ae_up, OrgDb='org.Hs.eg.db', keyType="SYMBOL", ont = "BP")
lppp <- gglollipop(ggg)
lppp

# early down and advance up GO
gggg <- enrichGO(progress_g, OrgDb='org.Hs.eg.db', keyType="SYMBOL", ont = "BP")
lpppp <- gglollipop(gggg)
lpppp














