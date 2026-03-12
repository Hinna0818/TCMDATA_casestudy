## 4intersection
# 下载sepsis疾病靶点，与wgcna模块、degs取交集分析
rm(list = ls())
library(dplyr)
library(TCMDATA)
library(enrichplot)
library(aplot)

## 1. genecards
dir1 <- "E:/Yulab/pharmacology/写作/tcmdata_writing_skill/case_study/sepsis_targets/sepsis_genecards.csv"
genecards <- read.csv(dir1)
genecards <- genecards %>%
  select(Gene.Symbol) %>%
  unlist()

## 2. open target platform
dir2 <- "E:/Yulab/pharmacology/写作/tcmdata_writing_skill/case_study/sepsis_targets/sepsis_otp.tsv"
otp <- read.delim(dir2)
otp <- otp %>%
  select(symbol) %>%
  unlist()

## 3. CTD
dir3 <- "E:/Yulab/pharmacology/写作/tcmdata_writing_skill/case_study/sepsis_targets/sepsis_ctd.csv"
ctd <- read.csv(dir3)
ctd <- ctd %>%
  select(Gene.Symbol) %>%
  unlist()

## intersection with 3 sepsis_targets
sepsis_targets <- intersect(genecards, otp) |> intersect(ctd)

## intersection with sepsis_targets, degs, and wgcna hub genes
degs <- readRDS("./output/sig_degs.rds")$gene
hub_genes <- readRDS("./output/turquoise_module_genes.rds")

hub_nodes <- intersect(sepsis_targets, degs) |> intersect(hub_genes)

## herb ora analysis
herb_enrichres <- herb_enricher(hub_nodes)
dotplot(herb_enrichres)

lingzhi <- search_herb(herb = "lingzhi", type = "Herb_pinyin_name")
final_node <- intersect(hub_nodes, lingzhi$target)


## vennplot
venn_df1 <- getvenndata(genecards, otp, ctd,
                       set_names = c("GeneCards", "Open Targets Platform", "CTD"))
venn1 <- ggvenn_plot(venn_df1)
venn1

venn_df2 <- getvenndata(sepsis_targets, degs, hub_genes, lingzhi$target,
                        set_names = c("sepsis targets", "DEGs", "WGCNA hubs", "lingzhi targets"))
venn2 <- ggvenn_plot(venn_df2)
venn2

ggsave("output/fig5c_venn1.pdf", venn1, units = "mm", dpi = 600, device = cairo_pdf)
ggsave("output/fig5c_venn2.pdf", venn2, units = "mm", dpi = 600, device = cairo_pdf)

saveRDS(final_node, "./output/nodes_after_intersection.rds")



