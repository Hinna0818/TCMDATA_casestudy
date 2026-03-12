## ============================================================
## 5tcmanalysis.R — TCMDATA 中药网络药理学整合分析

rm(list = ls())

library(TCMDATA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(igraph)
library(ggplot2)
library(ggrepel)
library(ggtangle)
library(dplyr)

## ---- Lancet 学术主题 ----
theme_sci <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title    = element_text(face = "bold", size = base_size + 1),
      axis.text     = element_text(color = "black", size = base_size),
      axis.line     = element_line(color = "black", linewidth = 0.4),
      axis.ticks    = element_line(color = "black", linewidth = 0.35),
      legend.title  = element_text(face = "bold", size = base_size),
      legend.text   = element_text(size = base_size),
      legend.key    = element_blank(),
      plot.margin   = margin(4, 4, 4, 4)
    )
}

pal_lancet <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                "#F39B7F", "#8491B4", "#91D1C2", "#DC9A6C")
pal_edge   <- "grey65"

## 1. load data
intersection_genes <- readRDS("output/nodes_after_intersection.rds")

## 使用交集分析后的基因进行富集和 PPI
hub_targets <- unique(intersection_genes)

## fig5d中药-单体-靶点的桑基图
lingzhi <- search_herb(herb = "lingzhi", type = "Herb_pinyin_name")
sankey_data <- subset(lingzhi, lingzhi$target %in% hub_targets)

p5d <- tcm_sankey(sankey_data, font_face = NULL)
p5d

ggsave("output/fig5d_sankey.pdf", p5d,
       width = 200, height = 200, units = "mm", dpi = 600,
       device = cairo_pdf)

## ############################################################
## PART 2 — 富集分析
## ############################################################

## kegg
gene_id <- bitr(hub_targets,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

kegg_enrich <- enrichKEGG(gene     = gene_id$ENTREZID,
                          organism = "hsa",
                          keyType  = "kegg")

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

library(ggplot2)
p5e <- dotplot(kegg_enrich, showCategory = 30)
p5e
ggsave("output/fig5e_kegg.pdf", p5e, units = "mm", dpi = 600, height = 240,
       device = cairo_pdf)


## go
go_all <- enrichGO(gene = gene_id$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont = "ALL", keyType = "ENTREZID", readable = TRUE)

p5f <- barplot(
  go_all,
  x = "GeneRatio",
  showCategory = 10,
  split = "ONTOLOGY"
) +
  aes(fill = ONTOLOGY) +
  scale_fill_brewer(palette = "Set2") +
  enrichplot::autofacet(by = "row", scales = "free") + 
  scale_y_discrete()

p5f 
ggsave("output/fig5f_go.pdf", p5f, units = "mm", dpi = 600,
       device = cairo_pdf)



## ############################################################
## PART 3 — PPI 网络分析
## ############################################################

## ---- 6. 构建 PPI ----
## 使用 clusterProfiler::getPPI() 从 STRING 获取 PPI
ppi_raw <- getPPI(hub_targets, taxID = 9606)

cat("\nPPI 原始网络:", vcount(ppi_raw), "nodes,", ecount(ppi_raw), "edges\n")

## 过滤：score ≥ 0.5，移除 degree < 2 的孤立节点
ppi_filtered <- ppi_subset(ppi_raw, score_cutoff = 0.5)
low_deg <- V(ppi_filtered)[degree(ppi_filtered) < 2]$name
if (length(low_deg) > 0) {
  ppi_filtered <- delete_vertices(ppi_filtered, low_deg)
}
cat("PPI 过滤后:", vcount(ppi_filtered), "nodes,",
    ecount(ppi_filtered), "edges\n")

## ---- 7. 计算拓扑指标 + Hub 排序 ----
ppi_scored <- compute_nodeinfo(ppi_filtered, weight_attr = "score")
rank_res   <- rank_ppi_nodes(ppi_scored, use_weight = TRUE)
ppi_ranked <- rank_res$graph
rank_df    <- rank_res$table


V(ppi_ranked)$btw_raw <- igraph::betweenness(ppi_ranked, directed = FALSE,
                                              normalized = FALSE)


## ---- 8. Fig 5G: PPI Hub 网络 ----
set.seed(42)
p5g <- ggplot(ppi_ranked, layout = "kk") +
  geom_edge(alpha = 0.4, color = pal_edge, linewidth = 0.45) +
  geom_point(aes(size = btw_raw, color = degree), alpha = 0.9) +
  scale_color_gradientn(
    colours = c("#3C5488", "#4DBBD5", "#F7F7F7", "#F39B7F", "#E64B35"),
    name = "Degree"
  ) +
  scale_size_continuous(range = c(1.5, 7), name = "Betweenness") +
  geom_text_repel(
    aes(label = name),
    fontface = "italic",
    size = 2.2,
    max.overlaps = 25, segment.alpha = 0.4,
    segment.size = 0.25, box.padding = 0.35, point.padding = 0.2
  ) +
  theme_void() +
  theme(legend.title = element_text(face = "bold", size = 9),
        legend.text  = element_text(size = 8))

p5g
ggsave("output/fig5g_ppi_hub.pdf", p5g,
       units = "mm", dpi = 600,
       device = cairo_pdf)

## ---- 保存结果 ----
saveRDS(ppi_scored,  "output/ppi_scored.rds")
saveRDS(rank_res,    "output/ppi_rank.rds")


