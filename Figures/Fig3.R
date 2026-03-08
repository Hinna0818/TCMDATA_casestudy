## fig3a
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

theme_sci <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0, size = base_size),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(color = "black", size = base_size),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size),
      legend.key = element_blank(),
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = "grey75", linewidth = 0.3),
      strip.text = element_text(face = "bold", color = "black"),
      plot.margin = margin(6, 8, 6, 8)
    )
}

wrap_labels <- function(width = 32) {
  function(x) vapply(x, function(y) paste(strwrap(y, width = width), collapse = "\n"), character(1))
}

kegg_cols <- c("#B2182B", "#F4A582", "#92C5DE", "#2166AC")
go_cols <- c("#8EC5FC", "#E0C3FC", "#F9F586", "#F08A5D")
gsea_cols <- c("#C03A4A", "#2B6CB0", "#1F8A70")

## kegg preparation
set.seed(2026)
lz_targets <- search_herb("lingzhi", type = "Herb_pinyin_name")$target
query_genes <- sample(lz_targets, 200)

gene_id <- bitr(query_genes,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

kegg_enrich <- enrichKEGG(gene     = gene_id$ENTREZID,
                          organism = "hsa",
                          keyType  = "kegg")

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

p3a <- barplot(kegg_enrich, showCategory=15) + 
  set_enrichplot_color(type = "fill", transform = "log10", colors = kegg_cols) +
  labs(x = "Gene count", y = NULL, fill = expression(-log[10](adjusted~P))) +
  scale_y_discrete(labels = wrap_labels(34)) +
  theme_sci(base_size = 9) +
  theme(
    axis.text.y = element_text(hjust = 1),
    legend.position = "right",
    plot.margin = margin(6, 10, 6, 10)
  )

p3a

ggsave("output/fig3a_barplot.pdf", p3a,
       units = "mm", dpi = 600, device = cairo_pdf)


## fig3b
p3b <- ggdot_sankey(kegg_enrich, n = 8, dot_x_var = "GeneRatio") 
p3b

ggsave("output/fig3b_dotsankey.pdf", p3b, height = 280,
       units = "mm", dpi = 600, device = cairo_pdf)

## GO preparation
# BP
go_bp <- enrichGO(gene     = gene_id$ENTREZID,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP",
                  keyType  = "ENTREZID",
                  readable = TRUE)

# ALL
go_all <- enrichGO(gene = gene_id$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont = "ALL", keyType = "ENTREZID", readable = TRUE)


## fig3c
p3c <- gglollipop(go_bp, line.col = "skyblue", 
                         line.type = "solid", 
                         palette = "RdYlBu", top_n = 10) + 
  scale_y_discrete(labels = wrap_labels(34)) +
  theme_sci(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8.2)
  )

p3c
ggsave("output/fig3c_lollipop_bp.pdf", p3c,
       units = "mm", dpi = 600, device = cairo_pdf)

## fig3d
p3d <- go_barplot(go_all)
p3d
ggsave("output/fig3d_gobarplot.pdf", p3d, height = 250,
       units = "mm", dpi = 600, device = cairo_pdf)

## fig3e
p3e <- gocircle_plot(go_all, top = 8) 
p3e

## GSEA preparation
data(deg_earlydn)

## convert gene symbols to Entrez IDs
deg_ids <- bitr(deg_earlydn$names,
                fromType  = "SYMBOL",
                toType    = "ENTREZID",
                OrgDb     = org.Hs.eg.db)

## merge with DESeq2 results and build ranked list
deg_merged <- merge(deg_earlydn, deg_ids, by.x = "names", by.y = "SYMBOL")
gene_rank  <- sort(setNames(deg_merged$log2FoldChange, deg_merged$ENTREZID),
                   decreasing = TRUE)

gsea_kegg <- gseKEGG(geneList = gene_rank,
                     organism = "hsa",
                     eps      = 0,
                     pvalueCutoff = 0.05)

gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

p3f <- gseaplot2(
  gsea_kegg,
  geneSetID = 1:3,
  pvalue_table = TRUE,
  ES_geom = "line",
  color = gsea_cols,
  base_size = 10,
  rel_heights = c(1.6, 0.45, 1)
)

p3f
ggsave("output/fig3f_gseaplot.pdf", p3f,
       units = "mm", dpi = 600, device = cairo_pdf)
