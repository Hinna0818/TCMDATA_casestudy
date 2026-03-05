library(TCMDATA)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ivolcano)
library(aplotExtra)
library(igraph)
library(ggtangle)
library(ggrepel)
library(aplot)

dir.create("fig2", showWarnings = FALSE)

theme_sci <- function(base_size = 9, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0, size = base_size + 2),
      axis.title       = element_text(face = "bold", size = base_size + 1),
      axis.text        = element_text(color = "black", size = base_size),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.35),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4),
      legend.title     = element_text(face = "bold", size = base_size),
      legend.text      = element_text(size = base_size),
      legend.key       = element_blank(),
      legend.position  = "right",
      plot.margin      = margin(6, 6, 6, 6)
    )
}

## ---------- fig2a: Herb ORA dotplot ----------
data("dn_gcds")
enrich_res <- herb_enricher(genes = dn_gcds, type = "Herb_pinyin_name")

p1 <- dotplot(enrich_res, showCategory = 15) +
  labs(x = "Gene Ratio", y = NULL) +
  theme_sci() +
  theme(axis.text.y = element_text(hjust = 1),
        plot.margin = margin(6, 10, 6, 10))

ggsave("output/fig2a_dotplot.pdf", p1,
       units = "mm", dpi = 600, device = cairo_pdf)

## ---------- fig2b: Volcano plot ----------
data(deg_earlydn)

p2 <- ivolcano(deg_earlydn,
               logFC_col = "log2FoldChange", pval_col = "padj",
               gene_col = "names",
               pval_cutoff = 0.05, pval_cutoff2 = 0.01,
               logFC_cutoff = 1, logFC_cutoff2 = 2,
               size_by = "manual", top_n = 10,
               interactive = FALSE, title = NULL)

ggsave("output/fig2b_volcano.pdf", p2,
       units = "mm", dpi = 600, device = cairo_pdf)

## ---------- fig2c: Venn diagram ----------
degs <- deg_earlydn$names[deg_earlydn$g != "normal"]
data(dn_gcds)
data(dn_otp)
dn_ctd <- read.csv("CTD_DN.csv")
dn_ctd <- subset(dn_ctd, dn_ctd$Inference.Score >= 20)
dn_ctd <- dn_ctd$Gene.Symbol

venn_df <- getvenndata(degs, dn_gcds, dn_otp, dn_ctd,
                       set_names = c("DEGs", "GeneCards", "OpenTargets", "CTD"))

p3 <- ggvenn_plot(venn_df,
                  set.color = c("#FF8748", "#5BAA56", "#B8BB5B", "#6C7BD9"),
                  stroke.color = "white")

ggsave("output/fig2c_venn.pdf", p3,
       units = "mm", dpi = 600, device = cairo_pdf)

## ---------- fig2d: UpSet plot ----------
gene_list <- list(
  DEGs        = degs,
  GeneCards   = dn_gcds,
  OpenTargets = dn_otp,
  CTD         = dn_ctd
)

p4 <- upset_plot(gene_list, color.intersect.by = "Set2", color.set.by = "Dark2")

ggsave("output/fig2d_upset.pdf", p4,
       units = "mm", dpi = 600, device = cairo_pdf)

## ---------- fig2e: Sankey diagram ----------
herbs <- c("灵芝", "黄芪")
df <- search_herb(herb = herbs, type = "Herb_cn_name")

top_mol <- df %>%
  count(molecule, sort = TRUE) %>%
  head(5) %>%
  pull(molecule)

df_sub <- df %>% filter(molecule %in% top_mol)

set.seed(2026)
sampled_targets <- df_sub %>%
  distinct(target) %>%
  sample_n(min(30, n())) %>%
  pull(target)

df_sankey <- df_sub %>% filter(target %in% sampled_targets)

p5 <- tcm_sankey(df_sankey, font_face = NULL)

ggsave("output/fig2e_sankey.pdf", p5,
       units = "mm", dpi = 600, device = cairo_pdf)

## ---------- fig2f: Network graph ----------
g <- prepare_herb_graph(df, n = 60, compute_metrics = TRUE)

node_colors <- c("Herb" = "#C44E52", "Molecule" = "#4C72B0", "Target" = "#55A868")
node_shapes <- c("Herb" = 18, "Molecule" = 16, "Target" = 15)

set.seed(42)
p6 <- ggplot(g, layout = "kk") +
  geom_edge(alpha = 0.15, color = "grey55", linewidth = 0.3) +
  geom_point(aes(color = type, size = degree, shape = type), alpha = 0.85) +
  geom_text_repel(
    aes(label = name, color = type),
    size = 2.5, fontface = "italic",
    max.overlaps = 100, segment.alpha = 0.25,
    segment.size = 0.2, segment.color = "grey60",
    box.padding = 0.3, point.padding = 0.2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = node_colors, name = "Node Type") +
  scale_shape_manual(values = node_shapes, name = "Node Type") +
  scale_size_continuous(name = "Degree", range = c(2, 11)) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)),
    shape = guide_legend(order = 1),
    size  = guide_legend(order = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 10),
    legend.title    = element_text(size = 11, face = "bold"),
    plot.margin     = margin(10, 10, 10, 10)
  )

ggsave("output/fig2f_network.pdf", p6,
       units = "mm", dpi = 600, device = cairo_pdf)


