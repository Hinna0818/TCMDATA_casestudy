## fig4
rm(list = ls())
library(clusterProfiler)
library(TCMDATA)
library(ggtangle)
library(igraph)
library(aplot)
library(ggplot2)
library(ggrepel)
library(grid)
library(patchwork)

## ---- Academic theme ----
theme_sci <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text  = element_text(color = "black", size = base_size),
      axis.line  = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text  = element_text(size = base_size),
      legend.key   = element_blank(),
      plot.margin  = margin(4, 4, 4, 4)
    )
}

## Color palette
pal_node   <- "#4C72B0"
pal_edge   <- "grey65"
pal_clust  <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                "#F39B7F", "#8491B4", "#91D1C2", "#DC9A6C")

## ============================================================
## Prepare PPI demo data
## ============================================================
data("demo_ppi")
ppi_filtered <- ppi_subset(demo_ppi, score_cutoff = 0.7)

## Remove low-degree nodes (degree < 2) to reduce noise
low_deg <- V(ppi_filtered)[degree(ppi_filtered) < 2]$name
if (length(low_deg) > 0) {
  ppi_filtered <- delete_vertices(ppi_filtered, low_deg)
  message("Removed ", length(low_deg), " low-degree node(s): ", paste(low_deg, collapse = ", "))
}
cat("After filtering: ", vcount(ppi_filtered), " nodes, ", ecount(ppi_filtered), " edges\n")

ppi_scored   <- compute_nodeinfo(ppi_filtered, weight_attr = "score")
rank_res     <- rank_ppi_nodes(ppi_scored, use_weight = TRUE)
ppi_ranked   <- rank_res$graph
rank_df      <- rank_res$table
top_nodes    <- rank_df$name[1:2]   # INS, TNF

## Raw (un-normalized) betweenness for size mapping — better visual spread
V(ppi_ranked)$btw_raw <- igraph::betweenness(ppi_ranked, directed = FALSE, normalized = FALSE)

## ============================================================
## Fig 4A — Four layout algorithms comparison
## ============================================================
base_layers <- list(
  geom_edge(alpha = 0.25, color = pal_edge, linewidth = 0.55),
  geom_point(aes(size = btw_raw, color = degree), alpha = 0.9),
  scale_color_distiller(palette = "Spectral", direction = -1, name = "Degree"),
  scale_size_continuous(range = c(1.5, 7), name = "Betweenness"),
  geom_text_repel(
    aes(label = ifelse(degree >= quantile(degree, 0.85), name, "")),
    size = 2.2, fontface = "italic",
    max.overlaps = 20, segment.alpha = 0.4,
    segment.size = 0.25, box.padding = 0.3, point.padding = 0.2
  ),
  theme_void(),
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.title = element_text(face = "bold", size = 9))
)

set.seed(42)
p_fr <- ggplot(ppi_ranked, layout = "fr") + base_layers +
  ggtitle("Fruchterman\u2013Reingold")
p_fr

p_kk <- ggplot(ppi_ranked, layout = "kk") + base_layers +
  ggtitle("Kamada\u2013Kawai")

p_circle <- ggplot(ppi_ranked, layout = "circle") + base_layers +
  ggtitle("Circle")

p_dh <- ggplot(ppi_ranked, layout = "dh") + base_layers +
  ggtitle("Davidson\u2013Harel")

p4a <- plot_list(p_fr, p_kk, p_circle, p_dh, ncol = 2) 

p4a
ggsave("output/figure4/fig4a_layouts.pdf", p4a,
       width = 220, height = 200, units = "mm", dpi = 600,
       device = cairo_pdf)

## ============================================================
## Fig 4B — Radar plots for INS and TNF
## ============================================================
p_radar1 <- radar_plot(
  get_node_profile(rank_df, top_nodes[1]),
  fill_color = "#4DBBD5", line_color = "#00768A",
  title = top_nodes[1]
) + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 12))

p_radar2 <- radar_plot(
  get_node_profile(rank_df, top_nodes[2]),
  fill_color = "#E64B35", line_color = "#B83226",
  title = top_nodes[2]
) + theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 12))

p4b <- plot_list(p_radar1, p_radar2, ncol = 2)
p4b
ggsave("output/figure4/fig4b_radarplot.pdf", p4b,
       width = 180, height = 90, units = "mm", dpi = 600,
       device = cairo_pdf)

## ============================================================
## Fig 4C — Z-score heatmap of all nodes
## ============================================================
selected_cols <- c("degree", "betweenness", "closeness", "MCC", "MNC",
                   "DMNC", "EPC", "radiality", "Stress")

p4c <- plot_node_heatmap(
  rank_df,
  select_cols = selected_cols,
  colors = c("#3C5488", "white", "#E64B35"))
p4c
pdf("output/figure4/fig4c_heatmap.pdf", width = 6, height = 8)
draw(p4c)
dev.off()

## ============================================================
## Fig 4D — Network with node color/size mapped to composite score
## ============================================================
set.seed(42)
p4d <- ggplot(ppi_ranked, layout = "fr") +
  geom_edge(alpha = 0.25, color = "grey70", linewidth = 0.55) +
  geom_point(
    aes(size = btw_raw, color = Score_network),
    alpha = 0.9
  ) +
  scale_color_gradientn(
    colours = c("#3C5488", "#4DBBD5", "#F7F7F7", "#F39B7F", "#E64B35"),
    name = "Node Score",
    limits = c(0, 1)
  ) +
  scale_size_continuous(range = c(2, 10), name = "Betweenness") +
  guides(color = guide_colorbar(order = 1),
         size  = guide_legend(order = 2)) +
  geom_text_repel(
    aes(label = name),
    size = 2.5, fontface = "italic",
    max.overlaps = 50, segment.alpha = 0.35,
    segment.size = 0.2, segment.color = "grey50",
    box.padding = 0.3, point.padding = 0.2
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(6, 6, 6, 6)
  )

p4d
ggsave("output/figure4/fig4d_score_network.pdf", p4d,
       width = 160, height = 140, units = "mm", dpi = 600,
       device = cairo_pdf)

## ============================================================
## Fig 4E — Community detection (MCODE)
## ============================================================
ppi_mcode <- run_mcode(ppi_ranked)
mcode_df  <- get_mcode_res(ppi_mcode, only_clusters = FALSE)
print(mcode_df[!is.na(mcode_df$cluster), c("name", "cluster", "module_score", "is_seed")])

set.seed(42)
p4e <- ggplot(ppi_mcode, layout = "fr") +
  geom_edge(alpha = 0.25, color = "grey70", linewidth = 0.55) +
  geom_point(
    aes(color = mcode_cluster, size = degree),
    alpha = 0.9
  ) +
  scale_color_manual(
    values = pal_clust,
    name = "Module",
    na.value = "grey80"
  ) +
  scale_size_continuous(range = c(2, 9), name = "Degree") +
  geom_text_repel(
    aes(label = ifelse(degree >= quantile(degree, 0.80), name, "")),
    size = 2.5, fontface = "italic",
    max.overlaps = 25, segment.alpha = 0.4,
    segment.size = 0.25, box.padding = 0.35, point.padding = 0.25
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(6, 6, 6, 6)
  )

p4e
ggsave("output/figure4/fig4e_community.pdf", p4e,
       width = 160, height = 140, units = "mm", dpi = 600,
       device = cairo_pdf)

## ============================================================
## Fig 4F — PPI network after IL6 knockout
## ============================================================
ko_res <- ppi_knock(ppi_ranked, targets = "IL6", n_perm = 100, seed = 42)
print(ko_res$Summary)
cat("Total disruption score:", ko_res$Total_Score, "\n")
cat("Total P-value:", ko_res$Total_Pvalue, "\n")

## Build knocked-out igraph (remove IL6 and re-score)
ppi_ko <- delete_vertices(ppi_ranked, "IL6")

set.seed(42)
## Compute raw betweenness for knocked-out network
V(ppi_ko)$btw_raw <- igraph::betweenness(ppi_ko, directed = FALSE, normalized = FALSE)

p4f <- ggplot(ppi_ko, layout = "fr") +
  geom_edge(alpha = 0.25, color = "grey70", linewidth = 0.55) +
  geom_point(
    aes(size = btw_raw, color = degree),
    alpha = 0.9
  ) +
  scale_color_distiller(palette = "Spectral", direction = -1,
                        name = "Degree") +
  scale_size_continuous(range = c(2, 9), name = "Betweenness") +
  guides(color = guide_colorbar(order = 1),
         size  = guide_legend(order = 2)) +
  geom_text_repel(
    aes(label = name),
    size = 2.3, fontface = "italic",
    max.overlaps = 50, segment.alpha = 0.35,
    segment.size = 0.2, segment.color = "grey50",
    box.padding = 0.3, point.padding = 0.2
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = "IL6 knocked out", fontface = "bold.italic",
           size = 3.5, color = "#B2182B") +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(6, 6, 6, 6)
  )

ggsave("output/figure4/fig4f_knockout.pdf", p4f,
       width = 160, height = 140, units = "mm", dpi = 600,
       device = cairo_pdf)

