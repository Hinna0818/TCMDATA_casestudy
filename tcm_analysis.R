#' TCM analysis of DN progression

library(TCMDATA)
library(ggalluvial)
library(clusterProfiler)
library(dplyr)

load("../data/GSE142025_DEGs/deg_all_sig.rda")
df1_sig <- deg_all_sig[[1]]
df3_sig <- deg_all_sig[[3]]
ec_down <- subset(df1_sig, subset = df1_sig$g == "down")$names
ae_up <- subset(df3_sig, subset = df3_sig$g == "up")$names
progress_g <- intersect(ec_down, ae_up)

## search DN targets in GeneCards(https://www.genecards.org/)
DN_gc <- read.csv("../data/DN_GeneCards.csv") %>%
  select(Gene.Symbol, Relevance.score) %>%
  filter(Relevance.score > 1) %>%
  rename(target = Gene.Symbol)

## show overlap between quanxie targets and progression genes and DN targets in GeneCards
qx <- search_herb(herb = "quanxie", type = "Herb_pinyin_name")
a <- getvenndata(qx$target, progress_g, DN_gc$target, set_names = c("quanxie targets", 
                                                             "progression genes",
                                                             "Diabetic Nephropathy targets"))
ggvenn_plot(a)

## intersect qx targets, progress genes, and DN targets
qx_new <- subset(qx, subset = qx$target %in% intersect(progress_g, DN_gc$target))
ts <- TCM_sankey(qx_new)
ggsave("quanxie_sankey.pdf", ts, width = 12, height = 14)


## GO enrichment analysis for qx_new target
qx_g <- enrichGO(qx_new$target, OrgDb='org.Hs.eg.db', keyType="SYMBOL", ont = "BP")
qx_lolli <- gglollipop(qx_g, show_count = F)
qx_lolli

gs <- TCMDATA:::ggdot_sankey(qx_g)
ggsave("go_sankey.pdf", gs, width = 12, height = 14)

## PPI analysis
ppi <- getPPI(qx_new$target, taxID = "9606")

# PPI subset
ppi <- ppi_subset(ppi, score_cutoff = 0.7)
ppi <- compute_nodeinfo(ppi)

## PPI network
library(ggrepel)
library(ggplot2)
library(ggtangle)
library(grid)   

set.seed(42)
p1 <- ggtangle::ggplot(ppi, layout = "kk") +
  geom_edge(
    colour  = "grey80",
    alpha   = 0.75,
    linewidth = 0.15,
    lineend  = "round"
  ) +
  geom_point(
    aes(colour = degree, size = betweenness),
    alpha = 0.85
  ) +
  geom_text_repel(
    aes(label = name, colour = degree),
    size          = 2.7,
    fontface      = "bold",     
    force         = 1.8,
    max.overlaps  = 50,
    box.padding   = unit(0.2, "lines"),
    point.padding = unit(0.15, "lines"),
    segment.size  = 0.1,
    segment.color = "grey80",
    show.legend   = FALSE
  ) +
  scale_color_gradientn(
    name    = "Degree",
    colours = c("#4DBBD5FF", "#E64B35FF")
  ) +
  scale_size_continuous(
    name  = "Betweenness",
    range = c(1.5, 6)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 9),
    plot.margin     = unit(c(0.6, 0.8, 0.6, 0.8), "cm"),
    plot.background = element_rect(fill = "white", colour = NA),
    text            = element_text(face = "bold"))

print(p1)

## delete two outlier nodes
nodes_to_remove <- c("TOP2A", "HMMR")
ppi <- delete_vertices(ppi, nodes_to_remove)

p2 <- ggtangle::ggplot(ppi, layout = "kk") +
  geom_edge(
    colour  = "grey80",
    alpha   = 0.75,
    linewidth = 0.15,
    lineend  = "round"
  ) +
  geom_point(
    aes(colour = degree, size = betweenness),
    alpha = 0.85
  ) +
  geom_text_repel(
    aes(label = name, colour = degree),
    size          = 3.5,
    fontface      = "bold",     
    force         = 1.8,
    max.overlaps  = 50,
    box.padding   = unit(0.2, "lines"),
    point.padding = unit(0.15, "lines"),
    segment.size  = 0.1,
    segment.color = "grey80",
    show.legend   = FALSE
  ) +
  scale_color_gradientn(
    name    = "Degree",
    colours = c("#4DBBD5FF", "#E64B35FF")
  ) +
  scale_size_continuous(
    name  = "Betweenness",
    range = c(1.5, 6)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 9),
    plot.margin     = unit(c(0.6, 0.8, 0.6, 0.8), "cm"),
    plot.background = element_rect(fill = "white", colour = NA),
    text            = element_text(face = "bold"))

print(p2)

## rank hub nodes
ppi_res <- rank_ppi_nodes(ppi)
ppi_new <- ppi_res[["graph"]]
ppi_res <- ppi_res[["table"]]


## radar plot for hub nodes
top_nodes <- ppi_res$name |> head(5)
node_info <- get_node_profile(ppi_res, node_name = top_nodes[3])
p3 <- radar_plot(node_info, title = "Centrality profile of IL1B")
p3

node_info2 <- get_node_profile(ppi_res, node_name = top_nodes[5])
p4 <- radar_plot(node_info2, title = "Centrality profile of JUN", fill_color = "#D59390", line_color = "#D59380")
p4

aplot::plot_list(p3, p4)


## MCL clustering (optional)
ppi_new1 <- run_MCL(ppi_new, inflation = 2.5)

## louvain clustering (recommended)
set.seed(42)
ppi_new1 <- run_louvain(ppi_new1, resolution = 0.6)
p5 <- ggraph(ppi_new1, layout = 'kk') +
  geom_edge_link() +
  geom_point_interactive(
    mapping = aes(
      x = x, 
      y = y,
      color = louvain_cluster,
      size = betweenness,
      tooltip = name,
      data_id = name
    )
  ) +
  geom_text_repel(data = td_filter(degree > 60),
                  mapping = aes(x = x, y = y, label = name), bg.color ='white'
  ) +
  theme_graph()

girafe(ggobj = p5)





## subgraph for top nodes
top_nodes <- head(ppi_res$name, 10)
subg <- igraph::induced_subgraph(ppi_new, vids = top_nodes)
plot(subg)


p3 <- ggtangle::ggplot(subg, layout = "kk") +
  geom_edge(
    colour  = "grey80",
    alpha   = 0.75,
    linewidth = 0.15,
    lineend  = "round"
  ) +
  geom_point(
    aes(colour = Score_network, size = betweenness),
    alpha = 0.85
  ) +
  geom_text_repel(
    aes(label = name, colour = degree),
    size          = 3.5,
    fontface      = "bold",     
    force         = 1.8,
    max.overlaps  = 50,
    box.padding   = unit(0.2, "lines"),
    point.padding = unit(0.15, "lines"),
    segment.size  = 0.1,
    segment.color = "grey80",
    show.legend   = FALSE
  ) +
  scale_color_gradientn(
    name    = "Score_network",
    colours = c("#4DBBD5FF", "#E64B35FF")
  ) +
  scale_size_continuous(
    name  = "Betweenness",
    range = c(1.5, 6)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 9),
    plot.margin     = unit(c(0.6, 0.8, 0.6, 0.8), "cm"),
    plot.background = element_rect(fill = "white", colour = NA),
    text            = element_text(face = "bold"))

print(p3)



## top nodes enrichment analysis
topnodes_go <- enrichGO(top_nodes, OrgDb='org.Hs.eg.db', keyType="SYMBOL", ont = "BP")
tn_p <- gglollipop(topnodes_go, show_count = FALSE)
tn_p


## top nodes PPI networks visualization
p3 <- ggtangle::ggplot(subg, layout = "kk") +
  geom_edge(
    colour  = "grey70",
    alpha   = 0.75,
    linewidth = 0.15,
    lineend  = "round"
  ) +
  geom_point(
    aes(colour = Score_col, size = betweenness),
    alpha = 0.9
  ) +
  geom_text(
    aes(label = name),
    colour      = "black",    
    size        = 3,        
    fontface    = "bold",
    show.legend = FALSE      
  )+
  scale_color_gradientn(
    name    = "node score",
    colours = c("#4DBBD5FF", "#FFDC91", "#E64B35FF"),
    limits  = range(igraph::vertex_attr(subg)$Score_network)
  ) +
  scale_size_continuous(
    name  = "Betweenness",
    range = c(10, 20)   
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 9),
    plot.margin     = unit(c(0.6, 0.8, 0.6, 0.8), "cm"),
    plot.background = element_rect(fill = "white", colour = NA),
    text            = element_text(face = "bold")
  )

print(p3)








