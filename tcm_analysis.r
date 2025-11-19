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
ppi <- compute_nodeinfo(ppi)

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
