
rm(list = ls())

## ---- Load Libraries ----
library(DESeq2)
library(apeglm)
library(WGCNA)
library(glmnet)
library(randomForest)
library(TCMDATA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(igraph)
library(ggtangle)
library(ggplot2)
library(ggrepel)
library(aplot)
library(patchwork)
library(dplyr)

allowWGCNAThreads()

dir.create("output/figure5", showWarnings = FALSE, recursive = TRUE)

## ---- Academic theme (consistent with fig2–fig4) ----
theme_sci <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
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

## Lancet-style colour palette (consistent with fig2–fig4)
pal_lancet <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                "#F39B7F", "#8491B4", "#91D1C2", "#DC9A6C")
pal_edge   <- "grey65"


## ############################################################
## STEP 0 — LOAD PRE-PROCESSED DATA (from 1dataloader.R)
## ############################################################
## Pre-processed by case_study/1dataloader.R:
##   - Raw counts (GSE154918_raw_counts_GRCh38.p13_NCBI.tsv)
##   - GEOquery metadata → two-group filtering (Hlty vs Seps_P)
##   - GeneID → Symbol mapping (protein-coding only)
##   - Low-expression gene filtering
## Run 1dataloader.R first if output files do not exist.
## ############################################################

data_dir <- "../case_study/output"
if (!file.exists(file.path(data_dir, "dds.rds"))) {
  stop("Pre-processed data not found. Please run case_study/1dataloader.R first.")
}

counts <- readRDS(file.path(data_dir, "counts_filtered.rds"))
meta   <- readRDS(file.path(data_dir, "sample_meta.rds"))
dds    <- readRDS(file.path(data_dir, "dds.rds"))

cat("Samples retained:", nrow(meta),
    "  Genes after filtering:", nrow(counts), "\n")


## ############################################################
## STEP 1 — DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2)
## ############################################################

## dds already loaded from 1dataloader.R output; run DESeq2
dds <- DESeq(dds)

## LFC shrinkage with apeglm
res <- lfcShrink(dds, coef = "condition_Seps_P_vs_Hlty", type = "apeglm")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

## Significant DEGs
sig_degs <- res_df %>%
  filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) >= 1) %>%
  arrange(padj)

cat("Significant DEGs:", nrow(sig_degs),
    " (Up:", sum(sig_degs$log2FoldChange > 0),
    " Down:", sum(sig_degs$log2FoldChange < 0), ")\n")


## ---- Fig 5A: Volcano plot ----
res_plot <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj <= 0.05 & log2FoldChange >=  1 ~ "Up",
      padj <= 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    ),
    neg_log10p = -log10(pmax(padj, 1e-300))
  )

## Label top 10 up + top 10 down by padj
top_label <- bind_rows(
  res_plot %>% filter(sig == "Up")   %>% slice_min(padj, n = 10),
  res_plot %>% filter(sig == "Down") %>% slice_min(padj, n = 10)
)

p5a <- ggplot(res_plot, aes(log2FoldChange, neg_log10p)) +
  geom_point(aes(color = sig), size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c(Up = "#E64B35", Down = "#3C5488", NS = "grey70"),
                     name = NULL) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  geom_text_repel(data = top_label, aes(label = gene),
                  size = 2.2, fontface = "italic",
                  max.overlaps = 20, segment.size = 0.2) +
  labs(x = expression(log[2]~Fold~Change), y = expression(-log[10]~padj),
       title = "Differential expression (Seps_P vs Hlty)") +
  theme_sci(base_size = 9)

p5a
ggsave("output/figure5/fig5a_volcano.pdf", p5a,
       width = 150, height = 130, units = "mm", dpi = 600,
       device = cairo_pdf)


## ############################################################
## STEP 2 — WGCNA CO-EXPRESSION NETWORK
## ############################################################

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

## Select top 5000 most variable genes
var_genes <- apply(vsd_mat, 1, var)
top5k     <- names(sort(var_genes, decreasing = TRUE))[1:5000]
expr_wgcna <- t(vsd_mat[top5k, ])   # samples × genes

## Soft-threshold selection
powers <- c(1:20)
sft    <- pickSoftThreshold(expr_wgcna, powerVector = powers,
                            networkType = "signed", verbose = 0)
soft_power <- sft$powerEstimate
if (is.na(soft_power)) soft_power <- 12
cat("Soft-threshold power:", soft_power, "\n")

## Build network
set.seed(42)
net <- blockwiseModules(
  expr_wgcna,
  power            = soft_power,
  networkType      = "signed",
  TOMType          = "signed",
  minModuleSize    = 30,
  mergeCutHeight   = 0.25,
  numericLabels    = FALSE,
  saveTOMs         = FALSE,
  verbose          = 0
)

## Module–trait correlation
trait <- data.frame(
  sepsis = as.numeric(meta$condition == "Seps_P"),
  row.names = rownames(meta)
)
MEs <- net$MEs
moduleTraitCor  <- cor(MEs, trait, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(expr_wgcna))

## ---- Fig 5B: Module–trait correlation heatmap ----
## Build a data frame for ggplot-based heatmap
mt_df <- data.frame(
  module = gsub("ME", "", rownames(moduleTraitCor)),
  cor    = moduleTraitCor[, "sepsis"],
  pval   = moduleTraitPval[, "sepsis"],
  label  = paste0(round(moduleTraitCor[, "sepsis"], 2), "\n(",
                  signif(moduleTraitPval[, "sepsis"], 2), ")")
)
mt_df <- mt_df[order(mt_df$cor), ]
mt_df$module <- factor(mt_df$module, levels = mt_df$module)

p5b <- ggplot(mt_df, aes(x = "Sepsis", y = module, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.2) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Correlation") +
  labs(x = NULL, y = NULL, title = "Module–trait correlation") +
  theme_sci(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold"),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

p5b
ggsave("output/figure5/fig5b_module_trait.pdf", p5b,
       width = 100, height = 140, units = "mm", dpi = 600,
       device = cairo_pdf)

## Extract hub genes from significant modules (|r| > 0.5, p < 0.05)
sig_modules <- mt_df$module[abs(mt_df$cor) > 0.5 & mt_df$pval < 0.05]
sig_modules <- as.character(sig_modules)

wgcna_hub <- c()
for (mod in sig_modules) {
  mod_genes <- names(net$colors)[net$colors == mod]
  GS <- abs(cor(expr_wgcna[, mod_genes, drop = FALSE],
                trait$sepsis, use = "p"))
  MM <- abs(cor(expr_wgcna[, mod_genes, drop = FALSE],
                MEs[, paste0("ME", mod)], use = "p"))
  hub <- mod_genes[GS > 0.2 & MM > 0.8]
  wgcna_hub <- c(wgcna_hub, hub)
}
wgcna_hub <- unique(wgcna_hub)
cat("WGCNA hub genes:", length(wgcna_hub), "\n")


## ############################################################
## STEP 3 — MACHINE LEARNING TARGET SCREENING
## ############################################################
## Input: DEG ∩ WGCNA hub genes
## LASSO (10-fold CV) + Random Forest (Gini Top 20%)
## Final key targets = LASSO ∩ RF
## ############################################################

candidate_genes <- intersect(sig_degs$gene, wgcna_hub)
cat("Candidate genes (DEG ∩ WGCNA hub):", length(candidate_genes), "\n")

## Prepare feature matrix
X <- t(vsd_mat[candidate_genes, , drop = FALSE])
y <- as.numeric(meta$condition == "Seps_P")

## ---- LASSO (10-fold CV) ----
set.seed(42)
cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial",
                      nfolds = 10, type.measure = "auc")

lasso_coef  <- coef(cv_lasso, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]
cat("LASSO selected:", length(lasso_genes), "genes\n")

## ---- Random Forest ----
set.seed(42)
rf_model <- randomForest(x = X, y = factor(y),
                         ntree = 500, importance = TRUE)
rf_imp   <- importance(rf_model, type = 2)   # MeanDecreaseGini
rf_top   <- rownames(rf_imp)[rf_imp[, 1] >= quantile(rf_imp[, 1], 0.80)]
cat("RF top 20%:", length(rf_top), "genes\n")

## ---- Intersection → key targets ----
key_targets <- intersect(lasso_genes, rf_top)
cat("Key targets (LASSO ∩ RF):", length(key_targets), "\n")
cat(" ", paste(key_targets, collapse = ", "), "\n")

## ---- Fig 5C: Venn — LASSO ∩ RF ----
venn_ml <- getvenndata(lasso_genes, rf_top,
                       set_names = c("LASSO", "Random Forest"))

p5c <- ggvenn_plot(venn_ml,
                   set.color = pal_lancet[c(1, 4)],
                   stroke.size = 0.5,
                   text.size = 4) +
  ggtitle("Machine learning target screening") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

p5c
ggsave("output/figure5/fig5c_ml_venn.pdf", p5c,
       width = 120, height = 110, units = "mm", dpi = 600,
       device = cairo_pdf)


## ############################################################
## STEP 4 — TCMDATA INTEGRATION
## ############################################################
## key_targets from ML screening → herb_enricher → Venn
## → Sankey → PPI → enrichment → ppi_knock
## ############################################################

## Expand to include all DEGs that overlap with WGCNA hubs
## for a more informative herb enrichment and PPI
disease_targets <- candidate_genes   # DEG ∩ WGCNA hub gene set

## ---- Fig 5D: Herb ORA ----
herb_res <- herb_enricher(genes = disease_targets, type = "Herb_pinyin_name")

p5d <- dotplot(herb_res, showCategory = 15) +
  labs(x = "Gene Ratio", y = NULL,
       title = "Herb enrichment (sepsis targets)") +
  theme_sci(base_size = 9) +
  theme(axis.text.y = element_text(hjust = 1))

p5d
ggsave("output/figure5/fig5d_herb_ORA.pdf", p5d,
       width = 160, height = 140, units = "mm", dpi = 600,
       device = cairo_pdf)

## Select top 2 enriched herbs
herb_df    <- as.data.frame(herb_res)
top_herbs  <- head(herb_df$Description, 2)
cat("Top enriched herbs:", paste(top_herbs, collapse = ", "), "\n")

herb1_data <- search_herb(top_herbs[1], type = "Herb_pinyin_name")
herb2_data <- search_herb(top_herbs[2], type = "Herb_pinyin_name")
herb1_targets <- unique(herb1_data$target)
herb2_targets <- unique(herb2_data$target)

## ---- Fig 5E: Three-set Venn (disease ∩ herb targets) ----
venn_df <- getvenndata(disease_targets, herb1_targets, herb2_targets,
                       set_names = c("DEG-Hub", top_herbs[1], top_herbs[2]))

p5e <- ggvenn_plot(venn_df,
                   set.color = pal_lancet[c(4, 1, 2)],
                   stroke.size = 0.5,
                   text.size = 3.5) +
  ggtitle("Target intersection") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

p5e
ggsave("output/figure5/fig5e_venn.pdf", p5e,
       width = 140, height = 120, units = "mm", dpi = 600,
       device = cairo_pdf)

## Extract intersection genes
venn_result <- getvennresult(venn_df)
all_herb_targets <- unique(c(herb1_targets, herb2_targets))
intersect_genes  <- intersect(disease_targets, all_herb_targets)
cat("Intersection genes:", length(intersect_genes), "\n")


## ---- Fig 5F: Sankey (herb–compound–target) ----
sankey_data <- rbind(herb1_data, herb2_data)
sankey_data <- sankey_data[sankey_data$target %in% intersect_genes, ]

top_compounds <- names(sort(table(sankey_data$molecule), decreasing = TRUE))
top_compounds <- head(top_compounds, 10)
sankey_subset <- sankey_data[sankey_data$molecule %in% top_compounds, ]
top_tgts <- names(sort(table(sankey_subset$target), decreasing = TRUE))
top_tgts <- head(top_tgts, 20)
sankey_subset <- sankey_subset[sankey_subset$target %in% top_tgts, ]

p5f <- tcm_sankey(sankey_subset, font_size = 2.8, alpha = 0.35) +
  ggtitle("Herb\u2013Compound\u2013Target network") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

p5f
ggsave("output/figure5/fig5f_sankey.pdf", p5f,
       width = 200, height = 160, units = "mm", dpi = 600,
       device = cairo_pdf)


## ---- PPI construction ----
ppi_raw <- getPPI(intersect_genes, taxID = 9606)

## Filter and score
ppi_filtered <- ppi_subset(ppi_raw, score_cutoff = 0.7)
low_deg <- V(ppi_filtered)[degree(ppi_filtered) < 2]$name
if (length(low_deg) > 0) {
  ppi_filtered <- delete_vertices(ppi_filtered, low_deg)
}
cat("PPI after filtering:", vcount(ppi_filtered), "nodes,",
    ecount(ppi_filtered), "edges\n")

## Compute 19 topological metrics + composite hub ranking
ppi_scored <- compute_nodeinfo(ppi_filtered, weight_attr = "score")
rank_res   <- rank_ppi_nodes(ppi_scored, use_weight = TRUE)
ppi_ranked <- rank_res$graph
rank_df    <- rank_res$table

V(ppi_ranked)$btw_raw <- igraph::betweenness(ppi_ranked, directed = FALSE,
                                              normalized = FALSE)

cat("Top 10 hub genes:", paste(rank_df$name[1:min(10, nrow(rank_df))],
                               collapse = ", "), "\n")


## ---- Fig 5G: PPI network with composite hub score ----
set.seed(42)
p5g <- ggplot(ppi_ranked, layout = "fr") +
  geom_edge(alpha = 0.2, color = pal_edge, linewidth = 0.45) +
  geom_point(aes(size = btw_raw, color = Score_network), alpha = 0.9) +
  scale_color_gradientn(
    colours = c("#3C5488", "#4DBBD5", "#F7F7F7", "#F39B7F", "#E64B35"),
    name = "Hub score"
  ) +
  scale_size_continuous(range = c(1.5, 7), name = "Betweenness") +
  geom_text_repel(
    aes(label = ifelse(Score_network >= quantile(Score_network, 0.75),
                       name, "")),
    size = 2.2, fontface = "bold.italic",
    max.overlaps = 25, segment.alpha = 0.4,
    segment.size = 0.25, box.padding = 0.35, point.padding = 0.2
  ) +
  ggtitle("PPI network with composite hub score") +
  theme_void() +
  theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.title = element_text(face = "bold", size = 9),
        legend.text  = element_text(size = 8))

p5g
ggsave("output/figure5/fig5g_ppi_hub.pdf", p5g,
       width = 180, height = 160, units = "mm", dpi = 600,
       device = cairo_pdf)


## ---- Fig 5H: Network knockout ----
top1_hub <- rank_df$name[1]
cat("Knockout target:", top1_hub, "\n")

ko_res <- ppi_knock(ppi_scored, targets = top1_hub,
                    n_perm = 100, weight_attr = "score", seed = 42)

cat("Knockout Summary:\n")
print(ko_res$Summary)
cat("Z_Total =", ko_res$Total_Score, "  P =", ko_res$Total_Pvalue, "\n")

## Build residual graph for visualization
ppi_ko <- delete_vertices(ppi_scored, top1_hub)
iso <- V(ppi_ko)[degree(ppi_ko) == 0]$name
if (length(iso) > 0) ppi_ko <- delete_vertices(ppi_ko, iso)
V(ppi_ko)$btw_raw <- igraph::betweenness(ppi_ko, directed = FALSE,
                                          normalized = FALSE)

set.seed(42)
p5h <- ggplot(ppi_ko, layout = "fr") +
  geom_edge(alpha = 0.2, color = pal_edge, linewidth = 0.45) +
  geom_point(aes(size = btw_raw, color = degree), alpha = 0.9) +
  scale_color_distiller(palette = "Spectral", direction = -1, name = "Degree") +
  scale_size_continuous(range = c(1.5, 7), name = "Betweenness") +
  geom_text_repel(
    aes(label = name), size = 2.0, fontface = "bold.italic",
    max.overlaps = 30, segment.alpha = 0.4,
    segment.size = 0.25, box.padding = 0.3, point.padding = 0.2
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0(top1_hub, " knocked out"),
           fontface = "bold.italic", size = 3.5, color = "#E64B35") +
  ggtitle(paste0("Residual PPI after ", top1_hub, " removal")) +
  theme_void() +
  theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.title = element_text(face = "bold", size = 9),
        legend.text  = element_text(size = 8))

p5h
ggsave("output/figure5/fig5h_knockout.pdf", p5h,
       width = 180, height = 160, units = "mm", dpi = 600,
       device = cairo_pdf)

cat("\n=== Case study complete ===\n")
cat("All panels saved to output/figure5/\n")
