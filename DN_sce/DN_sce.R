#' single-cell analysis for GSE209781
library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(ggsci)
library(ggsc)

## preprocess
# load data and subset
dataloader <- function(dir){
  dkd <- readRDS(dir)
  dkd[["percent.mt"]] <- PercentageFeatureSet(dkd, pattern = "^MT-")
  dkd <- subset(dkd, subset = nFeature_RNA > 200 & nFeature_RNA < 6100 & percent.mt < 34)
  print(table(dkd$group))
  return(dkd)
}

dkd <- dataloader(dir = "../data/GSE209781_srt.rds")

# normalization
normalization <- function(srt){
  srt <- NormalizeData(srt) |> 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
    ScaleData()
  return(srt)
}

dkd <- normalization(dkd)

# PCA
dkd <- RunPCA(dkd)
ElbowPlot(dkd)

# clustering
dkd <- FindNeighbors(dkd, reduction = "pca", dims = 1:8)
dkd <- FindClusters(dkd, resolution = 0.38) 
dkd <- RunUMAP(dkd, reduction = "pca", dims = 1:8)
dkd <- RunTSNE(dkd, reduction = "pca", dims = 1:8)

p1 <- sc_dim(dkd, reduction = "tsne", geom = geom_bgpoint, pointsize = 1)
p1



## Findmarkers
markers <- FindAllMarkers(dkd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
d <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC)

print(d, n = 240) 

## cell annotation
# marker
new.cluster.ids <- c(
  "T_NK",             # Cluster 0
  "Proximal tubule",  # Cluster 1
  "Proximal tubule",  # Cluster 2
  "DCT-LOH",          # Cluster 3 (SLC12A1)
  "DCT",              # Cluster 4 (KRT19/MMP7)
  "EC",               # Cluster 5
  "LOH",              # Cluster 6 (CLDN14)
  "Fibroblast",       # Cluster 7 (RGS5/TAGLN)
  "B cells",          # Cluster 8
  "Mono",             # Cluster 9 (CD1C/LYZ)
  "EC",               # Cluster 10
  "Neutrophils",      # Cluster 11
  "CD-PC",            # Cluster 12
  "CD-IC",            # Cluster 13
  "T cells",             # Cluster 14
  "Fibroblast"        # Cluster 15
)


names(new.cluster.ids) <- levels(dkd)
dkd <- RenameIdents(dkd, new.cluster.ids)
dkd$cell_type <- Idents(dkd)

## tsne after annotation
p2 <- sc_dim(dkd, reduction = "tsne", geom = geom_bgpoint, pointsize = 1) + theme_classic()
p2


saveRDS(dkd, file = "../data/dkd_sce_final.rds")
