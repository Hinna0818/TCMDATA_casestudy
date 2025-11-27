#' single-cell analysis for GSE209781
library(Seurat)
library(ggplot2)
library(dplyr)

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

library(ggsci)
p1 <- DimPlot(dkd, 
                    reduction = "umap", 
                    group.by = "seurat_clusters",
                    label = T,          
                    label.size = 5,         
                    label.box = TRUE,       
                    repel = TRUE,           
                    pt.size = 0.05) +       
  scale_color_npg() +                      
  theme(legend.position = "none") +       
  ggtitle("Cell Clusters (UMAP)") +         
  NoAxes()                                

p1


## Findmarkers
markers <- FindAllMarkers(dkd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
d <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(d, n=100) 

## cell annotation
new_cluster_ids <- c(
  "NK/T Cells",             # 0
  "Proximal Tubule",        # 1
  "Proximal Tubule",        # 2
  "LOH/Distal Tubule",      # 3
  "Injured Tubule",         # 4
  "Endothelial Cells",      # 5
  "LOH (Thin Limb)",        # 6
  "Mesangial/Fibroblasts",  # 7
  "B Cells",                # 8
  "Macrophages",            # 9
  "Endothelial Cells",      # 10
  "Neutrophils",            # 11
  "Collecting Duct-PC",     # 12
  "Collecting Duct-IC",     # 13
  "NK/T Cells",             # 14
  "Fibroblasts/Pericytes"   # 15
)

Idents(dkd) <- "seurat_clusters"

names(new_cluster_ids) <- levels(dkd)
dkd <- RenameIdents(dkd, new_cluster_ids)

dkd$cell_type <- Idents(dkd)

DimPlot(dkd, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("Annotated Cell Types") +
  NoAxes()

p2 <- DimPlot(dkd, 
              reduction = "umap", 
              group.by = "cell_type",
              label = T,          
              label.size = 4,         
              label.box = TRUE,       
              repel = TRUE,           
              pt.size = 0.05) +       
  scale_color_npg() +                      
  theme(legend.position = "none") +       
  ggtitle("Cell Clusters (UMAP)") +         
  NoAxes()   
p2


