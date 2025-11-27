## load GSE209781
library(Seurat)
library(dplyr)

base_dir <- "GSE209781_RAW"

all_dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
tenx_dirs <- all_dirs[
  grepl("NM0|DKD0", basename(all_dirs)) &  
    !grepl("\\.tar.gz$", all_dirs)
]

tenx_dirs <- list.dirs("GSE209781_RAW", recursive = TRUE)

is_10x <- file.exists(file.path(tenx_dirs, "barcodes.tsv.gz")) &
          file.exists(file.path(tenx_dirs, "features.tsv.gz")) &
          file.exists(file.path(tenx_dirs, "matrix.mtx.gz"))

tenx_dirs <- tenx_dirs[is_10x]
tenx_dirs

sample_ids <- basename(tenx_dirs)  

sce_list <- lapply(seq_along(tenx_dirs), function(i) {
  message("Reading: ", tenx_dirs[i])
  
  counts <- Read10X(data.dir = tenx_dirs[i])
  
  sce <- CreateSeuratObject(
    counts = counts,
    project = sample_ids[i],
    min.cells = 5,
    min.features = 300
  )
  
  sce$sample_id <- sample_ids[i]
  sce$group <- ifelse(grepl("^NM", sample_ids[i]), "Control", "DKD")
  
  sce
})

dkd_all <- Reduce(function(x, y) {
  merge(x, y, add.cell.ids = c(x@meta.data$sample_id[1], y@meta.data$sample_id[1]))
}, sce_list)

dkd_all <- JoinLayers(dkd_all)
# saveRDS(dkd_all, "GSE209781_srt.rds")

## preprocess
dkd_all[["percent.mt"]] <- PercentageFeatureSet(dkd_all, pattern = "^MT-")
VlnPlot(dkd_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
