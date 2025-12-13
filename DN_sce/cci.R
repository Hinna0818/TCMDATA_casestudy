#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(Seurat)
library(ggplot2)
library(aplot)

## load data
dkd <- readRDS("../data/dkd.rds")
data.input <- GetAssayData(dkd, assay = "RNA", layer = "data")
labels <- Idents(dkd)
meta <- data.frame(group = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 50 * 1024^3)

## pipeline
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "../data/dkd_cci.rds")

## visu
groupSize <- as.numeric(table(cellchat@idents))
p1 <- plot_list(
  ~CellChat::netVisual_circle(cellchat@net$count, 
                              vertex.weight = groupSize, weight.scale = T, 
                              label.edge= F, title.name = "Number of interactions"),
  ~CellChat::netVisual_circle(cellchat@net$weight, 
                              vertex.weight = groupSize, weight.scale = T, 
                              label.edge= F, title.name = "Interaction weights/strength"))

p1


sources.use <- c("Macrophages", "Neutrophils", "NK/T Cells")
targets.use <- c("Fibroblasts/Pericytes", "Injured Tubule", "Endothelial Cells")
p2 <- CellChat::netVisual_bubble(cellchat, 
                 sources.use = sources.use, 
                 targets.use = targets.use, 
                 remove.isolate = FALSE)


p2

p3 <- CellChat::netVisual_aggregate(cellchat, 
                              signaling = c("CXCL"),  
                              layout = "circle")      

p3
