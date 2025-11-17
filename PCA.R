#' PCA plot for GSE142025
rm(list = ls())
options(stringsAsFactors = F)
library(FactoMineR)
library(factoextra)

dir.create("GSE142025-QC")

## PCA
load(file ='./GSE142025/step1-output.Rdata')

cg <- names(tail(sort(apply(dat, 1, sd)), 1000))
head(cg)

dat <- t(dat[cg, ])
dat <- as.data.frame(dat)
dat.pca <- PCA(dat , graph = FALSE)
fviz_pca_ind(dat.pca, geom.ind ="point", col.ind = group_list, addEllipses = T,
             legend.title ="Groups")
