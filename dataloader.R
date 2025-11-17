#' load RNA-seq data from GSE142025
rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 20)
library(data.table)
library(tinyarray)
library(AnnoProbe)
library(org.Hs.eg.db)
library(GEOquery)
library(edgeR)

setwd("E:/tcmdata_casestudy/data")
gse <-"GSE142025"
dir.create(gse)

if(T){
  # load counts table from GEO
  urld <-"https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&"
  gse <-"GSE142025"
  path <- paste0(urld,"acc=", gse,"&format=file&file=",gse,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
  file <- paste0(gse,"_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  download.file(url = path, destfile = file)
  
  tbl <- as.matrix(data.table::fread(file, header = T, colClasses = "integer"), rownames = 1)
  ensembl_matrix <- as.data.frame(tbl)

  e2s <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = rownames(ensembl_matrix),
                             columns ="SYMBOL",
                             keytype ="ENTREZID")
  ids <- na.omit(e2s)
  ids <- ids[!duplicated(ids$SYMBOL), ]
  ids <- ids[!duplicated(ids$ENTREZID), ]
  symbol_matrix <- ensembl_matrix[match(ids$ENTREZID, rownames(ensembl_matrix)), ]
  rownames(symbol_matrix) <- ids$SYMBOL
}

# get sample info
if(T){
  gset <- getGEO(gse, destdir ='.', getGPL = F)
  pd <- pData(gset[[1]])
  com <- intersect(rownames(pd), colnames(symbol_matrix))
  symbol_matrix <- symbol_matrix[, com]
  pd <- pd[com,]
  group_list <- pd$`group:ch1`
  table(group_list)
}

# save
dat <- log10(edgeR::cpm(symbol_matrix)+1)
save(symbol_matrix, dat, group_list, file ='./GSE142025/step1-output.Rdata')
