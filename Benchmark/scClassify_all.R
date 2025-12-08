# Split by dataset-level cross-validation

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly= TRUE)
print(args)

library(scClassify)
library(Seurat)
library(scater)
library(ggplot2)
library(gridExtra)
library(dplyr)

source("process.R")

species <- "ara"
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis
alldata = readRDS('Arabidopsis_clean.rds')

organ = args[1]
k = args[2]
k <- as.numeric(k)

ref_list = get_sample(info, organ, k, sheet="ara_dataset")[[1]]
que_list = get_sample(info, organ, k, sheet="ara_dataset")[[2]]

ref_cds = alldata[,which(alldata@meta.data$orig.ident %in% ref_list)]
que_cds = alldata[,which(alldata@meta.data$orig.ident %in% que_list)]


ref_cds <- ref_cds %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
ref_mat = as(ref_cds@assays$SCT@counts,"dgCMatrix")

que_cds = que_cds %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
que_mat = as(que_cds@assays$SCT@counts,"dgCMatrix")

ref_anno = ref_cds@meta.data$annotation
que_anno = que_cds@meta.data$annotation

#record time
startTime <- Sys.time()

trainRes = train_scClassify(ref_mat, ref_anno, tree = "HOPACH",
                            selectFeatures = "limma",
                            topN = 50,
                            hopach_kmax = 5,
                            pSig = 0.05,
                            cellType_tree = NULL,
                            weightsCal = FALSE,
                            parallel= FALSE,
                            BPPARAM = BiocParallel::SerialParam(),
                            verbose= TRUE,
                            returnList = TRUE)

scClassify_res <- predict_scClassify(exprsMat_test = que_mat,
                            trainRes,
                            cellTypes_test = que_anno,
                            algorithm = "WKNN",
                            similarity = c("pearson"),
                            verbose = FALSE)

#plotCellTypeTree(scClassify_res$trainRes@cellTypeTree)
#ggsave(paste('scClassify_result/scClassify_',organ,'_K',k,'.png',sep = ''),bg = 'white')
saveRDS(scClassify_res, paste('scClassify_result/scClassify_',organ,'_K',k,'.rds',sep = ''))
saveRDS(trainRes, paste('scClassify_result/scClassify_',organ,'_K',k,'_model.rds',sep = ''))

endTime <- Sys.time()
runtime = endTime - startTime
print(organ)
print(k)
print(startTime)
print(endTime)
print(paste("scClassify time",endTime - startTime, sep=":"))

###########################################################################################################
# Split by sample-level cross-validation

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly= TRUE)
print(args)

library(scClassify)
library(Seurat)
library(scater)
library(gridExtra)
library(dplyr)
library(readxl)

source("code/process.R")

species <- "ara"
organ = args[1]
n = args[2]


alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis
data <- read_excel(path="ref_query.xlsx", sheet="ara_sample")

alldataset <- unique(data$Dataset[which(data$Organ==organ)])
dataset <- alldataset[as.numeric(n)]
ks <- max(data$K[which(data$Dataset==dataset)])


for(k in 1:ks){
  
  ref_list = get_sample(info, organ, kvalue=k, dataset=dataset, sheet="ara_sample")[[1]]
  que_list = get_sample(info, organ, kvalue=k, dataset=dataset, sheet="ara_sample")[[2]]

  
  ref_cds = alldata[,which(alldata@meta.data$orig.ident %in% ref_list)]
  que_cds = alldata[,which(alldata@meta.data$orig.ident %in% que_list)]


  ref_cds <- ref_cds %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()
  ref_mat = as(ref_cds@assays$SCT@counts,"dgCMatrix")

  que_cds = que_cds %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()
  que_mat = as(que_cds@assays$SCT@counts,"dgCMatrix")

  ref_anno = ref_cds@meta.data$annotation
  que_anno = que_cds@meta.data$annotation

  #record time
  startTime <- Sys.time()

  trainRes = train_scClassify(ref_mat, ref_anno, tree = "HOPACH",
                              selectFeatures = "limma",
                              topN = 50,
                              hopach_kmax = 5,
                              pSig = 0.05,
                              cellType_tree = NULL,
                              weightsCal = FALSE,
                              parallel= FALSE,
                              BPPARAM = BiocParallel::SerialParam(),
                              verbose= TRUE,
                              returnList = TRUE)

  scClassify_res <- predict_scClassify(exprsMat_test = que_mat,
                              trainRes,
                              cellTypes_test = que_anno,
                              algorithm = "WKNN",
                              similarity = c("pearson"),
                              verbose = FALSE)

  scClassify_res$annotation <- que_cds@meta.data$annotation

  #plotCellTypeTree(scClassify_res$trainRes@cellTypeTree)
  #ggsave(paste('scClassify_result/scClassify_',organ,'_K',k,'.png',sep = ''),bg = 'white')
  saveRDS(scClassify_res, paste('scClassify_result/scClassify_sample_', organ, "_", dataset,'_K',k,'.rds',sep = ''))
  saveRDS(trainRes, paste('scClassify_result/scClassify_sample_', organ, "_", dataset, '_K',k,'_model.rds',sep = ''))

  endTime <- Sys.time()
  runtime = endTime - startTime
  print(organ)
  print(k)
  print(startTime)
  print(endTime)
  print(paste("scClassify time", endTime - startTime, sep=":"))

}
