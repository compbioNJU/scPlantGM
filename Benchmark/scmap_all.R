# Split by dataset-level cross-validation

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly= TRUE)
print(args)

library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(dplyr)

source("process.R")

species = "ara"
organ = args[1]
k = args[2]

alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

ref_list = get_sample(info, organ, k, sheet="ara_dataset")[[1]]
que_list = get_sample(info, organ, k, sheet="ara_dataset")[[2]]

ref_data = alldata[,which(alldata@meta.data$orig.ident %in% ref_list)]
que_data = alldata[,which(alldata@meta.data$orig.ident %in% que_list)]

#record time
startTime <- Sys.time()

#build reference
ref_matrix = ref_data@assays$RNA@counts
ref_ann = as.data.frame(ref_data@meta.data['annotation'])

ref_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(ref_matrix)), colData = ref_ann)
logcounts(ref_sce) <- log2(normcounts(ref_sce) + 1)
rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
#isSpike(ref_sce, "ERCC") <- grepl("^ERCC-", rownames(ref_sce))
ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]

#find feature gene
ref_sce <- selectFeatures(ref_sce, n_features = 2000, suppress_plot = FALSE)
#check feature gene
#table(rowData(ref_sce)$scmap_features)
ref_sce <- indexCluster(ref_sce, cluster_col='annotation')

#build query
que_matrix = que_data@assays$RNA@counts
que_sce = SingleCellExperiment(assays = list(normcounts = as.matrix(que_matrix)))
logcounts(que_sce) <- log2(normcounts(que_sce) + 1)
rowData(que_sce)$feature_symbol <- rownames(que_sce)
#isSpike(que_sce, "ERCC") <- grepl("^ERCC-", rownames(que_sce))
que_sce <- que_sce[!duplicated(rownames(que_sce)), ]

scmapCluster_results <- scmapCluster(
  projection = que_sce, 
  index_list = list(
  metadata(ref_sce)$scmap_cluster_index
  ),
  threshold = 0.1
)

scmapCluster_results$Cell <- colnames(que_matrix)
scmapCluster_results$annotation <- que_data@meta.data$anno
scmapCluster_results$prob <- get_scmap_probmatrix(que_sce, list(metadata(ref_sce)$scmap_cluster_index),threshold = 0.1)
saveRDS(scmapCluster_results, paste('scmap_result/scmap_',organ,'_K',k,'_model.rds',sep = ''))

endTime <- Sys.time()
print(startTime)
print(endTime)
runtime = endTime - startTime
print(sprintf("%s time %s: %s", organ, k, endTime - startTime))


##############################################################################################
# Split by sample-level cross-validation

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly= TRUE)
print(args)

library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(dplyr)
library(readxl)

source("code/process.R")

species <- "ara"
organ = args[1] # "Leaf", "SAA", "Inflorescence",  "Root", "Pollen"
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

  ref_data = alldata[,which(alldata@meta.data$orig.ident %in% ref_list)]
  que_data = alldata[,which(alldata@meta.data$orig.ident %in% que_list)]

  #record time
  startTime <- Sys.time()

  #build reference
  ref_matrix = ref_data@assays$RNA@counts
  ref_ann = as.data.frame(ref_data@meta.data['annotation'])

  ref_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(ref_matrix)), colData = ref_ann)
  logcounts(ref_sce) <- log2(normcounts(ref_sce) + 1)
  rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
  #isSpike(ref_sce, "ERCC") <- grepl("^ERCC-", rownames(ref_sce))
  ref_sce <- ref_sce[!duplicated(rownames(ref_sce)), ]

  #find feature gene
  ref_sce <- selectFeatures(ref_sce, n_features = 2000, suppress_plot = FALSE)
  #check feature gene
  #table(rowData(ref_sce)$scmap_features)
  ref_sce <- indexCluster(ref_sce, cluster_col='annotation')

  #build query
  que_matrix = que_data@assays$RNA@counts
  que_sce = SingleCellExperiment(assays = list(normcounts = as.matrix(que_matrix)))
  logcounts(que_sce) <- log2(normcounts(que_sce) + 1)
  rowData(que_sce)$feature_symbol <- rownames(que_sce)
  #isSpike(que_sce, "ERCC") <- grepl("^ERCC-", rownames(que_sce))
  que_sce <- que_sce[!duplicated(rownames(que_sce)), ]

  scmapCluster_results <- scmapCluster(
    projection = que_sce, 
    index_list = list(
    metadata(ref_sce)$scmap_cluster_index
    ),
    threshold = 0.1
  )

  scmapCluster_results$Cell <- colnames(que_matrix)
  scmapCluster_results$annotation <- que_data@meta.data$anno
  scmapCluster_results$prob <- get_scmap_probmatrix(que_sce, list(metadata(ref_sce)$scmap_cluster_index),threshold = 0.1)
  saveRDS(scmapCluster_results, paste('scmap_result/scmap_sample_', organ, "_", dataset, '_K',k,'_model.rds',sep = ''))

  endTime <- Sys.time()
  print(startTime)
  print(endTime)
  runtime = endTime - startTime
  print(sprintf("%s time %s: %s", organ, k, endTime - startTime))

}
