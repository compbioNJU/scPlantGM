setwd("plants")

# Save downsampled data
source("process.R")

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis
data("reference_Arabidopsis", package = "scPlantGM")
jaccard_mat <- reference_Arabidopsis$sim_mat

# Downsampling by cluster
organ = "Leaf" # "Leaf" 和 "Root"
for(percent in c(0.05, 0.1, 0.2, 0.5)){ # downsampling 5,10,20,50%
    for(k in 1:5){

        ref_list = get_sample(info, organ, k)[[1]]
        que_list = get_sample(info, organ, k)[[2]]

        cluster1 <- info %>% filter(Sample %in% ref_list) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
        cluster2 <- info %>% filter(Sample %in% que_list) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
        cluster2 <- intersect(cluster2,colnames(jaccard_mat))

        set.seed(123)
        downsample1 <- sample(1:length(cluster1), round(percent*(length(cluster1))))

        downsample <- list()
        downsample$ref <- cluster1[downsample1]
        saveRDS(downsample$ref, file = sprintf("downsample_ids_%s_%s_%s_ref.rds", organ, percent, k))
        downsample$query <- cluster2
        saveRDS(downsample$query, file = sprintf("downsample_ids_%s_%s_%s_query.rds", organ, percent, k))
        saveRDS(downsample, file = sprintf("downsample_ids_%s_%s_%s.rds", organ, percent, k))
    }
    print(percent)
}  



# Downsampling by dataset
organ = "Leaf" # "Leaf" 和 "Root"
for(percent in c(0.1, 0.2, 0.5)){ # downsampling 10,20,50%
    for(k in 1:5){

        ref_list = get_sample(info, organ, k)[[1]]
        que_list = get_sample(info, organ, k)[[2]]

        dataset1 <- info %>% filter(Sample %in% ref_list) %>% select(Dataset) %>% unique() %>% unlist() %>% as.character()
        dataset2 <- info %>% filter(Sample %in% que_list) %>% select(Dataset) %>% unique() %>% unlist() %>% as.character()

        set.seed(123)
        downsample1 <- sample(1:length(dataset1), max(round(percent*(length(dataset1))),1))
        # print(downsample)

        downsample <- list()
        downsample$ref <- unique(info$Sample[info$Dataset %in% dataset1[downsample1]])
        saveRDS(downsample$ref, file = sprintf("dataset_downsample_ids_%s_%s_%s_ref.rds", organ, percent, k))
        downsample$query <- unique(info$Sample[info$Dataset %in% dataset2])
        saveRDS(downsample$query, file = sprintf("dataset_downsample_ids_%s_%s_%s_query.rds", organ, percent, k))

        saveRDS(downsample, file = sprintf("dataset_downsample_ids_%s_%s_%s.rds", organ, percent, k))
    }
    print(percent)
}




####################################################################################################################################################
# Split by dataset-level cross-validation
# scPlantGM
library(dplyr)
library(openxlsx)
library(scPlantGM)

source("process.R")

species <- "ara"
organ <- "Root"
res <- 1
topn <- 100
pthres <- 0.8
xthres <- 0.2

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

data("layer_info_Arabidopsis", package = "scPlantGM")
# load("data/layer_info_Arabidopsis.rda")
layer_info <- layer_info_Arabidopsis %>% 
                dplyr::filter(Organ==organ)

data("reference_Arabidopsis", package = "scPlantGM")
gene_list <- reference_Arabidopsis$gene_list
jaccard_mat <- reference_Arabidopsis$sim_mat

data("acc_list_sta_all_Arabidopsis", package = "scPlantGM")
acc_list_sta <- acc_list_sta_all_Arabidopsis$layer0
acc_list_sta1 <- acc_list_sta_all_Arabidopsis$layer1
acc_list_sta2 <- acc_list_sta_all_Arabidopsis$layer2
acc_list_sta3 <- acc_list_sta_all_Arabidopsis$layer3

for(percent in c(0.05, 0.1, 0.2, 0.5)){

  print(organ)

  for(kvalue in 1:5){

    downsample <- readRDS(sprintf("downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refcluster <- downsample$ref
    querycluster <- downsample$query

    #按照module预测
    result <- predict_based_module_weight(info, jaccard_mat, layer=3,
                                acc_list_sta, acc_list_sta1, acc_list_sta2, acc_list_sta3,
                                refsample = refcluster, querysample = querycluster, 
                                gene_list, pthres, x_thres = xthres,
                                issample = FALSE, iscluster = TRUE)
    saveRDS(result, sprintf("scPlantGM_%s_%s_%s.rds", organ, percent, kvalue))
    }
}





####################################################################################
# singleR
library(SingleR)
library(openxlsx)

source("process.R")

alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

organ <- "Leaf"
percent = 0.05 
for(kvalue in 1:5){
    # kvalue = 5

    downsample <- readRDS(sprintf("downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refcluster <- downsample$ref
    querycluster <- downsample$query

    cells1 <- info$Cell[which(info$Cluster %in% refcluster)]
    cells2 <- info$Cell[which(info$Cluster %in% querycluster)]

    train_rds = subset(alldata, cells = cells1)
    test_rds = subset(alldata, cells = cells2)

    train_sce = sceasy::convertFormat(train_rds, from="seurat", to="sce")
    test_sce = sceasy::convertFormat(test_rds, from="seurat", to="sce")

    train_ann = train_rds$annotation
    test_ann = test_rds$annotation

    train_sce$celltype = train_ann
    test_sce$celltype = test_ann

    #train_result <- SingleR(train_sce, train_sce, labels = train_sce$celltype)
    test_result = SingleR(test_sce, train_sce, labels = train_sce$celltype)

    test_result$cell <- info$Cell[which(info$Cluster %in% querycluster)]
    test_result$annotation <- info$annotation[which(info$Cluster %in% querycluster)]

    saveRDS(test_result, paste('SingleR_', organ, "_", percent , '_K', kvalue,'.rds',sep = ''))
}



#########################################################################################################################
# scmap
library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(dplyr)
library(readxl)
library(openxlsx)

source("process.R")

organ <- "Leaf"
alldata = readRDS('methods/Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

percent = 0.05 
for(kvalue in 1:5){
    # kvalue <- 5

    downsample <- readRDS(sprintf("downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refcluster <- downsample$ref
    querycluster <- downsample$query

    cells1 <- info$Cell[which(info$Cluster %in% refcluster)]
    cells2 <- info$Cell[which(info$Cluster %in% querycluster)]

    ref_data = subset(alldata, cells = cells1)
    que_data = subset(alldata, cells = cells2)
  
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
    saveRDS(scmapCluster_results, paste('scmap_', organ, "_", percent, '_K', kvalue, '.rds', sep = ''))
}



########################################################################################################################
# scClassify
library(scClassify)
library(Seurat)
library(scater)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(readxl)
library(openxlsx)

source("process.R")

organ <- "Leaf"
alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

percent = 0.05 
for(kvalue in 1:5){
    # kvalue = 2

    downsample <- readRDS(sprintf("downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refcluster <- downsample$ref
    querycluster <- downsample$query

    cells1 <- info$Cell[which(info$Cluster %in% refcluster)]
    cells2 <- info$Cell[which(info$Cluster %in% querycluster)]

    ref_cds = subset(alldata, cells = cells1)
    que_cds = subset(alldata, cells = cells2)

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

    ref_anno = ref_cds$annotation
    que_anno = que_cds$annotation

    # train model
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

    scClassify_res$annotation <- que_anno

    saveRDS(scClassify_res, paste('scClassify_', organ, "_", percent, '_K', kvalue,'.rds', sep = ''))
    saveRDS(trainRes, paste('scClassify_', organ, "_", percent, '_K', kvalue, '_model.rds', sep = ''))

    print(kvalue)
}



####################################################################################################################
####################################################################################################################
# Split by sample-level cross-validation
# scPlantGM

source("process.R")

library(openxlsx)
library(scPlantGM)

species <- "ara"
res <- 1
topn <- 100
pthres <- 0.8
xthres <- 0.2

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

data("layer_info_Arabidopsis", package = "scPlantGM")
# load("data/layer_info_Arabidopsis.rda")
layer_info <- layer_info_Arabidopsis %>% 
                dplyr::filter(Organ==organ)

data("reference_Arabidopsis", package = "scPlantGM")
gene_list <- reference_Arabidopsis$gene_list
jaccard_mat <- reference_Arabidopsis$sim_mat

data("acc_list_sta_all_Arabidopsis", package = "scPlantGM")
acc_list_sta <- acc_list_sta_all_Arabidopsis$layer0
acc_list_sta1 <- acc_list_sta_all_Arabidopsis$layer1
acc_list_sta2 <- acc_list_sta_all_Arabidopsis$layer2
acc_list_sta3 <- acc_list_sta_all_Arabidopsis$layer3

organ <- "Root"
for(percent in c(0.1, 0.2, 0.5)){

  for(kvalue in 1:5){

    downsample <- readRDS(sprintf("dataset_downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refsample <- downsample$ref
    querysample <- downsample$query

    #按照module预测
    result <- predict_based_module_weight(info, jaccard_mat, layer=3,
                                acc_list_sta, acc_list_sta1, acc_list_sta2, acc_list_sta3,
                                refsample = refsample, querysample = querysample, 
                                gene_list, pthres, x_thres = xthres)
    saveRDS(result, sprintf("dataset_scPlantGM_%s_%s_%s.rds", organ, percent, kvalue))

    print(kvalue)
    }
    print(percent)
}


##################################################################################################################################
# singleR
library(SingleR)
library(Seurat)
library(magrittr)
library(readxl)
library(openxlsx)

source("code/process.R")

alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis


organ <- "Leaf"
percent = 0.05 
for(kvalue in 1:5){
    # kvalue = 5

    downsample <- readRDS(sprintf("dataset_downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refsample <- downsample$ref
    querysample <- downsample$query

    cells1 <- info$Cell[which(info$Sample %in% refsample)]
    cells2 <- info$Cell[which(info$Sample %in% querysample)]

    train_rds = subset(alldata, cells = cells1)
    test_rds = subset(alldata, cells = cells2)

    train_sce = sceasy::convertFormat(train_rds, from="seurat", to="sce")
    test_sce = sceasy::convertFormat(test_rds, from="seurat", to="sce")

    train_ann = train_rds$annotation
    test_ann = test_rds$annotation

    train_sce$celltype = train_ann
    test_sce$celltype = test_ann

    #train_result <- SingleR(train_sce, train_sce, labels = train_sce$celltype)
    test_result = SingleR(test_sce, train_sce, labels = train_sce$celltype)

    saveRDS(test_result, paste('dataset_SingleR_', organ, "_", percent , '_K', kvalue,'.rds',sep = ''))

}



########################################################################################################################
# scmap
library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(dplyr)
library(readxl)
library(openxlsx)

source("process.R")

organ <- "Leaf"
alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis


percent = 0.05 
for(kvalue in 1:5){
    # kvalue <- 5

    downsample <- readRDS(sprintf("dataset_downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refsample <- downsample$ref
    querysample <- downsample$query

    cells1 <- info$Cell[which(info$Sample %in% refsample)]
    cells2 <- info$Cell[which(info$Sample %in% querysample)]

    ref_data = subset(alldata, cells = cells1)
    que_data = subset(alldata, cells = cells2)
  
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
    saveRDS(scmapCluster_results, paste('dataset_scmap_', organ, "_", percent, '_K', kvalue, '.rds', sep = ''))
}



##########################################################################################################################
# scClassify
library(scClassify)
library(Seurat)
library(scater)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(readxl)
library(openxlsx)

source("process.R")

organ <- "Leaf"
alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

percent = 0.1
for(kvalue in 1:5){
    # kvalue = 2

    downsample <- readRDS(sprintf("dataset_downsample_ids_%s_%s_%s.rds", organ, percent, kvalue))
    refsample <- downsample$ref
    querysample <- downsample$query

    cells1 <- info$Cell[which(info$Sample %in% refsample)]
    cells2 <- info$Cell[which(info$Sample %in% querysample)]

    ref_cds = subset(alldata, cells = cells1)
    que_cds = subset(alldata, cells = cells2)

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

    ref_anno = ref_cds$annotation
    que_anno = que_cds$annotation

    # train model
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

    scClassify_res$annotation <- que_anno

    saveRDS(scClassify_res, paste('dataset_scClassify_', organ, "_", percent, '_K', kvalue,'.rds', sep = ''))
    saveRDS(trainRes, paste('dataset_scClassify_', organ, "_", percent, '_K', kvalue, '_model.rds', sep = ''))

    print(kvalue)
}
