# rejection experiments
# scPlantGM 
library(openxlsx)
library(scPlantGM)

# scPlantGM
source("process.R")

species <- "ara"
res <- 1
topn <- 100
organ <- "Leaf" # c("Root", "Leaf")

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis
info1 <- info %>% filter(Organ==organ)

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

datasets <- unique(info1$Dataset)
allscore <- c()
for(id in 1:length(datasets)){

    seldataset <- datasets[id]
    info2 <- info1 %>% filter(Dataset %in% seldataset)
    table(info2$annotation)

    data <- read_excel(path="dataset_info/Cell_Type1.xlsx", sheet="ara")
    layer_info <- data.frame(data) %>% dplyr::filter(Organ==organ)

    if(organ == "Leaf"){
        refdata <- layer_info[which(layer_info$Celltype1 %in% c("Mesophyll","Leaf epidermis")),c("Celltype1","Celltype2","Celltype3")]
        reftype <- setdiff(unlist(refdata),"/")

        querydata <- layer_info[which(layer_info$Celltype1 %in% c("Vascular tissue")),c("Celltype1","Celltype2","Celltype3")]
        querytype <- setdiff(unlist(querydata),"/")
    }

    if(organ == "Root"){
        refdata <- layer_info[which(layer_info$Celltype1 %in% c("Root epidermis","Root cortex","Root cap")),c("Celltype1","Celltype2","Celltype3")]
        reftype <- setdiff(unlist(refdata),"/")

        querydata <- layer_info[which(layer_info$Celltype1 %in% c("Root stele")),c("Celltype1","Celltype2","Celltype3")]
        querytype <- setdiff(unlist(querydata),"/")
    }

    refcluster1 <- names(acc_list_sta)[which(names(acc_list_sta[[length(acc_list_sta)]]) %in% reftype)]
    refcluster <- intersect(refcluster1, info2$Cluster)
    querycluster1 <- names(acc_list_sta)[which(names(acc_list_sta[[length(acc_list_sta)]]) %in% querytype)]
    querycluster <- intersect(querycluster1, info2$Cluster)

    result <- predict_based_module_weight(info, jaccard_mat, layer=3,
                                acc_list_sta, acc_list_sta1, acc_list_sta2, acc_list_sta3,
                                refsample=refcluster, querysample=querycluster, 
                                gene_list, pthres = 0.8, x_thres=0.2,
                                issample = FALSE, iscluster = TRUE)

    saveRDS(result,sprintf("scPlantGM_%s_%s.rds", organ, seldataset))

    score <- length(which(result$prediction1=="Unknown"))/length(result$prediction1)
    allscore <- c(allscore,score)
}

allscore <- data.frame(Value = allscore)
write.xlsx(allscore, file = sprintf("scPlantGM_dataset_%s.xlsx", organ), rowNames = FALSE)



##############################################################################################################
# scmap
library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(dplyr)
library(readxl)
library(openxlsx)


source("process.R")

species <- "ara"
organ <- "Leaf"
alldata = readRDS('Arabidopsis_clean.rds')

data("reference_Arabidopsis", package = "scPlantGM")
jaccard_mat <- reference_Arabidopsis$sim_mat

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis
info1 <- info %>% filter(Organ==organ)

datasets <- unique(info1$Dataset)
allscore <- c()
for(id in 1:length(datasets)){

    seldataset <- datasets[id]
    info2 <- info1 %>% filter(Dataset %in% seldataset)
    table(info2$annotation)

    if(organ == "Leaf"){
        refdata <- layer_info[which(layer_info$Celltype1 %in% c("Mesophyll","Leaf epidermis")),c("Celltype1","Celltype2","Celltype3")]
        reftype <- setdiff(unlist(refdata),"/")

        querydata <- layer_info[which(layer_info$Celltype1 %in% c("Vascular tissue")),c("Celltype1","Celltype2","Celltype3")]
        querytype <- setdiff(unlist(querydata),"/")
    }

    if(organ == "Root"){
        refdata <- layer_info[which(layer_info$Celltype1 %in% c("Root epidermis","Root cortex","Root cap")),c("Celltype1","Celltype2","Celltype3")]
        reftype <- setdiff(unlist(refdata),"/")

        querydata <- layer_info[which(layer_info$Celltype1 %in% c("Root stele")),c("Celltype1","Celltype2","Celltype3")]
        querytype <- setdiff(unlist(querydata),"/")
    }

    refcluster <- info2 %>% filter(annotation %in% reftype) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- info2 %>% filter(annotation %in% querytype) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- intersect(querycluster, colnames(jaccard_mat))

    cells1 <- info$Cell[which(info$Cluster%in%refcluster)]
    cells2 <- info$Cell[which(info$Cluster%in%querycluster)]

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
        threshold = 0.7
    )

    scmapCluster_results$Cell <- colnames(que_matrix)
    scmapCluster_results$annotation <- que_data@meta.data$anno
    scmapCluster_results$prob <- get_scmap_probmatrix(que_sce, list(metadata(ref_sce)$scmap_cluster_index),threshold = 0.1)
    saveRDS(scmapCluster_results, paste('scmap_',organ, "_", seldataset, '_model.rds',sep = ''))

    if(T){
        pos1 <- which(scmapCluster_results$annotation%in%c("Palisade mesophyll","Spongy mesophyll"))
        if(length(pos1)>0){scmapCluster_results$annotation[pos1]<-"Mesophyll"}
        pos2 <- which(scmapCluster_results$scmap_cluster_labs%in%c("Palisade mesophyll","Spongy mesophyll"))
        if(length(pos2)>0){scmapCluster_results$scmap_cluster_labs[pos2]<-"Mesophyll"}
    }

    score <- length(which(scmapCluster_results$scmap_cluster_labs=="unassigned"))/length(scmapCluster_results$scmap_cluster_labs)
    allscore <- c(allscore,score)
}

allscore <- data.frame(Value = allscore)
write.xlsx(allscore, file = sprintf("scmap_dataset_%s.xlsx", organ), rowNames = FALSE)



##############################################################################################################
# scClaasify
library(scClassify)
library(Seurat)
library(scater)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(readxl)
library(openxlsx)

source("process.R")

species <- "ara"
organ <- "Leaf"
alldata = readRDS('Arabidopsis_clean.rds')

data("reference_Arabidopsis", package = "scPlantGM")
jaccard_mat <- reference_Arabidopsis$sim_mat

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis
info1 <- info %>% filter(Organ==organ)

datasets <- unique(info1$Dataset)
allscore <- c()
for(id in 1:length(datasets)){

    seldataset <- datasets[id]
    info2 <- info1 %>% filter(Dataset %in% seldataset)
    table(info2$annotation)

    if(organ == "Leaf"){
        refdata <- layer_info[which(layer_info$Celltype1 %in% c("Mesophyll","Leaf epidermis")),c("Celltype1","Celltype2","Celltype3")]
        reftype <- setdiff(unlist(refdata),"/")

        querydata <- layer_info[which(layer_info$Celltype1 %in% c("Vascular tissue")),c("Celltype1","Celltype2","Celltype3")]
        querytype <- setdiff(unlist(querydata),"/")
    }

    if(organ == "Root"){
        refdata <- layer_info[which(layer_info$Celltype1 %in% c("Root epidermis","Root cortex","Root cap")),c("Celltype1","Celltype2","Celltype3")]
        reftype <- setdiff(unlist(refdata),"/")

        querydata <- layer_info[which(layer_info$Celltype1 %in% c("Root stele")),c("Celltype1","Celltype2","Celltype3")]
        querytype <- setdiff(unlist(querydata),"/")
    }

    refcluster <- info2 %>% filter(annotation %in% reftype) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- info2 %>% filter(annotation %in% querytype) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- intersect(querycluster, colnames(jaccard_mat))

    cells1 <- info$Cell[which(info$Cluster%in%refcluster)]
    cells2 <- info$Cell[which(info$Cluster%in%querycluster)]

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

    ref_anno = ref_cds@meta.data$annotation
    que_anno = que_cds@meta.data$annotation

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

    saveRDS(scClassify_res, paste('scClassify_', organ, "_", seldataset, '.rds',sep = ''))
    saveRDS(trainRes, paste('scClassify_', organ, "_", seldataset, '_model.rds',sep = ''))

    score <- length(which(scClassify_res$pearson_WKNN_limma$predRes=="unassigned"))/length(scClassify_res$pearson_WKNN_limma$predRes)
    allscore <- c(allscore,score)
}

allscore <- data.frame(Value = allscore)
write.xlsx(allscore, file = sprintf("scClassify_dataset_%s.xlsx", organ), rowNames = FALSE)




###########################################################################################################################################
############################################################################################################################################
# rejection experiments (cross-tissue)
organ1 <- "Root"
organ2 <- "Leaf"
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis
info1 <- info %>% filter(Organ==organ1)
info2 <- info %>% filter(Organ==organ2)

datasets1 <- unique(info1$Dataset)
datasets2 <- unique(info2$Dataset)

set.seed(123)
id1 <- sample(1:length(datasets1), 10)
set.seed(456)
id2 <- sample(1:length(datasets2), 10)

# saveRDS(id1, "methods/rejection/Root_ref_id.rds")
# saveRDS(id2, "methods/rejection/Leaf_query_id.rds")


#---------------------------------------------------------------------------------------------------------------------------------------
# scPlantGM
library(openxlsx)

source("process.R")

res <- 1
topn <- 100
species <- "ara"

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


organ1 <- "Leaf"
organ2 <- "Root"
info1 <- info %>% filter(Organ==organ1)
info2 <- info %>% filter(Organ==organ2)

datasets1 <- unique(info1$Dataset)
datasets2 <- unique(info2$Dataset)

id1 <- readRDS(sprintf("%s_id.rds", organ1))
id2 <- readRDS(sprintf("%s_id.rds", organ2))

allscore <- c()
for(k in 1:10){

    seldataset1 <- datasets1[id1[k]]
    info11 <- info1 %>% filter(Dataset %in% seldataset1)
    table(info11$annotation)

    seldataset2 <- datasets2[id2[k]]
    info22 <- info2 %>% filter(Dataset %in% seldataset2)
    table(info22$annotation)

    refcluster <- info11 %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- info2 %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- intersect(querycluster, colnames(jaccard_mat))

    result <- predict_based_module_weight(info, jaccard_mat, layer=3,
                                acc_list_sta, acc_list_sta1, acc_list_sta2, acc_list_sta3,
                                refsample=refcluster, querysample=querycluster, 
                                gene_list, pthres = 0.8, x_thres=0.2,
                                issample = FALSE, iscluster = TRUE)
    saveRDS(result, sprintf("scPlantGM_%s_%s_%s_%s.rds", organ1, seldataset1, organ2, seldataset2))

    score <- length(which(result$prediction1=="Unknown"))/length(result$prediction1)
    allscore <- c(allscore,score)

    print(k)

}
allscore <- data.frame(Value = allscore)
write.xlsx(allscore, file = sprintf("scPlantGM_tissue_%s_%s.xlsx", organ1, organ2), rowNames = FALSE)



#####################################################################################################################################
# scmap
library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(dplyr)
library(readxl)
library(openxlsx)


source("process.R")

species <- "ara"
alldata = readRDS('Arabidopsis_clean.rds')

data("reference_Arabidopsis", package = "scPlantGM")
jaccard_mat <- reference_Arabidopsis$sim_mat

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

organ1 <- "Root"
organ2 <- "Leaf"
info1 <- info %>% filter(Organ==organ1)
info2 <- info %>% filter(Organ==organ2)

datasets1 <- unique(info1$Dataset)
datasets2 <- unique(info2$Dataset)

id1 <- readRDS("Root_ref_id.rds")
id2 <- readRDS("Leaf_query_id.rds")

allscore <- c()
for(k in 1:10){

    seldataset1 <- datasets1[id1[k]]
    info11 <- info1 %>% filter(Dataset %in% seldataset1)
    table(info11$annotation)

    seldataset2 <- datasets2[id2[k]]
    info22 <- info2 %>% filter(Dataset %in% seldataset2)
    table(info22$annotation)

    refcluster <- info11 %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- info2 %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- intersect(querycluster, colnames(jaccard_mat))

    cells1 <- info$Cell[which(info$Cluster%in%refcluster)]
    cells2 <- info$Cell[which(info$Cluster%in%querycluster)]

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
        threshold = 0.7
    )

    scmapCluster_results$Cell <- colnames(que_matrix)
    scmapCluster_results$annotation <- que_data@meta.data$anno
    scmapCluster_results$prob <- get_scmap_probmatrix(que_sce, list(metadata(ref_sce)$scmap_cluster_index),threshold = 0.1)
    saveRDS(scmapCluster_results, paste('scmap_', organ1, "_", seldataset1, "_", organ2, "_", seldataset2, '_model.rds',sep = ''))

    if(T){
        pos1 <- which(scmapCluster_results$annotation%in%c("Palisade mesophyll","Spongy mesophyll"))
        if(length(pos1)>0){scmapCluster_results$annotation[pos1]<-"Mesophyll"}
        pos2 <- which(scmapCluster_results$scmap_cluster_labs%in%c("Palisade mesophyll","Spongy mesophyll"))
        if(length(pos2)>0){scmapCluster_results$scmap_cluster_labs[pos2]<-"Mesophyll"}
    }

    score <- length(which(scmapCluster_results$scmap_cluster_labs=="unassigned"))/length(scmapCluster_results$scmap_cluster_labs)
    allscore <- c(allscore,score)
}

allscore <- data.frame(Value = allscore)
write.xlsx(allscore, file = sprintf("scmap_tissue_%s_%s.xlsx", organ1, organ2), rowNames = FALSE)



###############################################################################################################################
# scClaasify
library(scClassify)
library(Seurat)
library(scater)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(readxl)
library(openxlsx)

source("process.R")

species <- "ara"
alldata = readRDS('Arabidopsis_clean.rds')

data("reference_Arabidopsis", package = "scPlantGM")
jaccard_mat <- reference_Arabidopsis$sim_mat

data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

organ1 <- "Root"
organ2 <- "Leaf"
info1 <- info %>% filter(Organ==organ1)
info2 <- info %>% filter(Organ==organ2)

datasets1 <- unique(info1$Dataset)
datasets2 <- unique(info2$Dataset)

id1 <- readRDS("Root_ref_id.rds")
id2 <- readRDS("Leaf_query_id.rds")

allscore <- c()
for(k in 1:10){

    seldataset1 <- datasets1[id1[k]]
    info11 <- info1 %>% filter(Dataset %in% seldataset1)
    table(info11$annotation)

    seldataset2 <- datasets2[id2[k]]
    info22 <- info2 %>% filter(Dataset %in% seldataset2)
    table(info22$annotation)

    refcluster <- info11 %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- info2 %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    querycluster <- intersect(querycluster, colnames(jaccard_mat))

    cells1 <- info$Cell[which(info$Cluster%in%refcluster)]
    cells2 <- info$Cell[which(info$Cluster%in%querycluster)]

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

    ref_anno = ref_cds@meta.data$annotation
    que_anno = que_cds@meta.data$annotation

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

    saveRDS(scClassify_res, paste('scClassify_', organ1, "_", seldataset1, "_", organ2, "_", seldataset2, '.rds',sep = ''))
    saveRDS(trainRes, paste('scClassify_', organ1, "_", seldataset1, "_", organ2, "_", seldataset2, '_model.rds',sep = ''))

    score <- length(which(scClassify_res$pearson_WKNN_limma$predRes=="unassigned"))/length(scClassify_res$pearson_WKNN_limma$predRes)
    allscore <- c(allscore,score)
}

allscore <- data.frame(Value = allscore)
write.xlsx(allscore, file = sprintf("scClassify_tissue_%s_%s.xlsx", organ1, organ2), rowNames = FALSE)


