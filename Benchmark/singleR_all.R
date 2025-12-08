# Split by dataset-level cross-validation

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly= TRUE)
print(args)

library(SingleR)
library(Seurat)
library(magrittr)
library(ggplot2)

source("process.R")

species <- "ara"
organ = args[1]
k = args[2]

alldata = readRDS('Arabidopsis_clean.rds')
data("info_Arabidopsis", package = "scPlantGM")
info <- info_Arabidopsis

ref_list = get_sample(info, organ, k, sheet="ara_dataset")[[1]]
que_list = get_sample(info, organ, k, sheet="ara_dataset")[[2]]

train_rds = alldata[, alldata@meta.data$orig.ident %in% ref_list]
test_rds = alldata[, alldata@meta.data$orig.ident %in% que_list]

train_sce = sceasy::convertFormat(train_rds, from="seurat", to="sce")
test_sce = sceasy::convertFormat(test_rds, from="seurat", to="sce")

train_ann = train_rds@meta.data$annotation
test_ann = test_rds@meta.data$annotation

train_sce$celltype = train_ann
test_sce$celltype = test_ann


start_time <- Sys.time()

#train_result <- SingleR(train_sce, train_sce, labels = train_sce$celltype)
test_result = SingleR(test_sce, train_sce, labels = train_sce$celltype)

end_time <- Sys.time()
runtime <- end_time - start_time
run_time <- sprintf("%s %s runtime: %s", organ, k, runtime)
print(start_time)
print(end_time)
print(run_time)

saveRDS(test_result, paste('SingleR_result/SingleR_',organ,'_K',k,'_model.rds',sep = ''))


#result record
#train_sce$singler <- train_result[match(colnames(train_sce), rownames(train_result)), 'labels']
test_sce$singler <- test_result[match(colnames(test_sce), rownames(test_result)), 'labels']
saveRDS(test_sce, paste('SingleR_result/SingleR_',organ,'_K',k,'.rds',sep = ''))



#########################################################################################################
# Split by sample-level cross-validation
options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly= TRUE)
print(args)

library(SingleR)
library(Seurat)
library(magrittr)
library(readxl)

source("process.R")

species <- "ara"
organ = args[1]  # organ = "Root"
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

    train_rds = alldata[, alldata@meta.data$orig.ident %in% ref_list]
    test_rds = alldata[, alldata@meta.data$orig.ident %in% que_list]

    train_sce = sceasy::convertFormat(train_rds, from="seurat", to="sce")
    test_sce = sceasy::convertFormat(test_rds, from="seurat", to="sce")

    train_ann = train_rds@meta.data$annotation
    test_ann = test_rds@meta.data$annotation

    train_sce$celltype = train_ann
    test_sce$celltype = test_ann


    start_time <- Sys.time()

    #train_result <- SingleR(train_sce, train_sce, labels = train_sce$celltype)
    test_result = SingleR(test_sce, train_sce, labels = train_sce$celltype)

    end_time <- Sys.time()
    runtime <- end_time - start_time
    run_time <- sprintf("%s %s runtime: %s", organ, k, round(runtime,4))
    print(start_time)
    print(end_time)
    print(run_time)

    saveRDS(test_result, paste('SingleR_result/SingleR_sample_', organ, "_", dataset, '_K',k,'_model.rds',sep = ''))


    #result record
    #train_sce$singler <- train_result[match(colnames(train_sce), rownames(train_result)), 'labels']
    test_sce$singler <- test_result[match(colnames(test_sce), rownames(test_result)), 'labels']
    saveRDS(test_sce, paste('SingleR_result/SingleR_sample_', organ, "_", dataset,'_K', k, '.rds',sep = ''))

}
