#' @title Predict cell types based on modules
#'
#' @param query a list of query Seurat object
#' @param reference a list of reference Seurat object (can be ignored if custom='None')
#' @param species species of scRNA-seq data (built-in species: Arabidopis, Maize, Rice)
#' @param organ organ of species [Built-in organs: Arabidopis(Root,Leaf,Pollen,Inflorescence), Maize(Root,Leaf), Rice(Root,Leaf)]
#' @param custom custom use built-in reference or not. Users can input 'All', 'Semi' or 'None'
#' @param layer_info a format-fixed data.frame includes layers information (can be ignored if custom='None')
#' @param layer layer of prediction: must be an int (0,1,2...)
#' @param p_thres a threshold for determining what level of module purity is acceptable. The higher the value, the purer the module are
#' @param m_thres a threshold for determining how many modules can be accepted. The higher the value, the more modules can be accepted 
#' @param x_thres a threshold for determining whether the cell type should be ‘unknown’
#' @param cores numbers of CPU for running
#'
#' @return prediction result in data.frame
#' @export
#'
#' @examples result -> scPlantGM(query, reference, species = 'Maize', organ = 'Root', layer = 0, custom = 'None')

scPlantGM <- function(query, reference, species, organ,
                                custom = 'None', layer_info=NA, layer=0, 
                                p_thres=0.8, m_thres=0.9, x_thres=0.1, cores=NA){
    require(dplyr)
    require(tidyr)
    require(doParallel)
    require(parallel)
    require(doParallel)

    allcores <- detectCores()
    if (is.na(cores)==TRUE){
        cores <- 1
    } else if (0<cores & cores<=allcores) {
        cores <- cores
        message("Using ", cores, " cores for parallel computation")
    } else if(cores>allcores){
        warning('Input cores exceed available system cores. Default cores will be used!')
        cores <- ceiling(allcores/2)
    } else {
        stop('Please input a correct value for cores!')
    }
    
    info_query <- get_info(query, 'query')

    if(custom == 'None'){
        data(list=sprintf('jaccard.mat_%s', species), package = "scPlantGM")
        # load(paste('scPlantGM/data/jaccard.mat_', species, '.rda', sep = ''))
        jaccard_mat_ref <- get(paste('jaccard.mat_', species, sep = '')[1])
        data(list=paste('info_', species, sep = ''), package = "scPlantGM")
        # load(paste('scPlantGM/data/info_', species,'.rda', sep = ''))
        info_reference <- get(paste('info_', species, sep = ''))
        info_reference <- info_reference %>% filter(Organ==organ)

        # load layer info
        data("layer_info_Arabidopis", package = "scPlantGM")
        layer_info <- layer_info_Arabidopis %>% dplyr::filter(Organ==organ) %>% dplyr::select(Celltype1,Celltype2,Celltype3)

    } else if(custom == 'All') {
        info_reference <- get_info(reference, type ='reference')
        if (!all(is.na(layer_info))){
            info_reference <- get_layers(info_reference=info_reference, info_layer=layer_info)
        }
        jaccard_mat_ref <- get_jaccardmat(reference,type = 'reference', cores)
    } else if (custom == 'Semi') {
        info_reference1 <- get_info(reference, type='reference')
        jaccard_mat_ref1 <- get_jaccardmat(reference,'reference',cores)

        data(list=paste('jaccard.mat_', species, sep = ''), package = "scPlantGM")
        # load(paste('scPlantGM/data/jaccard.mat_', species, '.rda', sep = ''))
        jaccard_mat_ref2 <- get(paste('jaccard.mat_', species, sep = '')[1])
        data(list=paste('info_', species, sep = ''), package = "scPlantGM")
        # load(paste('scPlantGM/data/info_', species,'.rda', sep = ''))
        info_reference2 <- get(paste('info_', species, sep = ''))
        info_reference2 <- info_reference2 %>% filter(Organ==organ) %>% 
                           select('Sample', 'Annotation', 'Cell','Cluster') %>% 
                           filter(Annotation %in% setdiff(unique(unlist(layer_info)),'/'))
        info_reference <- rbind(info_reference1,info_reference2)
        if (!all(is.na(layer_info))){
            info_reference <- get_layers(info_reference,layer_info)
        }

        jaccard_mat_ref <- fuse_ref_jm(jaccard_mat_ref1,jaccard_mat_ref2,info_reference2,cores)
    } else {
        stop('Please input a correct value for custom!')
    }
    message('Identifying top markers for reference clusters... done')

    jaccard_mat_que <- get_jaccardmat(query, 'query',cores)
    message('Identifying top markers for query clusters... done')
    #jaccard_mat_que = readRDS('/public/workspace/jiangz/xuankun/wrapR/test/newdata_jm.rds')
    jaccard_mat <- fuse_jaccardmat(jaccard_mat_que, jaccard_mat_ref, cores)
    message('Calculating Jaccard matrix... done')

    #refsample <- info_reference$Sample
    cluster1 <- info_reference %>% filter(Annotation %in% as.character(unlist(layer_info))) %>% 
                select(Cluster) %>% unique() %>% unlist() %>% as.character() 
    cluster2 <- info_query %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
    cluster2 <- intersect(cluster2,rownames(jaccard_mat))

    cells2 <- info_query %>% filter(Cluster %in% cluster2) %>% select(Cell) %>% unlist() %>% as.character()

    Cluster1 <- info_query$Cluster[match(cells2,info_query$Cell)]

    if (custom=='None'){
        data(list=paste('acc_list_sta_all_', species, sep = ''), package = "scPlantGM")
        # load(paste('scPlantGM/data/acc_list_sta_all_', species, '.rda', sep = ''))
        acc_list_sta_all <- get(paste('acc_list_sta_all_', species, sep = ''))
    } else {
        acc_list_sta_all <- get_cluster_ratio(jaccard_mat_ref, info_reference, layer)
    }
    acc_list_sta <- acc_list_sta_all[[1]]

    jaccard_mat_ref1 <- jaccard_mat_ref[[1]]
    jaccard_mat_remain1 <- jaccard_mat_ref1[which(colnames(jaccard_mat_ref1) %in% cluster1), which(colnames(jaccard_mat_ref1) %in% cluster1), drop = FALSE]
    bestmodule <- get_module(jaccard_mat=jaccard_mat_remain1,acc_list_sta,p_thres,m_thres)
    out.id <- bestmodule[[1]]
    k <- bestmodule[[2]]
    message("Projecting ", k, " modules to cells... done")

    if (layer==0){
        acc_list_sta_new0 <- module_celltype(acc_list_sta_all[[1]], out.id, k)
        layer_init = 0
    } else {
        acc_list_sta_new0 <- module_celltype(acc_list_sta_all[[2]], out.id, k)
        layer_init = 1
    }
    acc_list_sta_new0 <- c(acc_list_sta_new0[[1]],list(acc_list_sta_new0[[2]]))

    jaccard_mat_re <- jaccard_mat[which(rownames(jaccard_mat) %in% cluster2), which(colnames(jaccard_mat) %in% cluster1), drop = FALSE]
    jaccard_mat0 <- c()
    for(j in 1:length(table(out.id))){
        tmpmat <- jaccard_mat_re[,names(which(out.id==j)), drop = FALSE]
        jaccard_mat0 <- cbind(jaccard_mat0,rowMeans(tmpmat))
    }
    colnames(jaccard_mat0) <- names(acc_list_sta_new0)[-length(acc_list_sta_new0)]
    
    accuracy <- acc_list_sta_new0[[length(acc_list_sta_new0)]]
    types0 <- names(accuracy)

    pred_tmp <- cbind(rownames(jaccard_mat0), types0[apply(jaccard_mat0, 1, which.max)], accuracy[apply(jaccard_mat0, 1, which.max)])
    pred_tmp <- data.frame(pred_tmp)
    colnames(pred_tmp) <- c("Cluster","prediction","probility")

    result <- data.frame(Cell=cells2, Cluster=Cluster1)
    result$prediction0 <- pred_tmp$prediction[match(result$Cluster,pred_tmp$Cluster)]
    result$probility0 <- pred_tmp$probility[match(result$Cluster,pred_tmp$Cluster)]
    colnames(result) <- c("Cell","Cluster",paste('prediction',as.character(layer_init),sep=''),paste('probility',as.character(layer_init),sep=''))
    print(paste('Layer ',layer_init, ' DONE!',sep = ''))     

    if (layer>1){
      jaccard_mat1 <- jaccard_mat0
      for (layer_num in 2:layer){
        acc_list_sta2 <- acc_list_sta_all[[layer_num+1]]
        loc_celltype <- which(colnames(info_reference)==paste('Celltype',as.character(layer_num),sep=''))
        info_celltype <- info_reference[[loc_celltype]]
        loc_predition <- which(colnames(result)==paste('prediction',as.character(layer_num-1),sep=''))
        info_prediction <- result[[loc_predition]]

        type1 <- unique(names(table(info_prediction)))
        result2 <- c()

        acc_list_sta_new2 <- module_celltype(acc_list_sta2, out.id, k)
        acc_list_sta_new2 <- c(acc_list_sta_new2[[1]],list(acc_list_sta_new2[[2]]))

        accuracy2 <- acc_list_sta_new2[[length(acc_list_sta_new2)]]
        types2 <- names(accuracy2)

        for(t in 1:length(type1)){
            cells22 <- result$Cell[which(info_prediction==type1[t])]
            if(length(cells22)==0){next}

            Cluster2 <- info_query$Cluster[match(cells22,info_query$Cell)]

            jaccard_mat2 <- jaccard_mat1[Cluster2,which(types0==type1[t]), drop = FALSE]
            types22 <- types2[which(types0==type1[t])]
            accuracy22 <- accuracy2[which(types0==type1[t])]

            pred_tmp <- cbind(rownames(jaccard_mat2),types22[apply(jaccard_mat2, 1, which.max)],accuracy22[apply(jaccard_mat2, 1, which.max)])
            pred_tmp <- data.frame(pred_tmp)
            colnames(pred_tmp) <- c("Cluster","prediction","probility")

            result_tmp <- data.frame(Cell=cells22, Cluster=Cluster2)
            result_tmp$prediction_new <- pred_tmp$prediction[match(result_tmp$Cluster,pred_tmp$Cluster)]
            result_tmp$probility_new <- pred_tmp$probility[match(result_tmp$Cluster,pred_tmp$Cluster)]
            colnames(result_tmp) <- c("Cell","Cluster",paste('prediction',as.character(layer_num),sep=''),paste('probility',as.character(layer_num),sep=''))

            result2 <- rbind(result2, result_tmp)
        }
        types0 <- types2
        result <- merge(result, result2, by = c("Cell","Cluster"), all = TRUE)
        print(paste('Layer',as.character(layer_num),'DONE!',sep=' '))
      }
    }
    result <- left_join(info_query['Cell'], result, by = 'Cell')
    result <- result %>% tidyr::separate(Cell, into=c('Sample','Cell'),sep=':') %>% select(-Sample)

    prob_names = colnames(result %>% select(starts_with('probility')))
    for (prob_name in prob_names){
        loc_prob = which(colnames(result)==prob_name)
        result[which(result[loc_prob]<x_thres),(loc_prob-1)] <- 'Unknown'
    }

    result <- result %>% select(-starts_with('probility'))
    result$prediction <- result$prediction2
    result$prediction[which(result$prediction2=="/")] <- result$prediction1[which(result$prediction2=="/")]

    return(result)                        
}