#' @title Predict cell types based on modules
#'
#' @param query A list of query Seurat object
#' @param reference A list of reference Seurat object (can be ignored if custom='None')
#' @param topmarkers_ref A list, top markers of each cluster obtained from Seurat pipeline for reference Seurat object
#' @param topmarkers_query A data.frame, top markers of each cluster obtained from Seurat pipeline for query Seurat object
#' @param species Species of scRNA-seq data (built-in species: Arabidopis, Maize, Rice)
#' @param organ Organ of species [Built-in organs: Arabidopis(Root,Leaf,Pollen,Inflorescence), Maize(Root,Leaf), Rice(Root,Leaf)]
#' @param custom Custom use built-in reference or not. Users can input 'All', 'Semi' or 'None'
#' @param layer_info A format-fixed data.frame includes layers information (can be ignored if custom='None')
#' @param layer Layer of prediction: must be an int (0,1,2...)
#' @param p_thres A threshold for determining what level of module purity is acceptable. The higher the value, the purer the module are
#' @param m_thres A threshold for determining how many modules can be accepted. The higher the value, the more modules can be accepted 
#' @param x_thres A threshold for determining whether the cell type should be ‘unassigned’
#' @param cores Number of CPU for running
#'
#' @return Prediction result in data.frame
#' @export
#'
#' @examples result -> scPlantGM(query, reference, species = 'Maize', organ = 'Root', layer = 0, custom = 'None')

scPlantGM <- function(query, reference, 
                        topmarkers_ref = NULL, topmarkers_query = NULL,
                        species, organ,
                        custom = 'None', layer_info = NA, layer = 0, 
                        p_thres = 0.8, m_thres = 0.8, x_thres = 0.2, 
                        ref_clustername = "seurat_clusters", query_clustername = "seurat_clusters",
                        topn = 100, cores = NA){
    require(dplyr)
    require(tidyr)
    require(doParallel)
    require(parallel)

    allcores <- detectCores()
    if (is.na(cores)){
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
    
    seuorder <- names(query)
    info_query <- get_info(seulist = query, type = 'query', query_clustername)

    if(!is.null(topmarkers_query)){

        if(is.data.frame(topmarkers_ref)){ 
            topmarkers_ref <- list(topmarkers_ref)
            names(topmarkers_ref) <- names(reference)
        }
        if(is.data.frame(topmarkers_query)){ 
            topmarkers_query <- list(topmarkers_query)
            names(topmarkers_query) <- names(query)
        }

        if(custom == 'None'){
            
            data(list=sprintf('reference_%s', species), package = "scPlantGM")
            # load(sprintf('data/reference_%s.rda', species))
            sim_mat_ref <- get(sprintf('reference_%s', species))$sim_mat
            gene_list_ref <- get(sprintf('reference_%s', species))$gene_list

            data(list=sprintf('info_%s', species), package = "scPlantGM")
            # load(sprintf('data/info_%s.rda', species))
            info_reference <- get(sprintf('info_%s', species))
            info_reference <- info_reference %>% filter(Organ==organ)

            # load layer info
            data(list=sprintf("layer_info_%s", species), package = "scPlantGM")
            # load(sprintf("data/layer_info_%s.rda", species))
            layer_info <- get(sprintf("layer_info_%s", species))
            layer_info <- layer_info_Arabidopsis %>% 
                        dplyr::filter(Organ==organ) %>% 
                        dplyr::select(Celltype1, Celltype2, Celltype3)

            allcluster <- intersect(unique(info_reference$Cluster), colnames(sim_mat_ref))
            sim_mat_ref <- sim_mat_ref[allcluster, allcluster]
            gene_list_ref <- gene_list_ref[allcluster]

        } else if(custom == 'All') {

            info_reference <- get_info(reference, type ='reference', clustername = ref_clustername)
            if (!all(is.na(layer_info))){
                info_reference <- get_layers(info_reference, info_layer=layer_info)
            }
            gene_list_ref <- get_gene_list(seulist = reference, seuorder = names(reference),
                                                 clustername = ref_clustername, topn, cores = cores, 
                                                 topmarkers_list = topmarkers_ref)
            sim_mat_ref <- get_sim_mat(gene_list1 = gene_list_ref, cores = cores) 

            # Ensure consistency of clusters
            sameclusters <- intersect(unique(info_reference$Cluster), colnames(sim_mat_ref))
            info_reference <- info_reference %>% filter(Cluster %in% sameclusters)
            gene_list_ref <- gene_list_ref[sameclusters]
            sim_mat_ref <- sim_mat_ref[sameclusters, sameclusters]

        } else if (custom == 'Semi') {

            info_reference1 <- get_info(reference, type ='reference', clustername = ref_clustername)
            gene_list_ref1 <- get_gene_list(seulist = reference, seuorder = names(reference),
                                                 clustername = ref_clustername,, topn, cores, 
                                                 topmarkers_list = topmarkers_ref)
            sim_mat_ref1 <- get_sim_mat(gene_list1 = gene_list_ref1, cores = cores) 

            # Ensure consistency of clusters
            sameclusters <- intersect(unique(info_reference1$Cluster), colnames(sim_mat_ref1))
            info_reference1 <- info_reference1 %>% filter(Cluster %in% sameclusters)
            gene_list_ref1 <- gene_list_ref1[sameclusters]
            sim_mat_ref1 <- sim_mat_ref1[sameclusters, sameclusters]

            data(list=paste('reference_', species, sep = ''), package = "scPlantGM")
            # load(paste('data/reference_', species, '.rda', sep = ''))
            sim_mat_ref2 <- get(paste('reference_', species, sep = ''))$sim_mat
            gene_list_ref2 <- get(paste('reference_', species, sep = ''))$gene_list

            data(list=paste('info_', species, sep = ''), package = "scPlantGM")
            # load(paste('scPlantGM/data/info_', species,'.rda', sep = ''))
            info_reference2 <- get(paste('info_', species, sep = ''))
            info_reference2 <- info_reference2 %>% filter(Organ==organ) %>% 
                              select('Sample', 'Annotation', 'Cell','Cluster') %>% 
                              filter(Annotation %in% setdiff(unique(unlist(layer_info)),'/'))
            info_reference <- rbind(info_reference1, info_reference2)
            if (!all(is.na(layer_info))){
                info_reference <- get_layers(info_reference, layer_info)
            }

            allcluster <- intersect(unique(info_reference$Cluster), colnames(sim_mat_ref2))
            sim_mat_ref2 <- sim_mat_ref2[allcluster, allcluster]
            sim_mat_ref3 <- get_sim_mat(gene_list_ref1, gene_list_ref2[allcluster], cores = cores)

            allcluster <- union(allcluster, colnames(sim_mat_ref1))       
            sim_mat_ref <- matrix(NA, nrow = length(allcluster), ncol = length(allcluster), 
                                    dimnames = list(allcluster, allcluster))
            sim_mat_ref[rownames(sim_mat_ref1), colnames(sim_mat_ref1)] <- sim_mat_ref1
            sim_mat_ref[rownames(sim_mat_ref2), colnames(sim_mat_ref2)] <- sim_mat_ref2
            sim_mat_ref[rownames(sim_mat_ref3), colnames(sim_mat_ref3)] <- sim_mat_ref3

            gene_list_ref <- c(gene_list_ref1,gene_list_ref2)[allcluster]

        } else {
            stop('Please input a correct value for custom!')
        }
        message('Identifying top markers for reference clusters... done')

        gene_list_query <- get_gene_list(seulist = query, seuorder = names(query), 
                                              clustername = query_clustername, topn, cores,
                                              topmarkers_list = topmarkers_query) 
        # sim_mat_que <- get_sim_mat(gene_list1 = gene_list_query, cores)
        message('Identifying top markers for query clusters... done')

    }else{
        if(custom == 'None'){
            
            data(list=sprintf('reference_%s', species), package = "scPlantGM")
            # load(sprintf('data/reference_%s.rda', species))
            sim_mat_ref <- get(sprintf('reference_%s', species))$sim_mat
            gene_list_ref <- get(sprintf('reference_%s', species))$gene_list

            data(list=sprintf('info_%s', species), package = "scPlantGM")
            # load(sprintf('data/info_%s.rda', species))
            info_reference <- get(sprintf('info_%s', species))
            info_reference <- info_reference %>% filter(Organ==organ)

            # load layer info
            data(list=sprintf("layer_info_%s", species), package = "scPlantGM")
            # load(sprintf("data/layer_info_%s.rda", species))
            layer_info <- get(sprintf("layer_info_%s", species))
            layer_info <- layer_info_Arabidopsis %>% 
                        dplyr::filter(Organ==organ) %>% 
                        dplyr::select(Celltype1, Celltype2, Celltype3)

            allcluster <- intersect(unique(info_reference$Cluster), colnames(sim_mat_ref))
            sim_mat_ref <- sim_mat_ref[allcluster, allcluster]
            gene_list_ref <- gene_list_ref[allcluster]

        } else if(custom == 'All') {

            info_reference <- get_info(reference, type ='reference', clustername = ref_clustername)
            if (!all(is.na(layer_info))){
                info_reference <- get_layers(info_reference, info_layer=layer_info)
            }
            gene_list_ref <- get_gene_list(seulist = reference, seuorder = names(reference), 
                                            clustername = ref_clustername, topn, cores) 
            sim_mat_ref <- get_sim_mat(gene_list1 = gene_list_ref, cores = cores) 

        } else if (custom == 'Semi') {

            info_reference1 <- get_info(reference, type ='reference', clustername = ref_clustername)
            gene_list_ref1 <- get_gene_list(seulist = reference, seuorder = names(reference), 
                                            clustername = ref_clustername, topn, cores) 
            sim_mat_ref1 <- get_sim_mat(gene_list1 = gene_list_ref1, cores = cores) 

            data(list=paste('reference_', species, sep = ''), package = "scPlantGM")
            # load(paste('data/reference_', species, '.rda', sep = ''))
            sim_mat_ref2 <- get(paste('reference_', species, sep = ''))$sim_mat
            gene_list_ref2 <- get(paste('reference_', species, sep = ''))$gene_list

            data(list=paste('info_', species, sep = ''), package = "scPlantGM")
            # load(paste('scPlantGM/data/info_', species,'.rda', sep = ''))
            info_reference2 <- get(paste('info_', species, sep = ''))
            info_reference2 <- info_reference2 %>% filter(Organ==organ) %>% 
                              select('Sample', 'Annotation', 'Cell','Cluster') %>% 
                              filter(Annotation %in% setdiff(unique(unlist(layer_info)),'/'))
            info_reference <- rbind(info_reference1, info_reference2)
            if (!all(is.na(layer_info))){
                info_reference <- get_layers(info_reference, layer_info)
            }

            allcluster <- intersect(unique(info_reference$Cluster), colnames(sim_mat_ref2))
            sim_mat_ref2 <- sim_mat_ref2[allcluster, allcluster]
            sim_mat_ref3 <- get_sim_mat(gene_list_ref1, gene_list_ref2[allcluster], cores = cores)

            allcluster <- union(allcluster, colnames(sim_mat_ref1))       
            sim_mat_ref <- matrix(NA, nrow = length(allcluster), ncol = length(allcluster), 
                                    dimnames = list(allcluster, allcluster))
            sim_mat_ref[rownames(sim_mat_ref1), colnames(sim_mat_ref1)] <- sim_mat_ref1
            sim_mat_ref[rownames(sim_mat_ref2), colnames(sim_mat_ref2)] <- sim_mat_ref2
            sim_mat_ref[rownames(sim_mat_ref3), colnames(sim_mat_ref3)] <- sim_mat_ref3

            gene_list_ref <- c(gene_list_ref1,gene_list_ref2)[allcluster]

        } else {
            stop('Please input a correct value for custom!')
        }
        message('Identifying top markers for reference clusters... done')

        gene_list_query <- get_gene_list(seulist = query, seuorder, clustername = query_clustername, topn, cores) 
        # sim_mat_que <- get_sim_mat(gene_list1 = gene_list_query, cores)
        message('Identifying top markers for query clusters... done')

    }

    sim_mat_cross <- get_sim_mat(gene_list1 = gene_list_query, gene_list2 = gene_list_ref, cores) 
    message('Calculating similarity matrix... done')


    cells2 <- info_query %>% filter(Cluster %in% rownames(sim_mat_cross)) %>% 
                select(Cell) %>% unlist() %>% as.character()
    Cluster1 <- info_query$Cluster[match(cells2, info_query$Cell)]

    if (custom == 'None'){ 
        data(list=paste('acc_list_sta_all_', species, sep = ''), package = "scPlantGM")
        # load(paste('data/acc_list_sta_all_', species, '.rda', sep = ''))
        acc_list_sta_all <- get(paste('acc_list_sta_all_', species, sep = ''))
    } else {
        acc_list_sta_all <- get_cluster_ratio(sim_mat_ref, info_reference, layer = layer)
    }
    acc_list_sta <- acc_list_sta_all$layer0

    ncelltype <- length(unique(info_reference$Annotation))
    bestmodule <- get_module(sim_mat = sim_mat_ref, acc_list_sta, m_thres, p_thres, min_nm = ncelltype)
    out.id <- bestmodule$ids
    k <- bestmodule$best_nm
    message("Projecting ", k, " modules to cells... done")

    
    if (layer==0){
        acc_list_sta_new1 <- module_celltype(acc_list_sta_all$layer0, out.id, k)
        layer_init = 0
    } else {
        acc_list_sta_new1 <- module_celltype(acc_list_sta_all$layer1, out.id, k)
        layer_init = 1
    }
    acc_list_sta_new1 <- c(acc_list_sta_new1[[1]], list(acc_list_sta_new1[[2]]))

    # Merge genes of the same module
    marker_list_module <- list()
    for(j in 1:length(table(out.id))){
        clustername <- names(which(out.id==j))
        gene_sets1 <- list()
        for(jj in 1:length(clustername)){
            gene_sets1[[jj]] <- gene_list_ref[[clustername[jj]]]
        }

        if(length(gene_sets1)==1){
            marker_list_module[[j]]  <- gene_sets1[[1]]
        }else{
            marker_list_module[[j]]  <- merged_weighted(gene_sets=gene_sets1)$order
        }
    }
    names(marker_list_module) <- names(acc_list_sta_new1)[-length(acc_list_sta_new1)]
    print("The organization of ordered marker list for gene modules has been finished!")

    if(k == ncol(sim_mat_cross)){
        sim_mat_cross1 <- sim_mat_cross
        colnames(sim_mat_cross1) <- names(marker_list_module)[1:length(marker_list_module)]
    }else{
        n1 <- length(marker_list_module)
        sim_mat_cross1 <- matrix(0, nrow = nrow(sim_mat_cross), ncol = n1, 
                dimnames = list(rownames(sim_mat_cross), names(acc_list_sta_new1)[-length(acc_list_sta_new1)]))

        # identical(names(out.id), colnames(sim_mat_cross11))
        for(i in 1:nrow(sim_mat_cross1)){
            x <- gene_list_query[[rownames(sim_mat_cross1)[i]]]
            o <- numeric(n1) 
            for(j in 1:n1) {
                y <- marker_list_module[[j]]
                z <- sim_score(x, y)
                o[j] <- z
            }
            sim_mat_cross1[i,] <- o
            # print(i)
        }
    }
    print("The calculation of jaccard matrix for gene modules has been finished!")
        
    accuracy <- acc_list_sta_new1[[length(acc_list_sta_new1)]]
    types1 <- names(accuracy)

    pred_tmp <- cbind(rownames(sim_mat_cross1), types1[apply(sim_mat_cross1, 1, which.max)], apply(sim_mat_cross1, 1, max))
    pred_tmp <- data.frame(pred_tmp)
    colnames(pred_tmp) <- c("Cluster","prediction", "similarity")

    result <- data.frame(Cell = cells2, Cluster = Cluster1)
    result$prediction1 <- pred_tmp$prediction[match(result$Cluster, pred_tmp$Cluster)]   
    result$similarity1 <- round(as.numeric(pred_tmp$similarity[match(result$Cluster, pred_tmp$Cluster)]),4) 
    print(paste('Layer ', layer_init, ' DONE!',sep = '')) 

    if(layer>1){

        type1 <- names(table(result$prediction1))
        result2 <- c()

        acc_list_sta_new2 <- module_celltype(acc_list_sta_all$layer2, out.id, k)
        acc_list_sta_new2 <- c(acc_list_sta_new2[[1]],list(acc_list_sta_new2[[2]]))

        accuracy2 <- acc_list_sta_new2[[length(acc_list_sta_new2)]]
        types2 <- names(accuracy2)

        if(length(setdiff(type1,"/")) > 0){
            for(t in 1:length(type1)){
                cells22 <- result$Cell[which(result$prediction1==type1[t])]
                Cluster2 <- info_query$Cluster[match(cells22,info_query$Cell)]

                sim_mat_cross2 <- sim_mat_cross1[Cluster2, which(types1==type1[t]), drop = FALSE]
                types22 <- types2[which(types1==type1[t])]

                pred_tmp <- cbind(rownames(sim_mat_cross2), types22[apply(sim_mat_cross2, 1, which.max)], apply(sim_mat_cross2, 1, max))
                pred_tmp <- data.frame(pred_tmp)
                colnames(pred_tmp) <- c("Cluster", "prediction", "similarity")

                result_tmp <- data.frame(Cell=cells22, Cluster=Cluster2)
                result_tmp$prediction2 <- pred_tmp$prediction[match(result_tmp$Cluster,pred_tmp$Cluster)]
                result_tmp$similarity2 <- round(as.numeric(pred_tmp$similarity[match(result_tmp$Cluster, pred_tmp$Cluster)]),4)

                result2 <- rbind(result2, result_tmp)
            }
        }else{
            result2 <- as.data.frame(matrix(nrow = 0, ncol = 5))
            colnames(result2) <- c("Cell", "Cluster", "annotation2", "prediction2", "similarity2")
        }
        result <- merge(result, result2, by = c("Cell","Cluster"), all = TRUE)
        print("Layer 2 DONE!")
    }
    

    if(layer>2){  
        
        # type2 <- unique(c(names(table(result$annotation2)),names(table(result$prediction2))))
        type2 <- unique(names(table(result$prediction2)))

        if(!is.null(type2)){

            result3 <- c()

            acc_list_sta_new3 <- module_celltype(acc_list_sta_all$layer3, out.id, k)
            acc_list_sta_new3 <- c(acc_list_sta_new3[[1]],list(acc_list_sta_new3[[2]]))

            accuracy3 <- acc_list_sta_new3[[length(acc_list_sta_new3)]]
            types3 <- names(accuracy3)
            
            # if(length(setdiff(types3,"/")) > 0){
            if(length(setdiff(type2,"/")) > 0){
                for(t in 1:length(type2)){
                    cells222 <- result$Cell[which(result$prediction2==type2[t])]
                    if(length(cells222)==0){next}
                    Cluster3 <- info_query$Cluster[match(cells222, info_query$Cell)]

                    sim_mat_cross3 <- sim_mat_cross1[Cluster3, which(types2==type2[t]), drop = FALSE]
                    types33 <- types3[which(types2==type2[t])]
                    accuracy33 <- accuracy3[which(types2==type2[t])]

                    pred_tmp <- cbind(rownames(sim_mat_cross3), types33[apply(sim_mat_cross3, 1, which.max)], apply(sim_mat_cross3, 1, max))
                    pred_tmp <- data.frame(pred_tmp)
                    colnames(pred_tmp) <- c("Cluster", "prediction", "similarity")

                    result_tmp <- data.frame(Cell=cells222, Cluster=Cluster3)
                    result_tmp$prediction3 <- pred_tmp$prediction[match(result_tmp$Cluster,pred_tmp$Cluster)]
                    result_tmp$similarity3 <- round(as.numeric(pred_tmp$similarity[match(result_tmp$Cluster, pred_tmp$Cluster)]),4)

                    result3 <- rbind(result3, result_tmp)
                }
                result <- merge(result, result3, by =  c("Cell","Cluster"), all = TRUE)

            }else{
                result$annotation3 <- "/"
                result$prediction3 <- "/"
                result$similarity3 <- result$similarity2              
            }
        
        }
        print("Layer 3 DONE!")
    }
    
    prob_names = colnames(result %>% select(starts_with("similarity")))
    for (prob_name in prob_names) {

        loc_prob = which(colnames(result) == prob_name)
        if(length(which(result[loc_prob]=="/")) > 0){
            pos1 = setdiff(which(result[loc_prob] < x_thres), which(result[loc_prob]=="/"))
        }else{
            pos1 = which(result[loc_prob] < x_thres)
        }
        result[pos1, (loc_prob - 1)] <- "unassigned"

        if(prob_name == "similarity3"){
            result[pos1, c("prediction3", "similarity3")] <- "unassigned"
        }

        if(prob_name == "similarity2"){
            result[pos1, c("prediction2", "similarity2")] <- "unassigned"
        }
    }
    
    predpos <- grep("prediction", colnames(result)) 
    if(length(predpos)==2){
        result$prediction <- result$prediction2
        result$prediction[which(result$prediction2=="/")] <- result$prediction1[which(result$prediction2=="/")]
    }
    if(length(predpos)==3){
        result$prediction <- result$prediction3
        result$prediction[which(result$prediction3=="/")] <- result$prediction2[which(result$prediction3=="/")] 
        result$prediction[which(result$prediction2=="/")] <- result$prediction1[which(result$prediction2=="/")] 
    }
    return(result)                        
}