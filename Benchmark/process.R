get_sample <- function(info, organ, kvalue, dataset = NULL, sheet = NULL){

    require(readxl)
    require(dplyr)
    data <- read_excel(path="ref_query.xlsx", sheet)
    data <- data.frame(data)

     
    if(is.null(dataset)){
      data1 <- data %>% filter(Organ==organ & K==kvalue)
    }else{
      data1 <- data %>% 
              filter(Dataset==dataset) %>% 
              filter(Organ==organ & K==kvalue) 
    }

    if(length(unique(data1$Dataset))>1){
        dataset1 <- data1$Dataset[which( data1$RefQuery=="ref")]
        dataset2 <- data1$Dataset[which( data1$RefQuery=="query")]
        sample1 <- info %>% filter(Organ==organ) %>% filter(Dataset %in% dataset1)  %>% select(Sample) %>% unique() %>% unlist() %>% as.character()
        sample2 <- info %>% filter(Organ==organ) %>% filter(Dataset %in% dataset2)  %>% select(Sample) %>% unique() %>% unlist() %>% as.character()
    }
    if(length(unique(data1$Dataset))==1){
        sample1 <- data1$Sample[which( data1$RefQuery=="ref")]
        sample2 <- data1$Sample[which( data1$RefQuery=="query")]
    }
    return(list(refsample=sample1,querysample=sample2))
}



predict_based_module_weight <- function(info, jaccard_mat, layer,
                                acc_list_sta, acc_list_sta1, acc_list_sta2, acc_list_sta3,
                                refsample, querysample, gene_list, pthres = 0.8, x_thres = 0.2,
                                issample = TRUE, iscluster = FALSE){
    require(dplyr)

    if(issample){
        cluster1 <- info %>% filter(Sample %in% refsample) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
        cluster2 <- info %>% filter(Sample %in% querysample) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
        cluster2 <- intersect(cluster2,colnames(jaccard_mat))
    }
    if(iscluster){
        cluster1 <- refsample
        cluster2 <- querysample
    }

    cells2 <- info %>% filter(Cluster %in% cluster2) %>% select(Cell) %>% unlist() %>% as.character()

    Cluster1 <- info$Cluster[match(cells2,info$Cell)]
    real1 <- info$Celltype1[match(cells2,info$Cell)]

    # purity_cluster <- names(acc_list_sta[-length(acc_list_sta)])
    jaccard_mat_remain1 <- jaccard_mat[which(colnames(jaccard_mat) %in% c(cluster1)), 
                                        which(colnames(jaccard_mat) %in% c(cluster1)), 
                                        drop = FALSE]

    ncelltype <- length(unique(info$annotation[which(info$Cluster %in% Cluster1)]))
    bestmodule <- get_module(sim_mat=jaccard_mat_remain1, acc_list_sta, thres=0.8, pthres, nm1=ncelltype)
    out.id <- bestmodule[[1]]
    k <- bestmodule[[2]]
    # print(k)

    acc_list_sta_new1 <- module_celltype(acc_list_sta1, out.id, k)
    acc_list_sta_new1 <- c(acc_list_sta_new1[[1]],list(acc_list_sta_new1[[2]]))

    if(T){
        marker_list_module <- list()
        for(j in 1:length(table(out.id))){
            clustername <- names(which(out.id==j))
            gene_sets1 <- list()
            for(jj in 1:length(clustername)){
                gene_sets1[[jj]] <- gene_list[[clustername[jj]]]
            }

            if(length(gene_sets1)==1){
                marker_list_module[[j]]  <- gene_sets1[[1]]
            }else{
                marker_list_module[[j]]  <- merged_weighted(gene_sets=gene_sets1)$order
            }
        }
        print("The organization of ordered marker list for gene modules has been finished!")

        if(k==nrow(jaccard_mat_remain1)){
            jaccard_mat1 <- jaccard_mat[which(colnames(jaccard_mat) %in% cluster2), which(colnames(jaccard_mat) %in% cluster1), drop = FALSE]
            colnames(jaccard_mat1) <- names(marker_list_module)[1:length(marker_list_module)]

        }else{
            jaccard_mat_re <- jaccard_mat[which(colnames(jaccard_mat) %in% cluster2), which(colnames(jaccard_mat) %in% cluster1), drop = FALSE]
            n1 <- length(marker_list_module)
            jaccard_mat1 <- matrix(0, nrow = nrow(jaccard_mat_re), ncol = n1, 
                    dimnames = list(rownames(jaccard_mat_re), names(acc_list_sta_new1)[-length(acc_list_sta_new1)]))

            # identical(names(out.id),colnames(jaccard_mat_re))
            for(i in 1:nrow(jaccard_mat_re)){
                x <- gene_list[[rownames(jaccard_mat_re)[i]]]
            
                o <- numeric(n1) 
                for(j in 1:n1) {
                    y <- marker_list_module[[j]]
                    z <- weighted_jaccard(x, y)
                    o[j] <- z
                }
                jaccard_mat1[i,] <- o
                # print(i)
            }
        }
        print("The calculation of jaccard matrix for gene modules has been finished!")
    }

    accuracy <- acc_list_sta_new1[[length(acc_list_sta_new1)]]
    types1 <- names(accuracy)

    pred_tmp <- cbind(rownames(jaccard_mat1), types1[apply(jaccard_mat1, 1, which.max)], apply(jaccard_mat1, 1, max))
    pred_tmp <- data.frame(pred_tmp)
    colnames(pred_tmp) <- c("Cluster","prediction", "similarity")

    result1 <- data.frame(Cell=cells2, Cluster=Cluster1, annotation1=real1, prediction1=rep(NA,length(cells2)), similarity1=rep(NA,length(cells2)))
    result1$prediction1 <- pred_tmp$prediction[match(result1$Cluster,pred_tmp$Cluster)]   
    result1$similarity1 <- round(as.numeric(pred_tmp$similarity[match(result1$Cluster, pred_tmp$Cluster)]),4) 
    result <- result1
    # length(which(result$annotation1==result$prediction1))/nrow(result)

    print("Layer1 has been finished!")

    if(layer>1){
        type1 <- names(table(result$prediction1))
        result2 <- c()

        acc_list_sta_new2 <- module_celltype(acc_list_sta2, out.id, k)
        acc_list_sta_new2 <- c(acc_list_sta_new2[[1]],list(acc_list_sta_new2[[2]]))

        accuracy2 <- acc_list_sta_new2[[length(acc_list_sta_new2)]]
        types2 <- names(accuracy2)

        if(length(setdiff(type1,"/"))>0){
            for(t in 1:length(type1)){
                cells22 <- result$Cell[which(result$prediction1==type1[t])]
                cluster22 <- info %>% filter(Cell %in% cells22) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
                cluster11 <- info %>% filter(Sample %in% refsample) %>% filter(Celltype1==type1[t]) %>% 
                                filter(layer>1) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
                cluster11 <- intersect(cluster11,colnames(jaccard_mat))

                Cluster2 <- info$Cluster[match(cells22,info$Cell)]
                real2 <- info$Celltype2[match(cells22,info$Cell)]

                jaccard_mat2 <- jaccard_mat1[Cluster2,which(types1==type1[t]), drop = FALSE]
                types22 <- types2[which(types1==type1[t])]

                pred_tmp <- cbind(rownames(jaccard_mat2), types22[apply(jaccard_mat2, 1, which.max)], apply(jaccard_mat2, 1, max))
                pred_tmp <- data.frame(pred_tmp)
                colnames(pred_tmp) <- c("Cluster", "prediction", "similarity")

                result_tmp <- data.frame(Cell=cells22, Cluster=Cluster2, annotation2=real2, prediction2=rep(NA,length(cells22)), similarity2=rep(NA,length(cells22)))
                result_tmp$prediction2 <- pred_tmp$prediction[match(result_tmp$Cluster,pred_tmp$Cluster)]
                result_tmp$similarity2 <- round(as.numeric(pred_tmp$similarity[match(result_tmp$Cluster, pred_tmp$Cluster)]),4)

                result2 <- rbind(result2, result_tmp)
            }
        }else{
            result2 <- as.data.frame(matrix(nrow = 0, ncol = 5))
            colnames(result2) <- c("Cell", "Cluster", "annotation2", "prediction2", "similarity2")
        }
        result <- merge(result, result2, by = c("Cell","Cluster"), all = TRUE)
    }
    
    print("Layer2 has been finished!")

    if(layer>2){  
        
        # type2 <- unique(c(names(table(result$annotation2)),names(table(result$prediction2))))
        type2 <- unique(names(table(result$prediction2)))

        if(!is.null(type2)){

            result3 <- c()

            acc_list_sta_new3 <- module_celltype(acc_list_sta3, out.id, k)
            acc_list_sta_new3 <- c(acc_list_sta_new3[[1]],list(acc_list_sta_new3[[2]]))

            accuracy3 <- acc_list_sta_new3[[length(acc_list_sta_new3)]]
            types3 <- names(accuracy3)
            
            # if(length(setdiff(types3,"/")) > 0){
            if(length(setdiff(type2,"/")) > 0){
                for(t in 1:length(type2)){
                    cells222 <- result$Cell[which(result$prediction2==type2[t])]
                    if(length(cells222)==0){next}
                    cluster222 <- info %>% filter(Cell %in% cells222) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()
                    cluster111 <- info %>% filter(Cluster %in% cluster1) %>% filter(Celltype2==type2[t]) %>% 
                                    filter(layer==3) %>% select(Cluster) %>% unique() %>% unlist() %>% as.character()

                    Cluster3 <- info$Cluster[match(cells222,info$Cell)]
                    real3 <- info$Celltype3[match(cells222,info$Cell)]

                    jaccard_mat3 <- jaccard_mat1[Cluster3,which(types2==type2[t]), drop = FALSE]
                    types33 <- types3[which(types2==type2[t])]
                    accuracy33 <- accuracy3[which(types2==type2[t])]

                    pred_tmp <- cbind(rownames(jaccard_mat3), types33[apply(jaccard_mat3, 1, which.max)], apply(jaccard_mat3, 1, max))
                    pred_tmp <- data.frame(pred_tmp)
                    colnames(pred_tmp) <- c("Cluster", "prediction", "similarity")

                    result_tmp <- data.frame(Cell=cells222, Cluster=Cluster3, annotation3=real3, prediction3=rep(NA,length(cells222)), similarity3=rep(NA,length(cells222)))
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
    }
    print("Layer3 has been finished!")

    prob_names = colnames(result %>% select(starts_with("similarity")))
    for (prob_name in prob_names) {

        loc_prob = which(colnames(result) == prob_name)
        if(length(which(result[loc_prob]=="/"))>0){
            pos1 = setdiff(which(result[loc_prob] < x_thres), which(result[loc_prob]=="/"))
        }else{
            pos1 = which(result[loc_prob] < x_thres)
        }
        result[pos1, (loc_prob - 1)] <- "Unknown"

        if(prob_name=="similarity1"){
            result[pos1, c("annotation2", "prediction2", "similarity2", "annotation3", "prediction3", "similarity3")] <- "/"
        }

        if(prob_name=="similarity2"){
            result[pos1, c("annotation3", "prediction3", "similarity3")] <- "/"
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