#' @title get markers from Seurat object
#'
#' @param object a Seurat object
#' @param clustername Cluster column name used
#' @param topn Number of top markers of each cluster
#' @param cores Number of CPU for running
#'
#' @return markers
#'
#' @export
get_markers <- function(object, clustername, topn, cores, avg.log2fc.threshold = 0.25){
  
    require(dplyr)
    require(doMC)

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    #object$seurat_clusters <- object$integrated_snn_res.0.8
    Idents(object) <- clustername
    csx <- table(unlist(object[[clustername]]))
    if(all(csx > 5)){
        top.markers <- foreach(i = sort(unique(unlist(object[[clustername]]))),
                              .combine=rbind, 
                              .packages = "Seurat",
                              .verbose = FALSE) %dopar% {
            if(csx[as.character(i)]>5){
            o <- FindMarkers(object, ident.1=i, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1,
                            verbose = FALSE)
            data.frame(o, gene=rownames(o), cluster=i)
            }else{
            NULL
            }
        }
    }else{
        top.markers <- FindAllMarkers(object, min.pct=0.1, logfc.threshold=0.1,
                                        return.thresh=0.1, only.pos=TRUE, 
                                        verbose = FALSE) ## return.thresh=0.01, test.use="wilcox"
    }
    stopCluster(cl)

    top.markers$pct.diff <- top.markers$pct.1 - top.markers$pct.2
    if(is.null(object@misc$geneName)){
      top.markers$name <- top.markers$gene
    }else{
      top.markers$name <- object@misc$geneName[top.markers$gene]
    }    

    ## keep top 100 marker
    topmarkers <- top.markers %>% 
      filter(p_val_adj < 0.05 & pct.1>0.2 & avg_log2FC > avg.log2fc.threshold & !grepl("^MT-", name)) %>% 
      dplyr::group_by(cluster)  %>% 
      dplyr::top_n(n=topn, wt=avg_log2FC) 
    
    return(topmarkers)
}


#' @title Calculate ordered gene list of clusters based on Seurat object
#'
#' @param seulist A list of Seurat objects or a Seurat object
#' @param seuorder A vector that show the order of Suerat in the list
#' @param clustername Cluster column name used
#' @param topn Number of top markers of each cluster
#' @param topmarkers A list, top markers of each cluster obtained from Seurat pipeline
#' 
#' @return Ordered gene list
#'
#' @export

get_gene_list <- function(seulist, seuorder, clustername, topn, cores, 
                          topmarkers_list = NULL, avg.log2fc.threshold = 0.25){

    require(Seurat)
    require(dplyr)
    require(data.table)

    if(is.null(topmarkers_list)){
        sdat <- lapply(seq_along(seulist), function(x) {
          topmarkers <- get_markers(object = seulist[[x]], clustername, topn, cores, avg.log2fc.threshold = 0.25)
          topmarkers$clusterID <- paste(seuorder[x], topmarkers$cluster, sep=":")
          topmarkers
        }
      )
    }else{
        sdat <- lapply(seq_along(seulist), function(x) {
          topmarkers <- topmarkers_list[[x]]
          topmarkers$clusterID <- paste(seuorder[x], topmarkers$cluster, sep=":")
          topmarkers
        }
      )
    }
    dat <- rbindlist(sdat)
    clusters <- unique(dat$clusterID)

    gene_list <- lapply(clusters, function(cid) {
      
      cluster_data <- dat %>% 
        dplyr::filter(clusterID == cid) %>%
        dplyr::arrange(desc(avg_log2FC)) %>%
        dplyr::slice_head(n = topn) %>%
        dplyr::pull(gene) %>%
        unique()
    })
    names(gene_list) <- clusters
    
    # 过滤marker数量少于3个的cluster
    delclus <- names(which(table(dat$clusterID)<3))
    if(length(delclus) > 0){
      gene_list <- gene_list[-which(clusters %in% delclus)]
    }
  
    return(gene_list)
}


#' @title Calculate the similarity score between two ordered gene sets
#'
#' @param set1 gene set 1
#' @param set2 gene set 2
#'
#' @return stat
#' @export
#'
sim_score <- function(set1, set2) {

  w1 <- 1 / seq_along(set1)
  w2 <- 1 / seq_along(set2)
  
  w_all <- c(setNames(w1, set1), setNames(w2, set2))
  union_set <- unique(c(set1, set2))
  intersection <- intersect(set1, set2)
  
  if(length(intersection)!=0){
    w_intersection <- sum(sapply(intersection, function(g) mean(w_all[g], na.rm = TRUE)))
    w_union <- sum(sapply(union_set, function(g) mean(w_all[g], na.rm = TRUE)))
    score <- w_intersection / w_union
  }else{
     score <- 0
  }

  return(score)
}


#' @title Calculate similarity matrix of two ordered gene list
#'
#' @param genelist1 Gene list 1
#' @param genelist2 Gene list 2
#' @param cores Number of CPU for running
#' 
#' @return Similarity matrix
#'
#' @export
get_sim_mat <- function(gene_list1, gene_list2 = NULL, cores){

    require(dplyr)
    require(foreach)
    require(data.table)
    require(parallel)
    require(doParallel)

    if(is.null(gene_list2)){
      gene_list2 <- gene_list1
    }

    clusters1 <- names(gene_list1)
    clusters2 <- names(gene_list2)

    # 开始计算相似矩阵
    cl <- makeCluster(cores)
    registerDoParallel(cl)

    sim_mat <- foreach(
      i = 1:length(clusters1),
      .combine = rbind,
      .export = c("sim_score")
    ) %dopar% {
      x <- gene_list1[[i]]
      o <- sapply(gene_list2, function(y) sim_score(x, y))
      o
    }
    stopCluster(cl)      

    rownames(sim_mat) <- clusters1
    colnames(sim_mat) <- clusters2

    return(sim_mat)
}


#' @title  Selecting an appropriate number of modules and defining the module identity for each cluster
#'
#' @param sim_mat A similarity matrix
#' @param acc_list_sta A list, where each element corresponds to cell types of a cluster
#' @param p_thres Proportion of modules meeting the m_thres threshold
#' @param m_thres The minimum purity threshold that each module satisfies
#' @param min_nm Minimum number of modules
#'
#' @return Number of modules and module identity for each cluster
#'
#' @export
get_module <- function(sim_mat, acc_list_sta, m_thres=0.8, p_thres=0.8, min_nm = NULL){ 

    min_nm <- max(round(dim(sim_mat)[1]/10,0),1,min_nm)
    max_nm <- round(dim(sim_mat)[1],0)
    if(min_nm>max_nm){
        min_nm <- max_nm
    }

    purity <- 0
    if (max_nm > 100) {
        step_size <- 100 
    } else {
        step_size <- 10 
    }

    best_nm <- min_nm
    best_purity <- 0

    # First stage (determining step size based on the size of min_nm and max_nm) 
    for (nm in c(seq(min_nm, max_nm, by = step_size),max_nm)) {
        clus <- hclust(as.dist(1 - sim_mat), method = "complete")
        out.id <- cutree(clus, nm)

        # Calculate cell type purity of each module
        acc_list <- try(module_celltype(acc_list_sta, out.id, nm)[[2]], silent = TRUE)
        if (class(acc_list) == "try-error") {
            acc_list = NA
        }

        purity1 <- length(which(acc_list > m_thres)) / length(acc_list)
        if (purity < p_thres && purity1 >= p_thres) {
            best_nm <- nm
            best_purity <- purity1
            search_steps <- seq(max(best_nm - step_size, min_nm), best_nm, by = step_size)
            step_size <- 10  
            break
        }
        if(nm==max_nm){
            best_nm <- max_nm
            best_purity <- purity1
            search_steps <- seq(max(max_nm - step_size, min_nm), max_nm, by = step_size)
            step_size <- 10
        }
        purity <- purity1
    }

    
    # Second stage (step size = 10)
    for (nm in search_steps) {
        clus <- hclust(as.dist(1 - sim_mat), method = "complete")
        out.id <- cutree(clus, nm)

        # Calculate cell type purity of each module
        acc_list <- try(module_celltype(acc_list_sta, out.id, nm)[[2]], silent = TRUE)
        if (class(acc_list) == "try-error") {
            acc_list = NA
        }

        purity1 <- length(which(acc_list > m_thres)) / length(acc_list)
        if (purity < p_thres && purity1 >= p_thres) {
            best_nm <- nm
            best_purity <- purity1
            search_steps <- seq(max(best_nm - step_size, min_nm), best_nm, by = step_size)
            step_size <- 1  # 满足条件时，细化步长为 1
            break
        }
        if(nm==max_nm){
            best_nm <- max_nm
            best_purity <- purity1
            search_steps <- seq(max(max_nm - step_size, min_nm), max_nm, by = step_size)
            step_size <- 1  
        }
        purity <- purity1
    }
    
    
    # Third stage (step size = 1)
    for (nm in search_steps) {
        clus <- hclust(as.dist(1 - sim_mat), method = "complete")
        out.id <- cutree(clus, nm)

        # Calculate cell type purity of each module
        acc_list <- try(module_celltype(acc_list_sta, out.id, nm)[[2]], silent = TRUE)
        if (class(acc_list) == "try-error") {
            acc_list = NA
        }

        purity1 <- length(which(acc_list > m_thres)) / length(acc_list)
        if (purity < p_thres && purity1 >= p_thres) {
            best_nm <- nm
            best_purity <- purity1
            break
        }
        if(nm==max_nm){
            best_nm <- max_nm
            best_purity <- purity1
        }
        purity <- purity1
    }
    
    clus <- hclust(as.dist(1-sim_mat), method = "complete")
    ids <- cutree(clus, best_nm)
    print(paste0("The recommended number of gene modules is ", best_nm))
    print(paste0("The best purity is ", round(best_purity,4)))
    return(list(ids=ids, best_nm=best_nm))
}


#' @title Calculate cell type composition of each cluster
#'
#' @param sim_mat A similarity matrix
#' @param clusters Clusters
#' @param types Cell types
#'
#' @return  List of composition of cell types for each cluster, and proportion of the largest cell type in the cluster
#'
#' @export
get_cluster_acc <- function(sim_mat, clusters, types){

    acc_list <- c()
    acc_list_sta <- list()
    for(j in 1:dim(sim_mat)[1]){
        all_celltype <- types[which(clusters %in% colnames(sim_mat)[j])]
        sta_celltype <- table(all_celltype)
        acc_list_sta[[j]] <- sta_celltype
        celltype <- sort(sta_celltype , decreasing = T)[1]
        celltype_acc <- celltype/length(all_celltype)
        acc_list <- c(acc_list, celltype_acc)
        #print(j)
    }
    acc_list_sta[[length(acc_list_sta)+1]] <- acc_list
    names(acc_list_sta) <- colnames(sim_mat)
    names(acc_list_sta)[length(acc_list_sta)] <- "acc_list"
    return(acc_list_sta)
}


#' @title alculate cell type composition of each cluster for all layers
#'
#' @param sim_mat A similarity matrix
#' @param info A information dataframe from Suerat object list
#' @param layer Layer
#'
#' @return List of composition of cell types for each cluster, and proportion of the largest cell type in the cluster for all layers
#'
#' @export
get_cluster_ratio <- function(sim_mat, info, layer){
  
    require(dplyr)

    info_celltype <- info %>% select(starts_with("Celltype"))

    if (!(class(layer) %in% c('NULL', 'numeric'))){
        stop('Please provide a numeric number or NULL as layer')
    } else if (class(layer) == 'NULL'){
      layer = 0
    } else {
       layer = layer
    }

    types <- info$Annotation
    clusters <- info$Cluster
    acc_list_sta <- list(get_cluster_acc(sim_mat, clusters, types))

    if (layer>0){
      for (lay in 1:layer){
        types <- info_celltype[[lay]]
        clusters <- info$Cluster
        acc_list_sta <- c(acc_list_sta, list(get_cluster_acc(sim_mat, clusters, types)))
      }
    }
    names(acc_list_sta) <- paste("layer", 0:layer, sep="")

    return(acc_list_sta)
}


#' @title Calculate cell type composition of each module
#'
#' @param acc_list_sta Cell type composition of each cluster
#' @param out.id Identity of the cluster to which each cell belongs
#' @param k Number of clusters
#'
#' @return List of composition of cell types for each module, and proportion of the largest cell type in the module
#'
#' @export
module_celltype <- function(acc_list_sta, out.id, k){
    acc_list1<-c()
    acc_list_sta1 <- list()
    for(target_clus in 1:k){
        cts <- c()
        for(j in which(out.id==target_clus)){
            cts <- c(cts, acc_list_sta[[names(out.id)[j]]])
        }
        if(length(cts)!=0){
            cts_sum <- as.vector(aggregate(as.numeric(cts), by=list(type=factor(names(cts))),sum))
            sta_celltype <- cts_sum$x
            names(sta_celltype) <- cts_sum$type
            acc_list_sta1[[target_clus]] <- sta_celltype
            celltype <- sort(sta_celltype , decreasing = T)[1]
            celltype_acc <- celltype/sum(cts_sum$x)
            acc_list1<-c(acc_list1, celltype_acc)
        }else{
            acc_list1<-c(acc_list1, NA)
        }

    }
    names(acc_list_sta1) <- paste("m",1:length(acc_list_sta1),sep="")
    return(list(acc_list_sta1, acc_list1))
}




#' @title Get a information dataframe from Suerat object list
#'
#' @param seulist A list of Seurat objects
#' @param type The type of Seurat: 'query' or 'reference'
#'
#' @return A information dataframe of Suerat object list
#' @export
#'
get_info <- function(seulist, type, clustername){
  require(dplyr)
  info <- data.frame()
  for (rds_num in 1:length(seulist)){
    rds_info <- seulist[[rds_num]]@meta.data %>% select(starts_with('scPlantGM'))
    rds_info$scPlantGM.cellname <- rownames(rds_info)
    # rds_info$scPlantGM.cellname <- paste(rds_info$scPlantGM.sample, rownames(rds_info),sep=':')
    rds_info$scPlantGM.cluster <- paste(rds_info$scPlantGM.sample, unlist(seulist[[rds_num]][[clustername]]), sep=':')
    rds_info <- rds_info %>% select(starts_with('scPlantGM'))
    info <- rbind(info,rds_info)
  }

  if (type == 'query'){
    colnames(info) <- c('Sample', 'Cell', 'Cluster')
  } else {
    info <- info[,c("scPlantGM.sample","scPlantGM.refanno","scPlantGM.cellname","scPlantGM.cluster")]
    colnames(info) <- c('Sample', 'Annotation', 'Cell', 'Cluster')
  }

  return(info)

}


#' @title Get layer information
#'
#' @param info_reference A dataframe
#' @param info_layer A dataframe of layer information
#'
#' @return A dataframe includes layer information
#' @export
#'
get_layers <- function(info_reference, info_layer){
    require(dplyr)
    cellanno <- info_reference[,c('Cell','Annotation')]
    spare_cellanno <- setdiff(unique(cellanno$Annotation),unique(unlist(info_layer)))
    if (length(spare_cellanno)!=0){
      warning(paste('These cell types can not be found in layer information:', paste(spare_cellanno, collapse=', ')), sep=' ')
    }
    layer_form <- data.frame()
    slash_form <- as.data.frame(matrix(rep('/', dim(info_reference)[1] * dim(info_layer)[2]), nrow = dim(info_reference)[1]))
    for (layer_num in dim(info_layer)[2]:1){
      layer_strc <- info_layer  %>% filter(info_layer[[layer_num]]!='/')
      layer_strc <- as.data.frame(unique(layer_strc[,1:layer_num]))
      cellanno_layer <- cellanno  %>% filter(Annotation %in% layer_strc[[layer_num]])
      layer_strc$Annotation <- layer_strc[[layer_num]]
      cellanno_layer <- left_join(cellanno_layer,layer_strc,by = 'Annotation')
      if(nrow(cellanno_layer)>0){
        if (layer_num!=dim(info_layer)[2]){
          cellanno_layer <- cbind(cellanno_layer, as.data.frame(slash_form[1:dim(cellanno_layer)[1],dim(slash_form)[2]:(layer_num+1)]))
        }
        colnames(cellanno_layer) <- c('Cell','Annotation',paste0(c('Celltype'),1:dim(info_layer)[2],sep=''))
        layer_form <- rbind(layer_form,cellanno_layer)
      }
    }
    layer_result <- left_join(info_reference,layer_form,by=c("Cell",'Annotation'))
    layer_result[is.na(layer_result)] <- '/'

    return(layer_result)
}


#' @title Calculate weight of genes afer merging clusters into module
#'
#' @param gene_sets List of ordered genes of all clusters
#' @param a A hyperparameter that adjusts the weights of gene order and frequency.
#' @param score TRUE outputs the score of each gene after merging clusters into module. By default FALSE. 
#'
#' @return Ordered gene vector
#' @export
#'
merged_weighted <- function(gene_sets, a = 0.5, score = FALSE) {

    all_genes <- unique(unlist(gene_sets))  
    gene_ranks <- matrix(NA, nrow = length(all_genes), ncol = length(gene_sets))
    gene_freqs <- numeric(length(all_genes)) 

    # 为每个基因集合中的基因分配排名
    for (i in seq_along(gene_sets)) {
        set <- gene_sets[[i]]
        ranks <- match(all_genes, set) 
        gene_ranks[, i] <- ranks
        gene_freqs[!is.na(ranks)] <- gene_freqs[!is.na(ranks)] + 1  # 统计基因频次
    }

    rank_df <- data.frame(gene = all_genes,
                        mean_rank = apply(gene_ranks, 1, function(x) {sum(x, na.rm=T)}),
                        freq = gene_freqs)
    rank_df$mean_rank <- rank_df$mean_rank / rank_df$freq


    # 排名越小越好 → normalize后要变成越大越好
    rank_df$rank_norm <- 1 - (rank_df$mean_rank - min(rank_df$mean_rank, na.rm = TRUE)) /
                            (max(rank_df$mean_rank, na.rm = TRUE) - min(rank_df$mean_rank, na.rm = TRUE))

    # 频次越大越好 → normalize到 [0,1]
    rank_df$freq_norm <- (rank_df$freq - min(rank_df$freq)) /
                        (max(rank_df$freq) - min(rank_df$freq))


    rank_df$score <- a * rank_df$rank_norm + (1 - a) * rank_df$freq_norm


    rank_df1 <- rank_df[order(-rank_df$score), ]
    gene_ranks1 <- gene_ranks[order(-rank_df$score), ]

    rownames(rank_df1) <- NULL
    rownames(gene_ranks1) <- NULL

    if(score){
        return(list(order = rank_df1$gene, scores = rank_df1$score))
    }else{
        return(list(order = rank_df1$gene))
    }
    
}
