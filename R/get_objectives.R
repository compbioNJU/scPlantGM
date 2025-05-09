#' @title get markers from Seurat object
#'
#' @param rds a Seurat object
#'
#' @return markers
#'
#' @export
get_markers <- function(rds, cores){
    require(dplyr)
    require(doMC)
    seurat_obj <- rds

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    #seurat_obj$seurat_clusters <- seurat_obj$integrated_snn_res.0.8
    Idents(seurat_obj) <- "seurat_clusters"
    csx <- table(seurat_obj$seurat_clusters)
    if(all(csx > 5)){
        top.markers <- foreach(i=sort(unique(seurat_obj$seurat_clusters)),.combine=rbind, .packages = "Seurat") %dopar% {
            if(csx[as.character(i)]>5){
            o <- FindMarkers(seurat_obj, ident.1=i, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1,verbose = FALSE)
            data.frame(o, gene=rownames(o), cluster=i)
            }else{
            NULL
            }
        }
    }else{
        top.markers <- FindAllMarkers(seurat_obj, min.pct=0.1, logfc.threshold=0.1,
                                        return.thresh=0.1, only.pos=TRUE, verbose = FALSE) ## return.thresh=0.01, test.use="wilcox"
    }
    stopCluster(cl)

    top.markers$pct.diff <- top.markers$pct.1 - top.markers$pct.2
    top.markers$name <- seurat_obj@misc$geneName[top.markers$gene]

    ## only keep top 100
    topmarkers <- top.markers[top.markers$p_val_adj < 0.05,] %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n=100, wt=avg_log2FC)
    top20 <- topmarkers %>% dplyr::group_by(cluster) %>%
        dplyr::top_n(n=20, wt=avg_log2FC)

    stat <- table(seurat_obj$seurat_clusters) %>% as.data.frame() %>%
            rename(cluster=Var1, number=Freq)
    rownames(stat) <- stat$cluster

    seurat_obj@misc$cellstat <- stat
    seurat_obj@misc$markerGenes <- top.markers
    seurat_obj@misc$topMarker <- topmarkers
    seurat_obj@misc$top20Marker <- top20

    return(seurat_obj@misc)
}

#' @title get top markers
#'
#' @param misc misc
#'
#' @return top markers
#'
#' @export
get_topmks <- function(misc){
  topmks <- misc$topMarker %>% 
    filter(p_val_adj < 0.05 & pct.1>0.2 & avg_log2FC > 0.25 & !grepl("^MT-", name)) %>%
    dplyr::group_by(cluster) %>% 
    dplyr::top_n(n=100, wt=avg_log2FC)

    return(topmks)
}

#' @title get dat
#'
#' @param seulist a list of seurat objects
#'
#' @return dat
#' @export
#'
get_dat <- function(seulist, cores){
  misc <- get_markers(seulist, cores)
  topmks <- get_topmks(misc)
  assay_sample <- seulist@meta.data$scPlantGM.sample[1:nrow(topmks)]
  topmks_df <- data.table(topmks, assay=assay_sample)

  return(topmks_df)
}

#' @title get cells
#'
#' @param seulist a list of seurat objects
#'
#' @return cells
#' @export
#'
get_cells <- function(seulist){
  stat <- get_cellstat(seulist)
  assay_sample <- seulist@meta.data$scPlantGM.sample[1:nrow(stat)]
  cells_df <- data.table(stat, assay=assay_sample)
  
  return(cells_df)
}


#' @title to get statistic of cells
#'
#' @param rds seurat object
#'
#' @return stat
#' @export
#'
get_cellstat <- function(rds){
    require(dplyr)
    stat <- table(rds$seurat_clusters) %>% as.data.frame() %>% 
      rename(cluster=Var1, number=Freq)
    rownames(stat) <- stat$cluster
    return(stat)
}

#' @title calculate jaccard matrix of Seurat object
#'
#' @param seulist a list of Seurat object
#' @param type query or reference
#' @param cores cores
#' 
#' @return jaccard matrix
#'
#' @export
get_jaccardmat <- function(seulist,type,cores){

    require(Seurat)
    require(dplyr)
    require(foreach)
    require(data.table)
    require(parallel)
    require(doParallel)


    sdat <- lapply(seulist, function(x) get_dat(x, cores))
    scells <- lapply(seulist, get_cells) 

    dat <- rbindlist(sdat)
    cells <- rbindlist(scells)

    cellx <- cells %>% group_by(assay) %>% summarize(n=sum(number))
    cellx <- setNames(cellx$n, cellx$assay)
    cells <- cells %>% 
      mutate(proportion=number/cellx[assay]) %>% 
      mutate(clusterID=sprintf("%s:%s", assay, cluster))
    meta <- dat %>% select(assay, cluster) %>% unique()
    rownames(meta) <- sprintf("%s:%s", meta$assay, meta$cluster) 


    if (type=='query'){
      out = c()
      return(list(out=out, dat=dat, cells=cells, meta=meta))
    } else {
      meta_assays <- meta$assay  
      meta_clusters <- meta$cluster  

      meta_list <- lapply(1:nrow(meta), function(i) {  
        dat %>% filter(assay == meta_assays[i] & cluster == meta_clusters[i]) %>% pull(gene)  
      })

      cl <- makeCluster(cores)
      registerDoParallel(cl)

      out <- foreach(i = 1:nrow(meta), .combine = 'rbind') %dopar% {
        x <- meta_list[[i]]
        o <- sapply(1:nrow(meta), function(j) {
          y <- meta_list[[j]]
          z <- length(union(x, y)) 
          if(z > 0) {
            length(intersect(x, y)) / z
          } else {
            0
          }
        })
        o
      }
      stopCluster(cl)      

      rownames(out) <- rownames(meta)
      colnames(out) <- rownames(meta)

      meta <- meta %>% ## mutate(accession=gsub("_.*","",assay)) %>%
        mutate(clusterN=table(meta$assay)[assay]) %>%
        as.data.frame()
      rownames(meta) <- sprintf("%s:%s", meta$assay, meta$cluster) 
      dat <- dat %>% mutate(clusterID=sprintf("%s:%s", assay, cluster)) %>%
        mutate(markerN=table(clusterID)[clusterID]) %>%
        as.data.frame()
      
      # 过滤marker数量少于3个的cluster
      delclus <- names(which(table(dat$clusterID)<3))
      if(length(delclus)>0){
        out <- out[-which(rownames(out) %in% delclus),-which(rownames(out) %in% delclus)]
        dat <- dat[-which(dat$clusterID %in% delclus),]
        cells <- cells[-which(cells$clusterID %in% delclus),]
        meta <- meta[-which(rownames(meta) %in% delclus),]
      }
    
      return(list(out=out,dat=dat,cells=cells,meta=meta))
    }
}

#' @title to get a whole reference combined input and built-in reference
#'
#' @param jaccard.mat1 jaccard matrix
#' @param jaccard.mat2 jaccard matrix
#' @param info_reference data.frame
#' @param cores cores
#'
#' @return jaccard matrix
#' @export
#'
fuse_ref_jm <- function(jaccard.mat1, jaccard.mat2, info_reference, cores){

  dat1 <- jaccard.mat1[[2]]
  cells1 <- jaccard.mat1[[3]]

  dat2 <- jaccard.mat2[[2]]  %>% filter(clusterID %in% unique(info_reference$Cluster))
  cells2 <- jaccard.mat2[[3]] %>% filter(clusterID %in% unique(info_reference$Cluster))

  dat <- rbind(dat1,dat2)
  cells <- rbind(cells1,cells2)

  meta <- dat %>% select(assay, cluster) %>% unique()
  rownames(meta) <- sprintf("%s:%s", meta$assay, meta$cluster)

  meta_assays <- meta$assay  
  meta_clusters <- meta$cluster  

  meta_list <- lapply(1:nrow(meta), function(i) {  
    dat %>% filter(assay == meta_assays[i] & cluster == meta_clusters[i]) %>% pull(gene)  
  })

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  out <- foreach(i = 1:nrow(meta), .combine = 'rbind') %dopar% {
    x <- meta_list[[i]]
    o <- sapply(1:nrow(meta), function(j) {
      y <- meta_list[[j]]
      z <- length(union(x, y)) 
      if(z > 0) {
        length(intersect(x, y)) / z
      } else {
        0
      }
    })
    o
  }

  stopCluster(cl)
  rownames(out) <- rownames(meta)
  colnames(out) <- rownames(meta)

  meta <- meta %>% ## mutate(accession=gsub("_.*","",assay)) %>%
    mutate(clusterN=table(meta$assay)[assay])
  dat <- dat %>% mutate(clusterID=sprintf("%s:%s", assay, cluster)) %>%
    mutate(markerN=table(clusterID)[clusterID])
    
  return(list(out,dat,cells,meta))
}



#' @title get a fusion matrix of two jaccard matrix
#'
#' @param jaccard.mat1 a jaccard matrix
#' @param jaccard.mat2 a jaccard matrix
#' @param cores cores
#'
#' @return a fusion jaccard matrix
#'
#' @export
fuse_jaccardmat <- function(jaccard.mat1, jaccard.mat2, cores){
    require(dplyr)

    meta1 <- jaccard.mat1[[4]]
    dat1 <- jaccard.mat1[[2]]

    meta2 <- jaccard.mat2[[4]]   
    dat2 <- jaccard.mat2[[2]]

    out <- matrix(0, nrow = nrow(meta1), ncol = nrow(meta2), 
                dimnames = list(rownames(meta1), rownames(meta2)))

    meta_assays1 <- meta1$assay
    meta_clusters1 <- meta1$cluster
    meta_list1 <- lapply(1:nrow(meta1), function(i) {  
      dat1 %>% filter(assay == meta_assays1[i] & cluster == meta_clusters1[i]) %>% pull(gene)  
    })  

    meta_assays2 <- meta2$assay
    meta_clusters2 <- meta2$cluster
    meta_list2 <- lapply(1:nrow(meta2), function(i) {  
      dat2 %>% filter(assay == meta_assays2[i] & cluster == meta_clusters2[i]) %>% pull(gene)  
    }) # 这里拖慢速度

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    out <- foreach(i = 1:nrow(meta1), .combine = 'rbind') %dopar% {
      x <- meta_list1[[i]]
      o <- sapply(1:nrow(meta2), function(j) {
        y <- meta_list2[[j]]
        z <- length(union(x, y)) 
        if(z > 0) {
          length(intersect(x, y)) / z
        } else {
          0
        }
      })
      o
    }

    stopCluster(cl)
    rownames(out) <- rownames(meta1)
    colnames(out) <- rownames(meta2)

    return(out)
}


#' @title  get modules
#'
#' @param jaccard_mat a jaccard matrix
#' @param acc_list_sta accuracy list
#' @param p_thres a threshold
#' @param m_thres a threshold
#'
#' @return module
#'
#' @export
get_module <- function(jaccard_mat,acc_list_sta,p_thres=0.8,m_thres=0.9){

    nm1 <- max(round(dim(jaccard_mat)[1]/15,0),1)
    nm2 <- round(dim(jaccard_mat)[1],0)

    purity <- 0
    for(nm in nm1:nm2){
        clus <- hclust(as.dist(1-jaccard_mat), method = "complete")
        out.id <- cutree(clus, nm)

      
        #查看每个module中cell type情况
        acc_list <- try(module_celltype(acc_list_sta, out.id, nm)[[2]], silent = TRUE)
        if(class(acc_list)=="try-error"){acc_list=NA} 
        purity1 <- length(which(acc_list>p_thres))/length(acc_list)
        if(purity<m_thres&purity1>=m_thres){break}
        purity <- purity1
        #print(purity1)
    }

    clus <- hclust(as.dist(1-jaccard_mat), method = "complete")
    out.id <- cutree(clus, nm)
    #print(nm)
    return(list(out.id,nm))
}

#' @title get cluster accuracy
#'
#' @param jaccard_mat a jaccard matrix
#' @param clusters clusters
#' @param types type
#'
#' @return cluster accuracy
#'
#' @export
get_cluster_acc <- function(jaccard_mat, clusters, types){

    acc_list <- c()
    acc_list_sta <- list()
    for(j in 1:dim(jaccard_mat)[1]){
        all_celltype <- types[which(clusters %in% colnames(jaccard_mat)[j])]
        sta_celltype <- table(all_celltype)
        acc_list_sta[[j]] <- sta_celltype
        celltype <- sort(sta_celltype , decreasing = T)[1]
        celltype_acc <- celltype/length(all_celltype)
        acc_list<-c(acc_list, celltype_acc)
        #print(j)
    }
    acc_list_sta[[length(acc_list_sta)+1]] <- acc_list
    names(acc_list_sta) <- colnames(jaccard_mat)
    names(acc_list_sta)[length(acc_list_sta)] <- "acc_list"
    return(acc_list_sta)
}

#' @title get accuracy list
#'
#' @param jaccard.mat a jaccard matrix
#' @param info a data.frame
#' @param layer layer
#'
#' @return accuracy list
#'
#' @export
get_cluster_ratio <- function(jaccard.mat, info, layer){
  
    require(dplyr)

    jaccard_mat <- jaccard.mat[[1]]
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
    acc_list_sta <- list(get_cluster_acc(jaccard_mat, clusters, types))

    if (layer>0){
      for (lay in 1:layer){
        types <- info_celltype[[lay]]
        clusters <- info$Cluster
        acc_list_sta <- c(acc_list_sta, list(get_cluster_acc(jaccard_mat, clusters, types)))
      }
    }

    return(acc_list_sta)
}

#' @title get connection between modules and cell types
#'
#' @param acc_list_sta accuracy list
#' @param out.id out.id
#' @param k k
#'
#' @return a list of accuracy list
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
    return(list(acc_list_sta1,acc_list1))
}

#' @title to get a information dataframe from Suerat object list
#'
#' @param seulist a list of seurat objects
#' @param type type
#'
#' @return a dataframe
#' @export
#'
get_info <- function(seulist, type){
  require(dplyr)
  info <- data.frame()
  for (rds_num in 1:length(seulist)){
    rds_info <- seulist[[rds_num]]@meta.data %>% select(starts_with('scPlantGM'))
    rds_info$scPlantGM.cellname <- paste(rds_info$scPlantGM.sample,rownames(rds_info),sep=':')
    rds_info$scPlantGM.cluster <- paste(rds_info$scPlantGM.sample,seulist[[rds_num]]@meta.data$seurat_clusters,sep=':')
    rds_info <- rds_info %>% select(starts_with('scPlantGM'))
    info <- rbind(info,rds_info)
  }

  if (type == 'query'){
    colnames(info) <- c('Sample', 'Cell','Cluster')
  } else {
    info <- info[,c("scPlantGM.sample","scPlantGM.refanno","scPlantGM.cellname","scPlantGM.cluster")]
    colnames(info) <- c('Sample', 'Annotation', 'Cell','Cluster')
  }

  return(info)

}

#' @title to get layer information
#'
#' @param info_reference a data.frame
#' @param info_layer a data.frame of layer information
#'
#' @return a data.frame includes layer information
#' @export
#'
get_layers <- function(info_reference,info_layer){
    require(dplyr)
    cellanno <- info_reference[,c('Cell','Annotation')]
    spare_cellanno <- setdiff(unique(cellanno$Annotation),unique(unlist(info_layer)))
    if (length(spare_cellanno)!=0){
      warning(paste('These cell types can not be found in layer information:', paste(spare_cellanno,collapse=', ')), sep=' ')
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
