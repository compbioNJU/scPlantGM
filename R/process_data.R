#' @title to pre-process Suerat objects into a fixed format
#'
#' @param seulist one/a list of Seurat objects
#' @param seuorder a vector that show the order of Suerat in the list. Strongly recommend using 'Sample' to present the order
#' @param type the type of Suerat: 'query' or 'reference'
#'
#' @return a list of Suerat that has been well processed
#' @export
#'
#' @examples query <- process_obj(query, seuorder, 'query')

process_obj <- function(seulist, seuorder, type='query'){

    require(Seurat)
    
    if(class(seulist) == 'Seurat'){
        seulist <- list(seulist)
    } else if(class(seulist) == 'list'){
        seulist <- seulist
    } else {
       stop('Please provide a list of all query or reference Seurat objectives!')
    }

    if (class(seuorder) == 'NULL'){
        stop('Please do not input an empty order factor!')
    }

    if (length(seulist)!=length(seuorder)){
        stop('The length of order vector should match the numbers of Seurat objects in seulist!')
    }

    new_seulist <- c()
    for (num in 1:length(seulist)){
        rds <- seulist[[num]]
        rds@meta.data$scPlantGM.sample <- rep(as.character(seuorder[num]),times=dim(rds)[2])
        new_seulist <- c(new_seulist,rds)

        if (!("seurat_clusters" %in% colnames(rds@meta.data))){
            print(paste('The meta.data of NO.',num,' Seurat object should include cluster info in a column named "seurat_clusters".',sep = ''))
        }
    }

    if (type=='reference'){
        print("Please remember to input celltypes of cells to the '@meta.data$scPlantGM.refanno' column of reference Seurat objects!")
    }

    return(new_seulist)
}
