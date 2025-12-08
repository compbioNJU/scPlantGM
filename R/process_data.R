#' @title Pre-process data and add corresponding order information to Seurat object
#'
#' @param seulist A list of Seurat objects or a Seurat object
#' @param seuorder A vector that show the order of Suerat in the list
#' @param type The type of Seurat: 'query' or 'reference'
#'
#' @return A list of Suerat object that has been well processed
#' @export
#'
#' @examples query <- process_obj(query, seuorder, 'query')

process_obj <- function(seulist, seuorder=NULL, type='query'){

    require(Seurat)
    
    if(class(seulist) == 'Seurat'){
        seulist <- list(seulist)
    } else if(class(seulist) == 'list'){
        seulist <- seulist
    } else {
       stop('Please provide a list of all query or reference Seurat objectives!')
    }

    if (is.null(seuorder)){
        seuorder <- paste0("Sample", length(seulist))
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
        print("Please remember to add celltypes of cells to the 'scPlantGM.refanno' column of meta data for reference Seurat objects!")
    }

    names(new_seulist) <- seuorder
    return(new_seulist)
}
