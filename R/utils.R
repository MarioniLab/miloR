## These are utility functions not meant to be exposed to the user

#' @importFrom methods slot
#' @importFrom Matrix rowSums
.check_empty <- function(x, attribute){
    # check if a Milo object slot is empty or not
    x.slot <- slot(x, attribute)

    if(is.list(x.slot) & names(slot(x, "graph")) == "graph"){
        return(length(x.slot[[1]]) > 0)
    } else if(is.list(x.slot) & is.null(names(x.slot))){
        return(length(x.slot))
    } else if(any(class(x.slot) %in% c("dgCMatrix", "dsCMatrix", "ddiMatrix", "matrix"))){
        return(sum(rowSums(x.slot)) == 0)
    }
}


.check_binary <- function(x){
    # check if a matrix is binary or not
    sum.zeros <- sum(x == 0)
    sum.ones <- sum(x == 1)
    n.comps <- nrow(x) * ncol(x)

    return(sum(c(sum.zeros, sum.ones)) == n.comps)
}


#' @importFrom igraph make_graph simplify
.neighborsToKNNGraph <- function(nn, directed=FALSE) {
    start <- as.vector(row(nn))
    end <- as.vector(nn)
    interleaved <- as.vector(rbind(start, end))

    if (directed) {
        g <- make_graph(interleaved, directed=TRUE)

    } else {
        g <- make_graph(interleaved, directed=FALSE)
        g <- simplify(g, edge.attr.comb = "first")
    }
    g
}

# setting internals for replacement methods that require multiple arguments - borrowed from SingleCellExperiment
#' @importFrom methods slot
.set_reduced_dims <- function(x, value, slot.x=NULL, rdim=NULL){
    x <- updateObject(x)
    content <- slot(x, slot.x)

    if(slot.x == "nhoodReducedDim"){

        if(!is.null(rdim)){
            content[[rdim]] <- value
            x@nhoodReducedDim <- content
        } else{
            stop("No reduced dimensionality slot provided")
        }
    }else{
        stop(paste0("replacement method not implemented for ", slot))
    }

    x
}

######################################
## neighbourhood grouping functions
######################################

#' @importFrom igraph graph_from_adjacency_matrix components
# .group_nhoods_by_overlap <- function(nhs, da.res, is.da, overlap=1,
#                                      lfc.threshold=NULL, merge.discord=FALSE,
#                                      subset.nhoods=NULL, cells=NULL){
# 
#     ## Build adjacency matrix for nhoods
#     nhood.adj <- .build_nhood_adjacency(nhs)
#     groups <- .group_nhoods_from_adjacency(nhs=nhs, nhood.adj=nhood.adj,
#                                            is.da=is.da, da.res=da.res,
#                                            subset.nhoods=subset.nhoods,
#                                            overlap=overlap,
#                                            lfc.threshold=lfc.threshold,
#                                            merge.discord=merge.discord)
# 
#     return(groups)
# }


#' @importFrom igraph graph_from_adjacency_matrix components cluster_louvain
# .group_nhoods_from_adjacency <- function(nhs, nhood.adj, da.res, is.da,
#                                          merge.discord=FALSE,
#                                          lfc.threshold=NULL,
#                                          overlap=1, subset.nhoods=NULL){
# 
#     if(is.null(colnames(nhs))){
#         warning("No names attributed to nhoods. Converting indices to names")
#         colnames(nhs) <- as.character(c(1:ncol(nhs)))
#     }
# 
#     # assume order of nhs is the same as nhood.adj
#     if(!is.null(subset.nhoods)){
#         if(mode(subset.nhoods) %in% c("character", "logical", "numeric")){
#             # force use of logicals for consistency
#             if(mode(subset.nhoods) %in% c("character", "numeric")){
#                 sub.log <- colnames(nhs) %in% subset.nhoods
#             } else{
#                 sub.log <- subset.nhoods
#             }
# 
#             nhood.adj <- nhood.adj[sub.log, sub.log]
# 
#             if(length(is.da) == ncol(nhs)){
#                 nhs <- nhs[sub.log]
#                 is.da <- is.da[sub.log]
#                 da.res <- da.res[sub.log, ]
#             } else{
#                 stop("Subsetting `is.da` vector length does not equal nhoods length")
#             }
#         } else{
#             stop(paste0("Incorrect subsetting vector provided:", class(subset.nhoods)))
#         }
#     } else{
#         if(length(is.da) != ncol(nhood.adj)){
#             stop("Subsetting `is.da` vector length is not the same dimension as adjacency")
#         }
#     }
# 
#     ## check for concordant signs - assume order is the same as nhoods
#     if(isFALSE(merge.discord)){
#         # nonz.nhs <- colSums(nhood.adj) > 0
#         discord.sign <- sign(da.res[, 'logFC'] %*% t(da.res[, 'logFC'])) < 0
#         nhood.adj[discord.sign] <- 0
#     }
# 
#     if(overlap > 1){
#         nhood.adj[nhood.adj < overlap] <- 0
#     }
# 
#     if(!is.null(lfc.threshold)){
#         nhood.adj[,which(da.res$logFC < lfc.threshold)] <- 0
#         nhood.adj[which(da.res$logFC < lfc.threshold),] <- 0
#     }
# 
#     # binarise
#     nhood.adj <- as.matrix((nhood.adj > 0) + 0)
# 
#     n.dim <- ncol(nhood.adj)
#     if(!isSymmetric(nhood.adj)){
#         stop("Overlap matrix is not symmetric")
#     }
# 
#     if(nrow(nhood.adj) != ncol(nhood.adj)){
#         stop("Non-square distance matrix - check nhood subsetting")
#     }
# 
#     g <- graph_from_adjacency_matrix(nhood.adj, mode="undirected", diag=FALSE)
#     groups <- cluster_louvain(g)$membership
#     names(groups) <- colnames(nhood.adj)
#     # groups <- components(g)$membership
# 
#     # only keep the groups that contain >= 1 DA neighbourhoods
#     keep.groups <- intersect(unique(groups[is.da]), unique(groups))
# 
#     return(groups[groups %in% keep.groups])
# }


#### nhood adjacency matrix function
# Build adjacency matrix of overlap between neighbourhoods
#' @importFrom gtools permutations
#' @importFrom Matrix crossprod
.build_nhood_adjacency <- function(nhoods, overlap=1){
    nh_intersect_mat <- Matrix::crossprod(nhoods)
    nh_intersect_mat[nh_intersect_mat < overlap] <- 0

    rownames(nh_intersect_mat) <- colnames(nhoods)
    colnames(nh_intersect_mat) <- colnames(nhoods)
    return(nh_intersect_mat)
}

