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
    } else if(any(class(x.slot) %in% c("dgCMatrix", "dsCMatrix", "matrix"))){
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
.group_nhoods_by_overlap <- function(nhs, da.res, is.da, overlap=1,
                                     lfc.threshold=NULL, merge.discord=FALSE,
                                     subset.nhoods=NULL, bits=FALSE, cells=NULL){

    ## Build adjacency matrix for nhoods
    if(isFALSE(bits)){
        nhood.adj <- .build_nhood_adjacency(nhs)
    }else if(isTRUE(bits)){
        nhood.adj <- .build_nhood_adjacency_bits(cells=cells, nhoods=nhs, overlap=1)
    }

    groups <- .group_nhoods_from_adjacency(nhs=nhs, nhood.adj=nhood.adj,
                                           is.da=is.da, da.res=da.res,
                                           subset.nhoods=subset.nhoods,
                                           overlap=overlap,
                                           lfc.threshold=lfc.threshold,
                                           merge.discord=merge.discord)

    return(groups)
}


#' @importFrom igraph graph_from_adjacency_matrix components
.group_nhoods_from_adjacency <- function(nhs, nhood.adj, da.res, is.da,
                                         merge.discord=FALSE,
                                         lfc.threshold=NULL,
                                         overlap=1, subset.nhoods=NULL){

    if(is.null(names(nhs))){
        warning("No names attributed to nhoods. Converting indices to names")
        names(nhs) <- as.character(c(1:length(nhs)))
    }

    # assume order of nhs is the same as nhood.adj
    if(!is.null(subset.nhoods)){
        if(mode(subset.nhoods) %in% c("character", "logical", "numeric")){
            # force use of logicals for consistency
            if(mode(subset.nhoods) %in% c("character", "numeric")){
                sub.log <- names(nhs) %in% subset.nhoods
            } else{
                sub.log <- subset.nhoods
            }

            nhood.adj <- nhood.adj[sub.log, sub.log]

            if(length(is.da) == length(nhs)){
                nhs <- nhs[sub.log]
                is.da <- is.da[sub.log]
                da.res <- da.res[sub.log, ]
            } else{
                stop("Subsetting `is.da` vector length does not equal nhoods length")
            }
        } else{
            stop(paste0("Incorrect subsetting vector provided:", class(subset.nhoods)))
        }
    } else{
        if(length(is.da) != ncol(nhood.adj)){
            stop("Subsetting `is.da` vector length is not the same dimension as adjacency")
        }
    }

    ## check for concordant signs - assume order is the same as nhoods
    if(isFALSE(merge.discord)){
        nonz.nhs <- colSums(nhood.adj) > 0
        ll_names <- expand.grid(c(1:length(nhs[nonz.nhs])), c(1:length(nhs[nonz.nhs])))

        concord.sign <- sapply(1:nrow(ll_names), function(x) sign(da.res[nonz.nhs, ][as.numeric(ll_names[x, 1]), ]$logFC) ==
                                   sign(da.res[nonz.nhs,][as.numeric(ll_names[x, 2]), ]$logFC))
        pairs_int <- sapply(which(concord.sign), function(x) length(intersect(nhs[nonz.nhs][[ll_names[x, 1]]], nhs[nonz.nhs][[ll_names[x, 2]]])))
        lintersect <- rep(0, nrow(ll_names))
        lintersect[concord.sign] <- pairs_int

        ## Count as connected only nhoods with at least n shared cells
        lintersect_filt <- ifelse(!concord.sign, 0, lintersect)
        ll_names <- cbind(ll_names, lintersect_filt)

        nhood.adj[nonz.nhs, nonz.nhs] <- ll_names[, 3]
    }

    if(overlap > 1){
        # loop over adj dimensions and mask out cells with insufficient overlapping cells
        nonz.nhs <- colSums(nhood.adj) > 0
        ll_names <- expand.grid(c(1:length(nhs[nonz.nhs])), c(1:length(nhs[nonz.nhs])))

        keep_pairs <- sapply(1:nrow(ll_names) , function(x) any(nhs[nonz.nhs][[ll_names[x, 1]]] %in% nhs[nonz.nhs][[ll_names[x, 2]]]))
        pairs_int <- sapply(which(keep_pairs), function(x) length(intersect(nhs[nonz.nhs][[ll_names[x, 1]]], nhs[nonz.nhs][[ll_names[x, 2]]])))
        lintersect <- rep(0, nrow(ll_names))
        lintersect[keep_pairs] <- pairs_int

        ## Count as connected only nhoods with at least n shared cells
        lintersect_filt <- ifelse(lintersect < overlap, 0, lintersect)
        ll_names <- cbind(ll_names, lintersect_filt)

        nhood.adj[nonz.nhs, nonz.nhs] <- ll_names[, 3]
    }

    if(!is.null(lfc.threshold)){
        nonz.nhs <- colSums(nhood.adj) > 0
        ll_names <- expand.grid(c(1:length(nhs[nonz.nhs])), c(1:length(nhs[nonz.nhs])))

        # set adjacency to 0 for nhoods with lfc < threshold
        lfc.pass <- sapply(1:nrow(ll_names), function(x) (abs(da.res[nonz.nhs, ][as.numeric(ll_names[x, 1]), ]$logFC) >= lfc.threshold) &
                               (abs(da.res[nonz.nhs, ][as.numeric(ll_names[x, 2]), ]$logFC) >= lfc.threshold))
        pairs_int <- sapply(which(lfc.pass), function(x) length(intersect(nhs[nonz.nhs][[ll_names[x, 1]]], nhs[nonz.nhs][[ll_names[x, 2]]])))
        lintersect[lfc.pass] <- pairs_int

        ## Count as connected only nhoods with at least n shared cells
        lintersect_filt <- ifelse(!lfc.pass, 0, lintersect)
        ll_names <- cbind(ll_names, lintersect_filt)

        nhood.adj[nonz.nhs, nonz.nhs] <- ll_names[, 3]
    }

    # binarise
    nhood.adj <- as.matrix((nhood.adj > 0) + 0)

    n.dim <- ncol(nhood.adj)
    if(!isSymmetric(nhood.adj)){
        stop("Overlap matrix is not symmetric")
    }

    if(nrow(nhood.adj) != ncol(nhood.adj)){
        stop("Non-square distance matrix - check nhood subsetting")
    }

    g <- graph_from_adjacency_matrix(nhood.adj, mode="undirected", diag=FALSE)
    groups <- components(g)$membership

    # only keep the groups that contain >= 1 DA neighbourhoods
    keep.groups <- intersect(unique(groups[is.da]), unique(groups))

    return(groups[groups %in% keep.groups])
}


#### nhood adjacency matrix function
# Build adjacency matrix of overlap between neighbourhoods
#' @importFrom gtools permutations
.build_nhood_adjacency <- function(nhoods, overlap=1){
    nms <- permutations(n = length(nhoods), v = names(nhoods), r = 2, repeats.allowed = TRUE)
    # keep_pairs <- sapply( 1:nrow(nms) , function(x) any(nhoods[[nms[x,1]]] %in% nhoods[[ nms[x,2] ]]))
    keep_pairs <- sapply( 1:nrow(nms), function(x) sum(nhoods[[nms[x,1]]] %in% nhoods[[ nms[x,2] ]]) > overlap)

    message("Calculating nhood adjacency")
    pairs_int <- sapply( which(keep_pairs), function(x) length( intersect( nhoods[[nms[x,1]]], nhoods[[ nms[x,2] ]]) ) )
    out <- rep(0, nrow(nms))
    out[which(keep_pairs)] <- pairs_int

    nh_intersect_mat <- matrix(out, nrow = length(nhoods), byrow = TRUE)
    rownames(nh_intersect_mat) <- unique(nms[,1])
    colnames(nh_intersect_mat) <- unique(nms[,2])
    return(nh_intersect_mat)
}


#' @importFrom bit as.bit
.build_nhood_adjacency_bits <- function(cells, nhoods, overlap=1){
    # each neighbourhood is a bit-vector with length(cells)
    # cell sharing is then determined using bit-wise operations
    bit.list <- sapply(1:length(nhoods), FUN=function(B) cells %in% c(as.numeric(nhoods[[B]]), names(nhoods)[B]))
    bit.adj.list <- lapply(bit.list, FUN=function(BX) unlist(lapply(bit.list, FUN=function(BK) sum(BX & BK))))
    bit.mat <- do.call(cbind, bit.adj.list)

    return(bit.mat)
}
