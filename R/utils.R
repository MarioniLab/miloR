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
                                     subset.nhoods=NULL){
    if(is.null(names(nhs))){
        warning("No names attributed to nhoods. Converting indices to names")
        names(nhs) <- as.character(c(1:length(nhs)))
    }

    if(!is.null(subset.nhoods)){
        if(mode(subset.nhoods) %in% c("character", "logical", "numeric")){
            sub.vec <- c(1:length(nhs))[subset.nhoods]
            #nhs <- nhs[sub.vec]
            ll_names <- expand.grid(sub.vec, sub.vec)
            n.dim <- length(sub.vec)
            if(length(is.da) == length(names(nhs))){
                is.da <- is.da[subset.nhoods]
            } else{
                stop("Subsetting `is.da` vector length does not equal nhoods length")
            }
        } else{
            stop(paste0("Incorrect subsetting vector provided:", class(subset.nhoods)))
        }
    } else{
        if(length(is.da) != length(names(nhs))){
            stop("Subsetting `is.da` vector length does not equal nhoods length")
        }

        ll_names <- expand.grid(c(1:length(nhs)), c(1:length(nhs)))
        n.dim <- length(names(nhs))
    }

    keep_pairs <- sapply(1:nrow(ll_names) , function(x) any(nhs[[ll_names[x, 1]]] %in% nhs[[ll_names[x, 2]]]))
    pairs_int <- sapply(which(keep_pairs), function(x) length(intersect(nhs[[ll_names[x, 1]]], nhs[[ll_names[x, 2]]])))
    lintersect <- rep(0, nrow(ll_names))

    if(length(pairs_int) != sum(keep_pairs)){
        stop("Incorrect number of neighbourhood pairs selected")
    }

    lintersect[keep_pairs] <- pairs_int

    ## Count as connected only nhoods with at least n shared cells
    lintersect_filt <- ifelse(lintersect < overlap, 0, lintersect)

    if(!is.null(lfc.threshold)){
        # set adjacency to 0 for nhoods with lfc < threshold
        lfc.pass <- sapply(1:nrow(ll_names), function(x) (abs(da.res[as.numeric(ll_names[x, 1]), ]$logFC) >= lfc.threshold) &
                               (abs(da.res[as.numeric(ll_names[x, 2]), ]$logFC) >= lfc.threshold))
        lintersect_filt <- ifelse(lfc.pass, 0, lintersect_filt)
    }

    ## check for concordant signs - assume order is the same as nhoods
    if(isFALSE(merge.discord)){
        concord.sign <- sapply(which(keep_pairs), function(x) sign(da.res[as.numeric(ll_names[x, 1]), ]$logFC) !=
                                   sign(da.res[as.numeric(ll_names[x, 2]), ]$logFC))
        lintersect_filt <- lintersect
        lintersect_filt[which(keep_pairs)] <- ifelse(concord.sign, 0, lintersect_filt[which(keep_pairs)])
    }

    ## Convert to adjacency matrix (values = no of common cells)
    ## This should be binary for an adjacency matrix
    d <- matrix(as.numeric(lintersect_filt > 0), nrow = n.dim, byrow = TRUE)

    if(!isSymmetric(d)){
        stop("Overlap matrix is not symmetric")
    }

    if(nrow(d) != ncol(d)){
        stop("Non-square distance matrix - check nhood subsetting")
    }

    g <- graph_from_adjacency_matrix(d, mode="undirected", diag=FALSE)
    groups <- components(g)$membership

    # only keep the groups that contain >= 1 DA neighbourhoods
    if(!is.null(subset.nhoods)){
        names(groups) <- names(nhs[subset.nhoods])
    } else{
        names(groups) <- names(nhs)
    }

    #keep.groups <- intersect(unique(groups[is.da]), unique(groups))
    # what about return ALL groups, and subset after?
    return(groups)#[groups %in% keep.groups])
}

#' @importFrom igraph graph_from_adjacency_matrix components
.group_nhoods_from_adjacency <- function(nhs, nhood.adj, da.res, is.da,
                                         merge.discord=FALSE,
                                         overlap=1, subset.nhoods=NULL){
    print(dim(nhood.adj))
    print(length(nhs))
    if(is.null(names(nhs))){
        warning("No names attributed to nhoods. Converting indices to names")
        names(nhs) <- as.character(c(1:length(nhs)))
    }

    # assume order of nhs is the same as nhood.adj
    if(!is.null(subset.nhoods)){
        if(mode(subset.nhoods) %in% c("character", "logical", "numeric")){
            nhood.adj <- nhood.adj[subset.nhoods]

            if(length(is.da) == length(nhs)){
                is.da <- is.da[subset.nhoods]
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
    if(!is.null(subset.nhoods)){
        names(groups) <- names(nhs[subset.nhoods])
    } else{
        names(groups) <- names(nhs)
    }

    keep.groups <- intersect(unique(groups[is.da]), unique(groups))

    return(groups[groups %in% keep.groups])
}




