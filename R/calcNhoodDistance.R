#' Calculate within neighbourhood distances
#'
#' This function will calculate Euclidean distances between single-cells in a
#' neighbourhood using the same dimensionality as was used to construct the graph.
#' This step follows the \code{makeNhoods} call to limit the number of distance
#' calculations required.
#'
#' @param x A \code{\linkS4class{Milo}} object with a valid \code{graph} slot.
#' If \code{reduced.dims} is not provided and there is no valid populated \code{reducedDim}
#' slot in \code{x}, then this is computed first with \code{d + 1} principal components.
#' @param d The number of dimensions to use for computing within-neighbourhood
#' distances. This should be the same value used construct the \code{graph}.
#' @param reduced.dim If x is an \code{\linkS4class{Milo}} object, a character indicating the name of the \code{reducedDim} slot in the
#' \code{\linkS4class{Milo}} object to use as (default: 'PCA'). Otherwise this should be an N X P matrix with rows in the same order as the
#' columns of the input Milo object \code{x}.
#'
#' @return A \code{\linkS4class{Milo}} object with the distance slots populated.
#'
#' @author
#' Mike Morgan
#'
#' #' @examples
#' library(SingleCellExperiment)
#' ux <- matrix(rpois(12000, 5), ncol=200)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#'
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'                             reducedDims=SimpleList(PCA=pca$x))
#'
#' milo <- Milo(sce)
#' milo <- buildGraph(milo, d=30, transposed=TRUE)
#' milo <- calcNhoodDistance(milo, d=30)
#'
#' milo
#' @name calcNhoodDistance
NULL

#' @export
#' @rdname calcNhoodDistance
#' @importFrom irlba prcomp_irlba
#' @importFrom SummarizedExperiment assay
#' @importFrom reshape2 melt
calcNhoodDistance <- function(x, d, reduced.dim=NULL, use.assay="logcounts", ...){
    if(class(x) == "Milo"){
        # check for reducedDims
        if(is.null(reducedDim(x)) & is.null(reduced.dim)){
            # assume logcounts is present?
            message("Computing PCA on input")
            x_pca <- prcomp_irlba(t(assay(x, use.assay)), n=min(d+1, ncol(x)-1),
                                  scale.=TRUE, center=TRUE)
            reducedDim(x, "PCA") <- x_pca$x
        } else if(is.null(reducedDim(x)) & is.character(reduced.dim)){
            stop(paste(reduce.dim, " not found in reducedDim slot"))
        }
    } else{
        stop("Input is not a valid Milo object")
    }

    if(any(names(reducedDims(x)) %in% c("PCA"))){
        nhood.dists <- sapply(names(nhoods(x)),
                              function(X) .calc_distance(reducedDim(x, "PCA")[c(as.numeric(X), unlist(nhoods(x)[[X]])), c(1:d),drop=FALSE]))
        names(nhood.dists) <- names(nhoods(x))

    } else if(is.character(reduced.dim)){
        nhood.dists <- sapply(names(nhoods(x)),
                              function(X) .calc_distance(reducedDim(x, reduced.dim)[c(as.numeric(X), unlist(nhoods(x)[[X]])), c(1:d),drop=FALSE]))
        names(nhood.dists) <- names(nhoods(x))
    }

    # ensure garbage collection
    sink(file="/dev/null")
    gc()
    sink(file=NULL)
    all.dist <- .aggregate_dists(nhoods(x), nhood.dists, colnames(x))
    nhoodDistances(x) <- all.dist

    return(x)
}


#' @importFrom Matrix rowSums sparseMatrix
#' @importFrom methods as
#' @export
.calc_distance <- function(in.x){

    dist.list <- list()
    for(i in seq_along(1:nrow(in.x))){
        i.diff <- t(apply(in.x, 1, FUN=function(P) P - in.x[i, ]))
        i.dist <- sqrt(rowSums(i.diff**2))
        dist.list[[paste0(i)]] <- list("rowIndex"=rep(i, nrow(in.x)), "colIndex"=c(1:length(i.dist)),
                                       "dist"=i.dist)
    }

    dist.df <- do.call(rbind.data.frame, dist.list)
    out.dist <- sparseMatrix(i=dist.df$rowIndex, j=dist.df$colIndex, x=dist.df$dist,
                             dimnames=list(rownames(in.x), rownames(in.x)),
                             giveCsparse=TRUE)
    return(out.dist)
}

#' @importFrom Matrix sparseMatrix
.aggregate_dists <- function(nhood.list, nhood.dists, cell.names){
    # we know the indices for each matrix, but these don't uniquely identify
    # entries across matrices <- replace these with the cell IDs and concatenate
    # the distance matrices into the appropriate order <- fill in missing values
    # as zeroes

    # this could be so much more elegant than this.
    dist.list <- list()
    index.df <- data.frame("cname"=cell.names, "idx"=c(1:length(cell.names)))

    for(x in seq_along(names(nhood.list))){
        x.ix <- names(nhood.list)[x]
        # summary on a triplet sparse matrix returns a matrix of i, j, x
        x.dists <- summary(nhood.dists[[x.ix]])
        x.dists$i.name <- dimnames(nhood.dists[[x.ix]])[[1]][x.dists$i]
        x.dists$j.name <- dimnames(nhood.dists[[x.ix]])[[1]][x.dists$j]
        dist.list[[paste0(x)]] <- x.dists
    }

    all.dists <- do.call(rbind.data.frame, dist.list)
    rownames(all.dists) <- NULL
    n.grid <- expand.grid(cell.names, setdiff(index.df$cname, all.dists$i.name))
    i.fake <- data.frame("i.name"=n.grid$Var1, "j.name"=n.grid$Var2, "x"=rep(0, nrow(n.grid)))
    i.fake <- merge(i.fake, index.df, by.x='i.name', by.y='cname')

    colnames(i.fake) <- c("i.name", "j.name", "x", "i")
    i.fake <- merge(i.fake, index.df, by.x='j.name', by.y='cname')
    colnames(i.fake) <- c("j.name", "i.name", "x", "i", "j")
    i.fake <- i.fake[, c("i", "j", "x", "i.name", "j.name")]

    n.grid <- expand.grid(cell.names, setdiff(index.df$cname, all.dists$j.name))
    j.fake <- data.frame("j.name"=n.grid$Var1, "i.name"=n.grid$Var2, "x"=rep(0, nrow(n.grid)))
    j.fake <- merge(j.fake, index.df, by.x='j.name', by.y='cname')
    colnames(j.fake) <- c("j.name", "i.name", "x", "j")
    j.fake <- merge(j.fake, index.df, by.x='i.name', by.y='cname')
    colnames(j.fake) <- c("i.name", "j.name", "x", "j", "i")
    j.fake <- j.fake[, c("i", "j", "x", "i.name", "j.name")]

    all.dists <- do.call(rbind.data.frame, list("in"=all.dists, "i.fake"=i.fake, "j.fake"=j.fake))

    # need to fill in the blanks
    sp.out <- sparseMatrix(i=all.dists$i, j=all.dists$j, x=all.dists$x,
                           dimnames=list(cell.names, cell.names))
    return(sp.out)

}
