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
#' @param use.assay A character scalar defining which \code{assay} slot in the
#' \code{\linkS4class{Milo}} to use
#'
#' @return A \code{\linkS4class{Milo}} object with the distance slots populated.
#'
#' @author
#' Mike Morgan, Emma Dann
#'
#' @examples
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
#' milo <- makeNhoods(milo)
#' milo <- calcNhoodDistance(milo, d=30)
#'
#' milo
#' @name calcNhoodDistance
NULL

#' @export
#' @rdname calcNhoodDistance
#' @importFrom irlba prcomp_irlba
#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix which
calcNhoodDistance <- function(x, d, reduced.dim=NULL, use.assay="logcounts"){
    if(is(x, "Milo")){
        # check for reducedDims
        if((length(reducedDimNames(x)) == 0) & is.null(reduced.dim)){
            # assume logcounts is present?
            message("Computing PCA on input")
            x_pca <- prcomp_irlba(t(assay(x, use.assay)), n=min(d+1, ncol(x)-1),
                                  scale.=TRUE, center=TRUE)
            reducedDim(x, "PCA") <- x_pca$x
        } else if((length(reducedDimNames(x)) == 0) & is.character(reduced.dim)){
            stop(paste(reduced.dim, " not found in reducedDim slot"))
        }
    } else{
        stop("Input is not a valid Milo object")
    }

    non.zero.nhoods <- which(nhoods(x)!=0, arr.ind = TRUE)

    if(any(names(reducedDims(x)) %in% c("PCA"))){
        nhood.dists <- sapply(1:ncol(nhoods(x)),
                              function(X) .calc_distance(reducedDim(x, "PCA")[non.zero.nhoods[non.zero.nhoods[,'col']==X,'row'], c(1:d),drop=FALSE]))
        names(nhood.dists) <- nhoodIndex(x)

    } else if(is.character(reduced.dim)){
        nhood.dists <- sapply(1:ncol(nhoods(x)),
                              function(X) .calc_distance(reducedDim(x, reduced.dim)[non.zero.nhoods[non.zero.nhoods[,'col']==X,'row'], c(1:d),drop=FALSE]))
        names(nhood.dists) <- nhoodIndex(x)
    }

    # ensure garbage collection
    sink(file="/dev/null")
    gc()
    sink(file=NULL)
    # all.dist <- .aggregate_dists.hard(nhoods(x), nhood.dists, colnames(x))
    # strictly this only actually needs to be a list of distance matrices
    nhoodDistances(x) <- nhood.dists

    return(x)
}


#' @importFrom Matrix rowSums sparseMatrix
#' @importFrom methods as
#' @export
.calc_distance <- function(in.x){

    dist.list <- lapply(1:nrow(in.x), FUN=function(i){
        i.dist <- apply(in.x, 1, FUN=function(P) sqrt(sum((P - in.x[i, ])**2)))
        list("rowIndex"=rep(i, nrow(in.x)), "colIndex"=c(1:length(i.dist)),
             "dist"=i.dist)
        })

    dist.df <- do.call(rbind.data.frame, dist.list)
    out.dist <- sparseMatrix(i=dist.df$rowIndex, j=dist.df$colIndex, x=dist.df$dist,
                             dimnames=list(rownames(in.x), rownames(in.x)),
                             repr="T")
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
    n.grid <- expand.grid(cell.names, cell.names)

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
                           dimnames=list(cell.names, cell.names),
                           repr="C")
    return(sp.out)

}


#' @importFrom Matrix sparseMatrix
.aggregate_dists.hard <- function(nhood.list, nhood.dists, cell.names){
    # we know the indices for each matrix, but these don't uniquely identify
    # entries across matrices <- replace these with the cell IDs and concatenate
    # the distance matrices into the appropriate order <- fill in missing values
    # as zeroes

    # this could be so much more elegant than this.
    ix <- 1
    for(x in seq_along(names(nhood.list))){
        x.ix <- names(nhood.list)[x]
        if(ix == 1){
            sp.out <- nhood.dists[[x.ix]]
        } else{
            int.id <- intersect(colnames(sp.out), colnames(nhood.dists[[x.ix]]))
            if(any(is.na(int.id))){
                stop("Corrupted dimension names")
            } else {
                # nice and easy join when no cells overlap between matrices
                i.n <- nrow(sp.out) - length(int.id)
                j.n <- ncol(nhood.dists[[x.ix]]) - length(int.id)
                if(length(int.id) > 0){
                null.mat.i <- matrix(0L, nrow=nrow(sp.out), ncol=j.n,
                                     dimnames=list(rownames(sp.out), setdiff(colnames(nhood.dists[[x.ix]]), colnames(sp.out))))

                null.mat.j <- matrix(0L, nrow=j.n, ncol=nrow(sp.out),
                                     dimnames=list(colnames(nhood.dists[[x.ix]]), rownames(sp.out)))

                .tmp.i <- cbind(sp.out, null.mat.i)
                .tmp.j <- cbind(null.mat.j, nhood.dists[[x.ix]])
                } else{
                    # remove the intersecting columns and rows
                    # construct the remaining matrices as normal
                    # extend the removed rows and columns
                    # concatenate at the endß

                }

                sp.out <- rbind(.tmp.i, .tmp.j)

                # check square
                if(ncol(sp.out) != nrow(sp.out)){
                    stop("Matrix corrupted - dimensions are not the same size")
                }
                sink(file="/dev/null")
                rm(list=c(".tmp.i", ".tmp.j"))
                gc()
                sink(file=NULL)
                print(dim(sp.out))
            }
        }
        ix <- ix + 1
    }

    sp.out <- as(sp.out, "dgCMatrix")
    return(sp.out)
}


