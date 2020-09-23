#' Build a k-nearest neighbour graph
#'
#' This function is borrowed from the old buildKNNGraph function in scran.
#' Instead of returning an igraph object it populates the graph and distance
#' slots in a Milo object. If the input is a SingleCellExperiment object or
#' a matrix then it will return a de novo Milo object with the same slots
#' filled.
#' @param x A matrix, \code{\linkS4class{SingleCellExperiment}} or Milo object
#' containing feature X cell gene expression data.
#' @param k An integer scalar that specifies the number of nearest-neighbours
#' to consider for the graph building.
#' @param d The number of dimensions to use if the input is a matrix of cells
#' X reduced dimensions. If this is provided, transposed should also be
#' set=TRUE.
#' @param transposed Logical if the input x is transposed with rows as cells.
#' @param BNPARAM refer to \code{\link[scran]{buildKNNGraph}} for details.
#' @param BSPARAM refer to \code{\link[scran]{buildKNNGraph}} for details.
#' @param BPPARAM refer to \code{\link[scran]{buildKNNGraph}} for details.
#' @param seed Seed number used for pseudorandom number generators.
#'
#' @details
#' This function computes a k-nearest neighbour graph. Each graph vertex is a
#' single-cell connected by the edges between its neighbours. Whilst a
#' kNN-graph is strictly directed, we remove directionality by forcing all
#' edge weights to 1; this behaviour can be overriden by providing
#' \code{directed=TRUE}.
#'
#' If you wish to use an
#' alternative graph structure, such as a shared-NN graph I recommend you
#' construct this separately and add to the relevant slot in the
#' \code{\link{Milo}} object.
#'
#' @return A \code{\linkS4class{Milo}} object with the graph and distance slots populated.
#'
#' @author
#' Mike Morgan, with KNN code written by Aaron Lun & Jonathan Griffiths.
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
#'
#' milo
#' @name buildGraph
NULL

#' @export
#' @rdname buildGraph
#' @importFrom irlba prcomp_irlba
#' @importFrom BiocSingular bsparam
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
buildGraph <- function(x, k=10, d=50, transposed=FALSE, BNPARAM=KmknnParam(),
                       BSPARAM=bsparam(), BPPARAM=SerialParam(), seed=42){
    set.seed(seed)
    # check class of x to determine which function to call
    # in all cases it must return a Milo object with the graph slot populated
    # what is a better design principle here? make a Milo object here and just
    # have one function, or have a separate function for input data type? I
    # think the former probably.

    if(class(x) == "Milo"){
        # check for reducedDims
        if(is.null(reducedDim(x))){
            # assume logcounts is present?
            x_pca <- prcomp_irlba(t(logcounts(x)), n=min(d+1, ncol(x)-1),
                                  scale.=TRUE, center=TRUE)
            reducedDim(x, "PCA") <- x_pca$x
        } else if(!any(names(reducedDims(x)) %in% c("PCA"))){
            # assume logcounts is present?
            x_pca <- prcomp_irlba(t(logcounts(x)), n=min(d+1, ncol(x)-1),
                                  scale.=TRUE, center=TRUE)
            reducedDim(x, "PCA") <- x_pca$x
        }
    } else if(class(x) == "matrix" & isTRUE(transposed)){
        # assume input are PCs - the expression data is non-sensical here
        SCE <- SingleCellExperiment(assays=list(counts=Matrix(0L, nrow=1, ncol=nrow(x))),
                                    reducedDims=SimpleList("PCA"=x))
        x <- Milo(SCE)
    } else if(class(x) == "matrix" & isFALSE(transposed)){
        # this should be a gene expression matrix
        SCE <- SingleCellExperiment(assays=list(logcounts=x))
        x_pca <- prcomp_irlba(t(logcounts(SCE)), n=min(d+1, ncol(x)-1),
                              scale.=TRUE, center=TRUE)
        reducedDim(SCE, "PCA") <- x_pca$x
        x <- Milo(SCE)
    } else if (class(x) == "SingleCellExperiment"){
        # test for reducedDims, if not then compute them
        # give me a Milo object
        if(is.null(reducedDim(x))){
            # assume logcounts is present - how dangerous is this?
            # better to check first, or have the user input the assay
            # to use?
            x_pca <- prcomp_irlba(t(logcounts(x)), n=min(d+1, ncol(x)-1),
                                  scale.=TRUE, center=TRUE)
            reducedDim(x, "PCA") <- x_pca$x
        }

        x <- Milo(x)
    }

    .buildGraph(x, k=k, d=d,
                BNPARAM=BNPARAM, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
}


#' @export
#' @importFrom Matrix Matrix
#' @importFrom BiocSingular bsparam
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
.buildGraph <- function(x, k=10, d=50,
                        BNPARAM=KmknnParam(), BSPARAM=bsparam(),
                        BPPARAM=SerialParam()){

    nn.out <- .setup_knn_data(x=reducedDim(x, "PCA"), d=d,
                              k=k, BNPARAM=BNPARAM, BSPARAM=BSPARAM,
                              BPPARAM=BPPARAM)
    sink(file="/dev/null")
    gc()
    sink(file=NULL)

    # separate graph and distances? At some point need to expand the distances
    # to the larger neighbourhood
    message(paste0("Constructing kNN graph with k:", k))
    zee.graph <- .neighborsToKNNGraph(nn.out$index, directed=FALSE)
    graph(x) <- zee.graph

    # adding distances
    message(paste0("Retrieving distances from ", k, " nearest neighbours"))
    # set this up as a dense matrix first, then coerce to a sparse matrix
    # starting with a sparse matrix requires a coercion at each iteration
    # which uses up lots of memory and unncessary CPU time
    old.dist <- matrix(0L, ncol=ncol(x), nrow=ncol(x))

    n.idx <- ncol(x)
    for(i in seq_along(1:n.idx)){
        sink(file="/dev/null")
        gc()
        sink(file=NULL)

        i.knn <- nn.out$index[i, ]
        i.dists <- nn.out$distance[i, ]
        old.dist[i, i.knn] <- i.dists
        old.dist[i.knn, i] <- i.dists
    }
    old.dist <- as(old.dist, "dgCMatrix")
    nhoodDistances(x) <- old.dist

    x
}


#' @importFrom BiocNeighbors findKNN
.setup_knn_data <- function(x, k, d=50,
                            BNPARAM, BSPARAM, BPPARAM) {

    # Finding the KNNs - keep the distances
    # input should be cells X dimensions
    findKNN(x[, c(1:d)], k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=TRUE)
}


