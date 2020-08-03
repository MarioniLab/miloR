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
#' @return A \code{\linkS4class{Milo}} object with the graph, adjacency and
#' distance slots populated.
#'
#' @author
#' Mike Morgan, with KNN code written by Aaron Lun & Jonathan Griffiths.
#'
#' @examples
#' m <- matrix(rnorm(10000), ncol=10)
#' milo <- buildGraph(m, d=10)
#'
#' milo
#' @name buildGraph
NULL

#' @export
#' @rdname buildGraph
#' @importFrom irlba prcomp_irlba
buildGraph <- function(x, k=10, d=50, transposed=FALSE, BNPARAM=KmknnParam(),
                       BSPARAM=bsparam(), BPPARAM=SerialParam()){

    # check class of x to determine which function to call
    # in all cases it must return a Milo object with the graph slot populated
    # what is a better design principle here? make a Milo object here and just
    # have one function, or have a separate function for input data type? I
    # think the former probably.

    if(class(x) == "Milo"){
        # check for reducedDims
        if(is.null(reducedDim(x, "PCA"))){
            # assume logcounts is present?
            x_pca <- prcomp_irlba(logcounts(x))
            reducedDim(x, "PCA") <- x_pca$x
        }
    } else if (class(x) == "matrix"){
        # assume input are PCs


    } else if (class(x) == "SingleCellExperiment"){
        # test for reducedDims, if not then compute them
        # give me a Milo object
        if(is.null(reducedDim(x))){
            # assume logcounts is present - how dangerous is this?
            # better to check first, or have the user input the assay
            # to use?
            x_pca <- prcomp_irlba(logcounts(x))
            reducedDim(x, "PCA") <- x_pca$x
        }

        x <- Milo(x)
    }

    .buildGraph(x, k=k, d=d, transposed=transposed,
                subset.row=NULL,
                BNPARAM=BNPARAM, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
}

#' @importFrom Matrix Matrix
.buildGraph <- function(x, k=10, d=50, transposed=transposed,
                        subset.row=subset.row,
                        BNPARAM=KmknnParam(), BSPARAM=bsparam(),
                        BPPARAM=SerialParam()){

    nn.out <- .setup_knn_data(x=reducedDim(x, "PCA"), subset.row=subset.row,
                              d=d, transposed=transposed,
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
    old.dist <- Matrix(0L, ncol=ncol(x), nrow=ncol(x), sparse=TRUE)

    n.idx <- ncol(x)
    for(i in seq_along(1:n.idx)){
        i.knn <- nn.out$index[i, ]
        i.dists <- nn.out$distance[i, ]
        old.dist[i, i.knn] <- i.dists
        old.dist[i.knn, i] <- i.dists
    }
    neighbourDistances(x) <- old.dist

    x
}


#' @importFrom BiocNeighbors findKNN
.setup_knn_data <- function(x, subset.row, d, transposed, k,
                            BNPARAM, BSPARAM, BPPARAM) {

    # Finding the KNNs - keep the distances
    # input should be cells X dimensions
    findKNN(x, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=TRUE)
}



#' @import BiocNeighbors
#' @import igraph
#' @importFrom reshape2 melt
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
