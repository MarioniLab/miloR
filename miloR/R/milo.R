#' The Milo class
#'
#' The Milo class extends the SingleCellExperiment class and is designed to
#' work with neighbourhoods of cells. Therefore, it inherits from the
#' \linkS4class{SingleCellExperiment} class and follows the same usage
#' conventions. There is additional support for cell-to-cell distances
#' via distance, and the KNN-graph used to define the neighbourhoods.
#'
#' @param ... Arguments passed to the Milo constructor to fill the slots of the
#' base class. This should be either a \linkS4class{SingleCellExperiment} or
#' matrix of features X cells
#' @param graph An igraph object or list of adjacent vertices that represents
#' the KNN-graph
#' @param nhoods A list of graph vertices, each containing the indices
#' of the constiuent graph vertices in the respective neighbourhood
#' @param nhoodDistances A sparse matrix of cell-to-cell distances for
#' cells in the same neighbourhoods
#' @param nhoodCounts A matrix of neighbourhood X sample counts of the
#' number of cells in each neighbourhood derived from the respective samples
#' @param nhoodIndex A list of cells that are the neighborhood index cells.
#' @param nhoodExpression A matrix of gene X neighbourhood expression.
#'
#' @details
#' In this class the underlying structure is the gene/feature X cell expression
#' data. The additional slots provide a link between these single cells and
#' the neighbourhood representation. This can be further extended by the use
#' of an abstracted graph for visualisation that preserves the structure of the
#' single-cell KNN-graph
#'
#' A Milo object can also be constructed by inputting a feature X cell gene
#' expression matrix. In this case it simply constructs a SingleCellExperiment
#' and fills the relevant slots, such as reducedDims.
#'
#' @return a Milo object
#'
#' @author Mike Morgan
#'
#' @examples
#'
#' library(SingleCellExperiment)
#' ux <- matrix(rpois(12000, 5), ncol=200)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#'
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'   reducedDims=SimpleList(PCA=pca$x))
#'
#' milo <- Milo(sce)
#' milo
#'
#' @docType class
#' @name Milo
NULL

#' @export
#' @rdname Milo
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix Matrix
Milo <- function(..., graph=list(), nhoodDistances=Matrix(0L, sparse=TRUE),
                 nhoods=list(),
                 nhoodCounts=Matrix(0L, sparse=TRUE),
                 nhoodIndex=list(),
                 nhoodExpression=Matrix(0L, sparse=TRUE)){
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    if(length(list(...)) == 0){
        milo <- .emptyMilo()
    } else if(class(unlist(...)) == "SingleCellExperiment"){
        milo <- .fromSCE(unlist(...))
    }

    milo
}


#' @importFrom S4Vectors SimpleList
#' @importFrom Matrix Matrix
#' @import SingleCellExperiment
.fromSCE <- function(sce, assayName="logcounts"){
    # make the distance and adjacency matrices the correct size
    out <- new("Milo", sce,
               graph=list(),
               nhoods=list(),
               nhoodDistances=Matrix(0L, sparse=TRUE),
               nhoodCounts=Matrix(0L, sparse=TRUE),
               nhoodIndex=list(),
               nhoodExpression=Matrix(0L, sparse=TRUE))

    reducedDims(out) <- reducedDims(sce)
    altExps(out) <- list()

    out
}

#' @importFrom Matrix Matrix
.emptyMilo <- function(...){
    # return an empty Milo object
    out <- new("Milo",
               graph=list(),
               nhoods=list(),
               nhoodDistances=Matrix(0L, sparse=TRUE),
               nhoodCounts=Matrix(0L, sparse=TRUE),
               nhoodIndex=list(),
               nhoodExpression=Matrix(0L, sparse=TRUE))

    reducedDims(out) <- list()
    altExps(out) <- list()

    out
}


## class validator
#' @importFrom igraph is_igraph
setValidity("Milo", function(object){
    if (class(object@nhoodCounts) != "matrixORMatrix"){
        "@nhoodCounts must be a matrix format"
    } else{
        TRUE
    }

    if(class(object@nhoodDistances) != "matrixORMatrix"){
        "@nhoodDistances must be a matrix format"
    } else{
        TRUE
    }

    if(class(object@nhoodExpression) != "matrixORMatrix"){
        "@nhoodExpression must be a matrix format"
    } else{
        TRUE
    }

    # can be a list or igraph object
    if (!is_igraph(object@graph)){
        if(typeof(object@graph) != "list"){
            "@graph must be of type list or igraph"
        }
        } else{
            TRUE
    }
})
