#' The Milo class
#'
#' The Milo class extends the SingleCellExperiment class and is designed to
#' work with neighbourhoods of cells. Therefore, it inherits from the
#' \linkS4class{SingleCellExperiment} class and follows the same usage
#' conventions. There is additional support for cell-to-cell distances
#' via \code{\link{distance}}, and the KNN-graph used to define the
#' neighbourhoods.
#'
#' @param ... Arguments passed to the \code{\link{Milo}} constructor to fill
#' the slots of the base class.
#' @param graph An igraph object or list of adjacent vertices that represents
#' the KNN-graph
#' @param adjacency A sparse matrix representation of the adjacency list
#' @param distance A sparse matrix of cell-to-cell distances for cells in the
#' same neighbourhoods
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
#' @export
#'
#' @importFrom SingleCellExperiment SingleCellExperiment


Milo <- function(...,
                 graph=NA_real_,
                 adjacency=Matrix::Matrix(0L, sparse=TRUE),
                 distance=Matrix::Matrix(0L, sparse=TRUE),
                 neighbourhoodCounts=Matrix::Matrix(0L, sparse=TRUE)){
    old <- S4Vectors:::disableValidity()
    if (!isTRUE(old)) {
        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(old))
    }

    if(class(unlist(...)) == "SingleCellExperiment"){
        milo <- .fromSCE(unlist(...))
    }

    milo
}

## class validator
setValidity("Milo", function(object){
    if (class(object@neighbourhoodCounts) != "matrixORdgCMatrixORdsCMatrix"){
        "@neighbourhoodCounts must be a matrix format"
    } else{
        TRUE
    }

    if (class(object@graph) != "numeric"){
        "@graph must be of type numeric"
    } else{
        TRUE
    }
})
