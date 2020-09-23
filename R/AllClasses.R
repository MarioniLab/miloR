#' The Milo container class
#'
#' @slot graph An igraph object that represents the kNN graph
#' @slot nhoods A list of neighbourhoods as graph indices and their constituent single cells
#' @slot nhoodDistances An NxN sparse matrix of Euclidean distances between vertices in each neighbourhood
#' @slot nhoodCounts An NxM sparse matrix of cells counts in each neighourhood across M samples
#' @slot nhoodIndex A list of the index vertices for each neighbourhood
#' @slot nhoodExpression An GxN matrix of genes X neighbourhoods containing average gene expression levels across cells in each neighbourhood
#' @slot nhoodReducedDim a list of reduced dimensional representations of neighbourhoods, including projections into lower dimension space
#' @slot nhoodGraph an igraph object that represents the graph of neighbourhoods
#'

#' @importClassesFrom Matrix dgCMatrix dsCMatrix dgTMatrix dgeMatrix sparseMatrix
setClassUnion("matrixORMatrix", c("matrix", "dgCMatrix", "dsCMatrix", "dgTMatrix", "dgeMatrix")) # is there a record for how long a virtual class can be?!
setClassUnion("characterORNULL", c("character", "NULL"))
#' @aliases Milo
#' @rdname Milo
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
setClass("Milo",
         contains = "SingleCellExperiment",
         slots=c(
             graph = "list", # this should be a list or an igraph object
             nhoods = "list", # this should be a list
             nhoodDistances = "matrixORMatrix", # this should be a matrix
             nhoodCounts = "matrixORMatrix", # this should be a matrix
             nhoodIndex = "list", # used to store nhood indices
             nhoodExpression = "matrixORMatrix", # this should be NA or a matrix
             nhoodReducedDim = "list", # this should be a list
             nhoodGraph = "list" # this should be an igraph object (I'm copying from the graph slot)
         ),
         prototype = list(
             graph = list(),
             nhoods = list(),
             nhoodDistances = Matrix::Matrix(0L, sparse=TRUE),
             nhoodCounts = Matrix::Matrix(0L, sparse=TRUE),
             nhoodIndex = list(),
             nhoodExpression = Matrix::Matrix(0L, sparse=TRUE),
             nhoodReducedDim = list(),
             nhoodGraph = list()
         )
)
