#' @export
#' @rdname Milo
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom Matrix dgCMatrix dsCMatrix
setClassUnion("matrixORdgCMatrixORdsCMatrix", c("matrix", "dgCMatrix", "dsCMatrix"))
setClassUnion("characterORNULL", c("character", "NULL"))

setClass("Milo",
         contains = "SingleCellExperiment",
         slots=c(
             graph = "ANY", # this should be NA or an igraph object
             adjacency = "matrixORdgCMatrixORdsCMatrix", # this should be NA or a matrix
             distance = "matrixORdgCMatrixORdsCMatrix", # this should be NA or a matrix
             neighbourhoodCounts = "matrixORdgCMatrixORdsCMatrix" # this should be NA or a matrix
         ),
         prototype = list(
             graph = NA_real_,
             adjacency = Matrix::Matrix(0L, sparse=TRUE),
             distance = Matrix::Matrix(0L, sparse=TRUE),
             neighbourhoodCounts = Matrix::Matrix(0L, sparse=TRUE)
         )
)
