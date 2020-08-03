#' @importClassesFrom Matrix dgCMatrix dsCMatrix
setClassUnion("matrixORdgCMatrixORdsCMatrix", c("matrix", "dgCMatrix", "dsCMatrix"))
setClassUnion("characterORNULL", c("character", "NULL"))


#' @rdname Milo
#' @aliases Milo
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList

setClass("Milo",
         contains = "SingleCellExperiment",
         slots=c(
             graph = "list", # this should be a list or an igraph object
             neighbourhoods = "list", # this should be a list
             neighbourDistances = "matrixORdgCMatrixORdsCMatrix", # this should be NA or a matrix
             neighbourhoodCounts = "matrixORdgCMatrixORdsCMatrix" # this should be NA or a matrix
         ),
         prototype = list(
             graph = list(),
             neighbourhoods = list(),
             neighbourDistances = Matrix::Matrix(0L, sparse=TRUE),
             neighbourhoodCounts = Matrix::Matrix(0L, sparse=TRUE)
         )
)
