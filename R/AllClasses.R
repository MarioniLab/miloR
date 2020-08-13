
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
             nhoods = "list", # this should be a list
             nhoodDistances = "matrixORdgCMatrixORdsCMatrix", # this should be NA or a matrix
             nhoodCounts = "matrixORdgCMatrixORdsCMatrix", # this should be NA or a matrix
             nhoodIndex = "list", # used to store nhood indices
             nhoodExpression = "matrixORdgCMatrixORdsCMatrix" # this should be NA or a matrix
         ),
         prototype = list(
             graph = list(),
             nhoods = list(),
             nhoodDistances = Matrix::Matrix(0L, sparse=TRUE),
             nhoodCounts = Matrix::Matrix(0L, sparse=TRUE),
             nhoodIndex = list(),
             nhoodExpression = Matrix::Matrix(0L, sparse=TRUE)
         )
)
