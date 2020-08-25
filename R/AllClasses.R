
#' @importClassesFrom Matrix dgCMatrix dsCMatrix dgTMatrix dgeMatrix sparseMatrix
setClassUnion("matrixORMatrix", c("matrix", "dgCMatrix", "dsCMatrix", "dgTMatrix", "dgeMatrix")) # is there a record for how long a virtual class can be?!
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
             nhoodDistances = "matrixORMatrix", # this should be a matrix
             nhoodCounts = "matrixORMatrix", # this should be a matrix
             nhoodIndex = "list", # used to store nhood indices
             nhoodExpression = "matrixORMatrix", # this should be NA or a matrix
             nhoodReducedDim = "list" # this should be a list
         ),
         prototype = list(
             graph = list(),
             nhoods = list(),
             nhoodDistances = Matrix::Matrix(0L, sparse=TRUE),
             nhoodCounts = Matrix::Matrix(0L, sparse=TRUE),
             nhoodIndex = list(),
             nhoodExpression = Matrix::Matrix(0L, sparse=TRUE),
             nhoodReducedDim = list()
         )
)
