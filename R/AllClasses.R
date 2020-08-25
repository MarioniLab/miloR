
#' @importClassesFrom Matrix sparseMatrix
#setClassUnion("matrixORdgCMatrixORdsCMatrixORdgTMatrix", c("matrix", "dgCMatrix", "dsCMatrix", "dgTMatrix")) # is there a record for how long a virtual class can be?!
setClassUnion("matrixORsparseMatrix", c("matrix", "sparseMatrix")) # can I generalise this much?
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
             nhoodDistances = "matrixORsparseMatrix", # this should be NA or a matrix
             nhoodCounts = "matrixORsparseMatrix", # this should be NA or a matrix
             nhoodIndex = "list", # used to store nhood indices
             nhoodExpression = "matrixORsparseMatrix", # this should be NA or a matrix
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
