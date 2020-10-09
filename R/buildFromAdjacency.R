#' Build a graph from an input adjacency matrix
#'
#' Construct a kNN-graph from an input adjacency matrix - either binary or distances between NNs.
#' @param x An n X n \code{matrix} of single-cells, where values represent edges between cells; 0 values
#' are taken to mean no edge between cells. If the matrix is not binary, then it is assumed the values
#' are distances; 0 retain the same meaning. This behaviour can be toggled using \code{is.binary=TRUE}.
#' @param k (optional) Scalar value that represents the number of nearest neighbours in the original graph. This
#' can also be inferred directly from the adjacency matrix \code{x}.
#' @param is.binary Logical scalar indicating if the input matrix is binary or not.
#'
#' @details
#' This function will take a matrix as input and construct the kNN graph that it describes. If
#' the matrix is not symmetric then the graph is assumed to be directed, whereas if the matrix
#' is not binary, i.e. all 0's and 1's then the input values are taken to be distances between
#' graph vertices; 0 values are assumed to represent a lack of edge between vertices.
#'
#' @return A \code{\linkS4class{Milo}} with the graph slot populated.
#' @author Mike Morgan
#' @examples
#' r <- 1000
#' c <- 1000
#' k <- 35
#' m <- floor(matrix(runif(r*c), r, c))
#' for(i in seq_along(1:r)){
#'     m[i, sample(1:c, size=k)] <- 1
#' }
#'
#' milo <- buildFromAdjacency(m)
#'
#' @name buildFromAdjacency
NULL

#' @export
#' @importFrom Matrix rowSums
#' @importFrom igraph graph_from_adjacency_matrix
buildFromAdjacency <- function(x, k=NULL, is.binary=NULL, ...){
    # test class? If it's not sparse then cast it to one
    if(!is(x, "sparseMatrix")){
        if(class(x) %in% c("matrix")){
            message("Casting to sparse matrix format")
            x <- as(x, "dgTMatrix")
        } else{
            stop(paste0("Input 'x' is not a recognisable matrix format. Class is ", class(x)))
        }
    }

    if(is.null(k)){
        message("Inferring k from matrix")
        k.val <- unique(rowSums(x > 0))
        if(length(k.val) > 1){
            warning("Row sums are not all equal, checking columns")
            k.val <- unique(colSums(x > 0))
            if(length(k.val) > 1){
                stop("Cannot infer k from matrix, please provide a value or check input")
            }
        }
    } else{
        k.val <- k
    }

    # check if matrix is binary
    if(is.null(is.binary)){
        is.binary <- .check_binary(x)
    }

    # check if square
    is.square <- ifelse(nrow(x) == ncol(x), TRUE, FALSE)

    # use igraph if it's square
    if(is.square){
        if(!is.binary){
            bin.x <- as(matrix(as.numeric(x > 0), nrow=nrow(x)), "dgCMatrix")
            nn.graph <- graph_from_adjacency_matrix(bin.x, mode="undirected",
                                                    weighted=NULL,
                                                    diag=FALSE)
        } else{
            nn.graph <- graph_from_adjacency_matrix(x, mode="undirected",
                                                    weighted=NULL,
                                                    diag=FALSE)
        }
    } else{
        if(is.binary){
            stop("Input matrix is binary but not square")
        }
        # assume the #ncols is k and the individual components are the NN indices
        graph <- .neighborsToKNNGraph(x, directed=FALSE)
    }

    mylo <- Milo()
    graph(mylo) <- nn.graph

    # if the matrix contains distances then also population the nhoodDistances slot
    if(!is.binary & is.square){
        message("Adding nhoodDistances to Milo object")
        nhoodDistances(mylo) <- NULL
    }

    return(mylo)
}
