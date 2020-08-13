#' Project median neighbourhood expression profiles on embedding of single-cells
#' 
#' This is for visualization of single-cells and neighbourhoods on the same embedding
#' 
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{neighbourhoods}
#' slot. 
#' @param d The number of reduced dimensions to use
#' @param reduced_dims a character indicating the name of the \code{reducedDim} slot in the
#'  \code{\linkS4class{Milo}} object to use  for projection (default: 'PCA'). \code{reducedDim(x,reduced_dim)}
#'  needs to have the attribute \code{rotation} storing the gene X dimensions loading matrix to use for projection.
#'  See details for more info.
#'
#' @details
#' This function projects neighbourhoods in the same reduced dimensions used for the full single-cell dataset.
#' First it calculates the median expression profile for all cells in a neighbourhood. Then the gene X nhood matrix is 
#' multiplied by the matrix of loadings of genes on reduced dimensions. The loading matrix should be stored as an attribute \code{rotation}
#' in \code{reducedDim(x, reduced_dims)}. This can be added by running \code{runPCA(x)} from the \code{scater} package or, using the
#' package \code{irlba} as in examples. The output of this function is used for visualization of test results.
#'
#' @return A \code{\linkS4class{matrix}} object of samples X reduced dimensions where samples are both
#' the single-cells and the neighbourhoods.
#' (This will need to become a slot in Milo object)
#' 
#' @author
#' Emma Dann
#'
#' @examples
#' 
#' requires("irlba")
#'
#' m <- matrix(rnorm(10000), ncol=10)
#' milo <- Milo(m)
#' 
#'
#' milo <- makeNeighbourhoods(milo, prop=0.1)
#' milo
#'
#' @export
#' @rdname makeNeighbourhoods
#' @importFrom BiocNeighbors findKNN
#' @importFrom igraph neighbors
#' @importFrom stats setNames
projectNeighbourhoodExpression <- function(x, d = 30, reduced_dims = "PCA"){
  if (class(x) != "Milo") {
    stop("Unrecognised input type - must be of class Milo")
  } else if (!isTRUE(.valid_neighbourhood(x))) {
    stop(
      "Not a valid Milo object - neighbourhoods are missing. Please run makeNeighbourhoods() first."
    )
  }
  
  X_reduced_dims  <- reducedDim(x, reduced_dims)
  if (d > ncol(X_reduced_dims)) {
    warning(
      paste(
        "Warning: specified d is higher than the total number of dimensions in reducedDim(x, reduced_dims). Falling back to using",
        ncol(X_reduced_dims),
        "dimensions\n"
      )
    )
    d <- ncol(X_reduced_dims)
  }
  
  if (is.null(attr(reducedDim(x), "rotation"))) {
    stop(
      "to project neighbourhoods on single-cell embedding the loading matrix for 'reduced_dim' needs to be stored as attribute 'rotation' in reducedDim(x,reduced_dim)."
    )
  }
  
  
  ## Calculate median profile of cells in a neighbourhood
  sim_milo <-
    calcNeighbourhoodExpression(sim_milo, assay = "logcounts")
  
  ## Project profiles to same PCA space of the single-cells
  X_reduced_dims <- reducedDim(x, reduced_dims)
  loadings <- attr(X_reduced_dims, "rotation")
  n.reducedDim <- t(neighbourhoodExpression(sim_milo)) %*% loadings
  
  ## Make one PC matrix including single-cells and neighbourhoods
  rownames(n.reducedDim) <-
    paste0("nh_", seq(1:nrow(n.reducedDim)))
  X_reduced_dims_merged <- rbind(n.reducedDim, X_reduced_dims)
  ## --> this will need to be added as a slot in the Milo object
  return(X_reduced_dims_merged)  
}
