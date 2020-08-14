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
#' @param scale logical indicating whether neighbourhoods expression profiles should be scaled before projection. 
#' Data should be scaled if the embedding was performed on scaled single-cell data. 
#' (Default: TRUE)
#' @param center logical indicating whether neighbourhoods expression profiles should be centered before projection
#' Data should be centered if the embedding was performed on scaled single-cell data.
#' (Default: TRUE)
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
#' require(SingleCellExperiment)
#' m <- matrix(rnorm(100000), ncol=100)
#' milo <- Milo(SingleCellExperiment(assays(logcounts=m)))
#' milo <- buildGraph(m, d=30, transposed=TRUE)
#' milo <- makeNeighbourhoods(milo)
#' milo <- projectNeighbourhoodExpression(milo)
#' 
#' @export
#' @rdname projectNeighbourhoodExpression
projectNeighbourhoodExpression <- function(x, d = 30, reduced_dims = "PCA", scale=TRUE){
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
  
  ## Calculate mean profile of cells in a neighbourhood
  if (is.null(neighbourhoodExpression(sim_milo))) {
    x <- calcNeighbourhoodExpression(x, assay = "logcounts")
  }
  
  ## Scale the neighbourhoods profile (if PCs were calculated on scaled data)
  if (isTRUE(scale) | isTRUE(center)) {
    neighbourhoodExpression(x) <- t(scale(t(neighbourhoodExpression(x)), scale=TRUE, center=TRUE))
  }
  
  ## Project profiles to same PCA space of the single-cells
  X_reduced_dims <- reducedDim(x, reduced_dims)
  loadings <- attr(X_reduced_dims, "rotation")
  n.reducedDim <- t(neighbourhoodExpression(x)) %*% loadings
  
  ## Make one PC matrix including single-cells and neighbourhoods
  rownames(n.reducedDim) <-
    paste0("nh_", seq(1:nrow(n.reducedDim)))
  X_reduced_dims_merged <- rbind(n.reducedDim, X_reduced_dims)
  ## --> this will need to be added as a slot in the Milo object
  return(X_reduced_dims_merged)  
}

#' Calculates the l2-norm of a vector
#'
#' Adapted from PMA package
#' @references Witten, Tibshirani, and Hastie, Biostatistics 2009
#' @references \url{https://github.com/cran/PMA/blob/master/R/PMD.R}
#'
#' @param vec numeric vector
#'
#' @return returns the l2-norm.
#'
.l2norm <- function(vec) {
  a <- sqrt(x = sum(vec ^ 2))
  if (a == 0) {
    a <- .05
  }
  return(a)
}
