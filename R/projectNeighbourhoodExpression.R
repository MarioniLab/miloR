#' Project mean nhood expression profiles on embedding of single-cells
#' 
#' This is for visualization of single-cells and nhoods on the same embedding
#' 
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{nhoods}
#' slot. 
#' @param d The number of reduced dimensions to use
#' @param reduced_dims a character indicating the name of the \code{reducedDim} slot in the
#'  \code{\linkS4class{Milo}} object to use  for projection (default: 'PCA'). \code{reducedDim(x,reduced_dim)}
#'  needs to have the attribute \code{rotation} storing the gene X dimensions loading matrix to use for projection.
#'  See details for more info.
#' @param scale logical indicating whether nhoods expression profiles should be scaled before projection. 
#' Data should be scaled if the embedding was performed on scaled single-cell data. 
#' (Default: TRUE)
#' @param center logical indicating whether nhoods expression profiles should be centered before projection
#' Data should be centered if the embedding was performed on scaled single-cell data.
#' (Default: TRUE)
#'
#' @details
#' This function projects nhoods in the same reduced dimensions used for the full single-cell dataset.
#' First it calculates the mean expression profile for all cells in a nhood. Then the gene X nhood matrix is 
#' multiplied by the matrix of loadings of genes on reduced dimensions. The loading matrix should be stored as an attribute \code{rotation}
#' in \code{reducedDim(x, reduced_dims)}. This can be added by running \code{runPCA(x)} from the \code{scater} package or, using the
#' package \code{irlba} as in examples. The output of this function is used for visualization of test results.
#'
#' @return A \code{\linkS4class{matrix}} object of samples X reduced dimensions where samples are both
#' the nhoods and the single-cells.
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
#' milo <- makeNhoods(milo)
#' milo <- projectNhoodExpression(milo)
#' 
#' @export
#' @rdname projectNhoodExpression
projectNhoodExpression <- function(x, d = 30, reduced_dims = "PCA", scale=TRUE){
  if (class(x) != "Milo") {
    stop("Unrecognised input type - must be of class Milo")
  } else if (!isTRUE(.valid_nhood(x))) {
    stop(
      "Not a valid Milo object - nhoods are missing. Please run makeNhoods() first."
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
      "to project nhoods on single-cell embedding the loading matrix for 'reduced_dim' needs to be stored as attribute 'rotation' in reducedDim(x,reduced_dim)."
    )
  }
  
  ## Calculate mean profile of cells in a nhood
  if (is.null(nhoodExpression(sim_milo))) {
    x <- calcNhoodExpression(x, assay = "logcounts")
  }
  
  ## Scale the nhoods profile (crucial if PCs were calculated on scaled data)
  if (isTRUE(scale) | isTRUE(center)) {
    nhoodExpression(x) <- t(scale(t(nhoodExpression(x)), scale=TRUE, center=TRUE))
  }
  
  ## Project profiles to same PCA space of the single-cells
  X_reduced_dims <- reducedDim(x, reduced_dims)
  loadings <- attr(X_reduced_dims, "rotation")
  n.reducedDim <- t(nhoodExpression(x)) %*% loadings
  
  ## Make one PC matrix including single-cells and nhoods
  rownames(n.reducedDim) <-
    paste0("nh_", seq(1:nrow(n.reducedDim)))
  X_reduced_dims_merged <- rbind(n.reducedDim, X_reduced_dims)
  
  ## Add to slot nhoodsReducedDim
  nhoodReducedDim(x) <- list(X_reduced_dims_merged)
  names(nhoodReducedDim(x)) <- reduced_dims
  return(x)  
}


