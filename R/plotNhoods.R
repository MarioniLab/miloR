### Plotting utils for neighbourhoods

#' Plot histogram of neighbourhood sizes
#'
#' This function plots the histogram of the number of cells belonging to
#' each neighbourhood
#'
#' @param milo A \code{\linkS4class{Milo}} object with a non-empty \code{nhoods}
#' slot.
#' @param bins number of bins for \code{geom_histogram}
#'
#' @return A \code{\linkS4class{ggplot}} object
#'
#' @author
#' Emma Dann
#'
#' @examples
#'
#' requires(igraph)
#' m <- matrix(rnorm(10000), ncol=10)
#' milo <- buildGraph(m, d=10)
#'
#' milo <- makeNhoods(milo, prop=0.1)
#' plotNhoodSizeHist(milo)
#'
#' @export
#' @rdname plotNhoodSizeHist
#' @importFrom ggplot2 ggplot geom_histogram xlab theme_classic
#' @importFrom igraph neighbors
plotNhoodSizeHist <- function(milo, bins=50){
  if (! isTRUE(.valid_nhood(milo))){
    stop("Not a valid Milo object - neighbourhoods are missing. Please run makeNhoods() first.")
  }
  df <- data.frame(nh_size=sapply(nhoods(milo), function(x) length(x)))
  ggplot(data=df, aes(nh_size)) + geom_histogram(bins=bins) +
    xlab("Neighbourhood size") +
    theme_classic(base_size = 16)
}


#' @importFrom igraph is_igraph
.valid_nhood <- function(milo){
  # check for a valid nhood slot
  n_neigh <- length(nhoods(milo))
  is_not_empty <- n_neigh > 0
  is_igraph_vx <- class(milo@nhoods[[sample(1:n_neigh, 1)]]) == "igraph.vs"
  if (isTRUE(is_igraph_vx & is_not_empty)){
    TRUE
  } else {
    FALSE
  }
}
