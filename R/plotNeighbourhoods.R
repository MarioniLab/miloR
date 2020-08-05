### Plotting utils for neighbourhoods

#' Plot histogram of neighbourhood sizes
#' 
#' This function plots the histogram of the number of cells belonging to 
#' each neighbourhood 
#' 
#' @param milo A \code{\linkS4class{Milo}} object with a non-empty \code{neighbourhoods}
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
#' milo <- makeNeighbourhoods(milo, prop=0.1)
#' plotNeighborhoodSizeHist(milo)
#' 
#' @export
#' @rdname plotNeighborhoodSizeHist
#' @importFrom ggplot2 ggplot geom_histogram xlab theme_classic
#' @importFrom igraph neighbors
plotNeighborhoodSizeHist <- function(milo, bins=50){
  if (! isTRUE(.valid_neighbourhood(milo))){
    stop("Not a valid Milo object - neighbourhoods are missing. Please run makeNeighbourhoods() first.")
  }
  df <- data.frame(nh_size=sapply(neighbourhoods(milo), function(x) length(x))) 
  ggplot(data=df, aes(nh_size)) + geom_histogram(bins=bins) +
    xlab("Neighbourhood size") +
    theme_classic(base_size = 16)
}


#' @importFrom igraph is_igraph
.valid_neighbourhood <- function(milo){
  # check for a valid neighbourhood slot
  n_neigh <- length(neighbourhoods(milo))
  is_not_empty <- n_neigh > 0
  is_igraph_vx <- class(milo@neighbourhoods[[sample(1:n_neigh, 1)]]) == "igraph.vs" 
  if (isTRUE(is_igraph_vx & is_not_empty)){
    TRUE
  } else {
    FALSE
  }
}
