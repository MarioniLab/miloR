### Neighbourhood abstracted graph ###

#' Build an abstracted graph of neighbourhoods for visualization
#'
#'
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{nhoods}
#' slot.
#'
#' @details
#' This constructs a graph where nodes represent neighbourhoods and edges represent the number of overlapping
#' cells between two neighbourhoods.
#'
#' @return A \code{\linkS4class{Milo}} object with filled
#'
#' @author
#' Emma Dann
#'
#' @examples
#'
#' NULL
#'
#' @import igraph
#' @export
#' @rdname buildNhoodGraph
buildNhoodGraph <- function(x){
  # are neighbourhoods calculated?
  if(length(nhoods(x)) == 0){
    stop("No neighbourhoods found - run makeNhoods first")
  }
  ## Build adjacency matrix for nhoods
  nh_intersect_mat <- .build_nhood_adjacency(nhoods(x))
  ## Make igraph object
  ig <- graph.adjacency(nh_intersect_mat, mode="undirected", weighted=TRUE)
  nhood_sizes <- sapply(nhoods(x), length)
  ig <- set_vertex_attr(ig, name = 'size', value = nhood_sizes[vertex.attributes(ig)$name])
  ## Add to nhoodGraph slot in milo object
  nhoodGraph(x) <- ig
  return(x)
}


# Build adjacency matrix of overlap between neighbourhoods
#' @importFrom gtools permutations
.build_nhood_adjacency <- function(nhoods){
  nms <- permutations(n = length(nhoods), v = names(nhoods), r = 2, repeats.allowed = T)
  print("Calculating nhood adjacency....")
  out <- sapply( 1:nrow(nms) , function(x) length( intersect( nhoods[[nms[x,1]]], nhoods[[ nms[x,2] ]]) ) )
  nh_intersect_mat <- matrix(out, nrow = length(nhoods), byrow = TRUE)
  rownames(nh_intersect_mat) <- unique(nms[,1])
  colnames(nh_intersect_mat) <- unique(nms[,2])
  return(nh_intersect_mat)
}

# ## A little test
# nh1 <- sample(names(nhoods))[1]
# nh2 <- sample(names(nhoods))[1]
# nh_intersect_mat[nh1, nh2] == length(intersect(nhoods(sim_milo)[[nh1]], nhoods(sim_milo)[[nh2]]))
