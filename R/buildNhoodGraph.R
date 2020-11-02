### Neighbourhood abstracted graph ###

#' Build an abstracted graph of neighbourhoods for visualization
#'
#'
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{nhoods}
#' slot.
#' @param overlap A numeric scalar that thresholds graph edges based on  the number
#' of overlapping cells between neighbourhoods.
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
#' @importFrom igraph graph.adjacency set_vertex_attr vertex.attributes
#' @export
#' @rdname buildNhoodGraph
buildNhoodGraph <- function(x, overlap=1){

    if(!is(x, "Milo")){
      stop("Not a valid Milo object")
    }

    # are neighbourhoods calculated?
    if(length(nhoods(x)) == 0){
    stop("No neighbourhoods found - run makeNhoods first")
    }

    ## Build adjacency matrix for nhoods
    nh_intersect_mat <- .build_nhood_adjacency(nhoods(x))

    # add to slot if empty
    nhoodAdjacency(x) <- nh_intersect_mat

    ## Make igraph object
    ig <- graph.adjacency(nh_intersect_mat, mode="undirected", weighted=TRUE)
    nhood_sizes <- sapply(nhoods(x), length)
    ig <- set_vertex_attr(ig, name = 'size', value = nhood_sizes[vertex.attributes(ig)$name])
    ## Add to nhoodGraph slot in milo object
    nhoodGraph(x) <- ig
    return(x)
}


