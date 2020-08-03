#' Define neighbourhoods on a graph
#'
#' This function randomly samples vertcies on a graph to define neighbourhoods.
#' These are then refined by computing the median profile for the neighbourhood
#' in reduced dimensional space and selecting the nearest vertex to this
#' position. Thus, multiple neighbourhoods may be collapsed down together to
#' prevent over-sampling the graph space.
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{graph}
#' slot. Alternatively an \code{igraph} object on which neighbourhoods will
#' be defined.
#' @param prop A double scalar that defines what proportion of graph vertices
#' to randomly sample. Must be 0 < prop < 1.
#' @param k An integer scalar - the same k used to construct the input graph.
#' @param reduced_dims A matrix of vertices X reduced dimensions. Only required
#' if x is an \code{igraph} object. This will be automatically extracted from
#' the \code{reducedDim} slot in the \code{\linkS4class{Milo}} object.
#' @param seed An integer scalar seed to initial the pseudorandom number
#' generator
#' @details
#' This function randomly samples graph vertices, then refines them to collapse
#' down the number of neighbourhoods to be tested. The refinement behaviour can
#' be turned off by setting \code{refine=FALSE}, however, we do not recommend
#' this as neighbourhoods will contain a lot of redundancy and lead to an
#' unncecessarily larger multiple-testing burden.
#'
#' @return A \code{\linkS4class{Milo}} object containing a list of vertices and
#' the indices of vertices that constitute the neighbourhoods in the
#' neighbourhoods slot. If the input is a \code{igraph} object then the output
#' is a list of vertices and the indices of vertices that constitute the
#' neighbourhoods.
#'
#' @author
#' Mike Morgan, refined sampling written by Emma Dann
#'
#' @examples
#'
#' requires(igraph)
#' m <- matrix(rnorm(10000), ncol=10)
#' milo <- buildGraph(m, d=10)
#'
#' milo <- makeNeighbourhoods(milo, prop=0.3)
#' milo
#' @name makeNeighbourhoods
NULL

#' @export
#' @rdname makeNeighbourhoods
#' @importFrom BiocNeighbors findKNN
makeNeighbourhoods <- function(x, prop=0.1, k=21, refined=TRUE, seed=42, reduced_dims=NULL){

    if(class(x) == "Milo"){
        message("Checking valid object")
        # check that a graph has been built
        if(!.valid_graph(graph(x))){
            stop("Not a valid Milo object - graph is missing. Please run buildGraph() first.")
        }

        if(isFALSE(refined)){
            message("No sampling refinement - returning randomly selected vertices")
            neighbourhoods(x) <- .sample_vertices(graph(x), prop=prop, seed=seed, return.vertices=FALSE)
            return(x)
        } else if(isTRUE(refined)){
            message("Using refined sampling")
            random.vertices <- .sample_vertices(graph(x), prop=prop, seed=seed, return.vertices=TRUE)
            vertex.knn <- findKNN(X=reducedDim(x, "PCA"), k=k, subset=as.vector(random.vertices))
            refined.vertices <- unique(V(graph(x))[sapply(1:nrow(vertex.knn$index), function(i) refine_vertex(vertex.knn, i, reducedDim(x, "PCA")))])
            neighbourhoods(x) <- sapply(1:length(refined.vertices), FUN=function(X) neighbors(graph(x), v=refined.vertices[X]))

            return(x)
        }
    } else if(class(x) == "igraph"){
        if(is.null(reduced_dims) & isTRUE(refined)){
            stop("No reduced dimensions provided - required for refined sampling")
        } else if(isTRUE(refined)){
            random.vertices <- .sample_vertices(x, prop=prop, seed=seed, return.vertices=TRUE)
            vertex.knn <- findKNN(X=reduced_dims, k=k, subset=as.vector(random.vertices))

            # this should be embarasingly parallelisable
            refined.vertices <- V(x)[sapply(1:nrow(vertex.knn$index), function(i) refine_vertex(vertex.knn, i, reduced_dims))]
            vertex.list.refined <- sapply(1:length(refined.vertices), FUN=function(X) neighbors(x, v=refined.vertices[X]))

            return(refined.vertices)
        } else if(isFALSE(refined)){
            return(.sample_vertices(x, prop=prop, seed=seed, return.vertices=FALSE))
        }

    } else{
        stop(paste0("Data format: ", class(x), " not recognised. Should be Milo or igraph"))
    }
}


#' @importFrom igraph is_igraph
.valid_graph <- function(x){
    # check for a valid graph
    if(isTRUE(is_igraph(x))){
        TRUE
    } else{
        FALSE
    }
}

#' @import igraph
.sample_vertices <- function(graph, prop, return.vertices=FALSE, seed=42){
    set.seed(seed)
    # define a set of vertices and neihbourhood centers - extract the neihbourhoods of these cells
    random.vertices <- sample(V(graph), size=floor(prop*length(V(graph))))

    if(isTRUE(return.vertices)){
        return(random.vertices)
    } else{
        message("Finding neighbours of sampled vertices")
        vertex.list <- sapply(1:length(random.vertices), FUN=function(X) neighbors(graph, v=random.vertices[X]))
        return(vertex.list)
    }
}

#' @export
#' @rdname makeNeighbourhoods
#' @importFrom BiocNeighbors findKNN
refine_vertex <- function(vertex.knn, v.ix, X_pca){
    # vertex.knn: KNN graph for randomly sampled points (output of BiocNeighbors::findKNN)
    # v.ix: index of vertex to refine in vertex.knn

    ## Calculate median profile of KNNs of vertex
    v.med <- apply(X_pca[vertex.knn$index[v.ix,],], 2, median)
    ## Find the closest point to the median and sample
    refined.vertex <- findKNN(rbind(v.med, X_pca), subset=1, k=1)[["index"]][1] - 1 ## -1 to remove the median
    return(refined.vertex)
}
