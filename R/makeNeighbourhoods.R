#' Define neighbourhoods on a graph (fast)
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
#' @param d The number of dimensions to use if the input is a matrix of cells
#' X reduced dimensions. 
#' @param reduced_dims If x is an \code{\linkS4class{Milo}} object, a character indicating the name of the \code{reducedDim} slot in the 
#' \code{\linkS4class{Milo}} object to use as (default: 'PCA'). If x is an \code{igraph} object, a 
#' matrix of vertices X reduced dimensions.
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
#' Emma Dann, Mike Morgan
#'
#' @examples
#'
#' requires(igraph)
#' m <- matrix(rnorm(10000), ncol=10)
#' milo <- buildGraph(m, d=10)
#'
#' milo <- makeNeighbourhoods(milo, prop=0.1)
#' milo
#' @name makeNeighbourhoods
makeNeighbourhoods <- function(x, prop=0.1, k=21, d=30, refined=TRUE, seed=42, reduced_dims="PCA") {
    if(class(x) == "Milo"){
        message("Checking valid object")
        # check that a graph has been built
        if(!.valid_graph(miloR::graph(x))){ ## Had to specify here for conflict w igraph
            stop("Not a valid Milo object - graph is missing. Please run buildGraph() first.")
        }
        graph <- miloR::graph(x)
        X_reduced_dims  <- reducedDim(x, reduced_dims)
        if (d > ncol(X_reduced_dims)) {
            warning(paste("Warning: specified d is higher than the total number of dimensions in reducedDim(x, reduced_dims). Falling back to using", ncol(X_reduced_dims),"dimensions\n"))
            d <- ncol(X_reduced_dims)
        }
        X_reduced_dims  <- X_reduced_dims[,1:d]
    } else if(class(x) == "igraph"){
        if(!is.matrix(reduced_dims) & isTRUE(refined)){
            stop("No reduced dimensions matrix provided - required for refined sampling")
        }
        graph <- x
        X_reduced_dims <- reduced_dims
    } else{
        stop(paste0("Data format: ", class(x), " not recognised. Should be Milo or igraph"))
    }
    random_vertices <- .sample_vertices(graph, prop, seed)
    if (isFALSE(refined)) {
        sampled_vertices <- random_vertices
    } else if (isTRUE(refined)) {
        vertex.knn <-
            findKNN(
                X = X_reduced_dims,
                k = k,
                subset = as.vector(random_vertices),
                get.index = TRUE,
                get.distance = FALSE
            )
        # knn_mat <-
        #     Matrix(0, nrow = length(as.vector(random_vertices)), ncol = ncol(milo), sparse=TRUE)
        # # knn_mat <- as(knn_mat, "sparseMatrix")
        # for (ix in 1:nrow(knn_mat)) {
        #     knn_mat[ix, vertex.knn$index[ix, ]] <- 1
        # }
        # ## Calculate avg profile of nearest neighbors
        # nh_reduced_dims <- knn_mat %*% X_reduced_dims
        # nh_reduced_dims <- t(apply(
        #     nh_reduced_dims,
        #     1,
        #     FUN = function(x)
        #         x / k
        # ))
       
        nh_reduced_dims <- t(apply(vertex.knn$index, 1, function(x) colMedians(X_reduced_dims[x,])))
        colnames(nh_reduced_dims) <- colnames(X_reduced_dims)
        rownames(nh_reduced_dims) <- paste0('nh_', 1:nrow(nh_reduced_dims))
        
        ## Search nearest cell to average profile
        # I have to do this trick because as far as I know there is no fast function to 
        # search for NN between 2 distinct sets of points (here I'd like to search NNs of 
        # nh_reduced_dims points among X_reduced_dims points). Suggestions are welcome
        all_reduced_dims <- rbind(nh_reduced_dims, X_reduced_dims)
        nn_mat <- BiocNeighbors::findKNN(all_reduced_dims,
                                         k = nrow(nh_reduced_dims) + 1,
                                         subset = rownames(nh_reduced_dims))[["index"]]
        nn_mat_names <- apply(nn_mat, c(1,2), function(x) rownames(all_reduced_dims)[x])
        i = 1
        sampled_vertices <- rep("nh_0", nrow(nn_mat))
        while (any(grepl(sampled_vertices, pattern = "nh_"))) {
            update_ix <- grep(sampled_vertices, pattern="nh_")
            sampled_vertices[update_ix] <- nn_mat_names[update_ix, i]
            i <- i + 1
        }
        sampled_vertices <- match( sampled_vertices, rownames(X_reduced_dims))
    } 
    
    sampled_vertices <- unique(sampled_vertices)
    nh_list <-
        sapply(
            1:length(sampled_vertices),
            FUN = function(X)
                neighbors(graph, v = sampled_vertices[X])
        )
    nh_list <- setNames(nh_list, sampled_vertices)
    if(class(x) == "Milo"){ 
        neighbourhoods(x) <- nh_list 
        return(x)
    } else {
        return(nh_list)
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
.sample_vertices <- function(graph, prop, seed=42){
    set.seed(seed)
    # define a set of vertices and neihbourhood centers - extract the neihbourhoods of these cells
    random.vertices <- sample(V(graph), size=floor(prop*length(V(graph))))
    return(random.vertices)
}

### PLotting utility function ###

plotNeighborhoodSizeHist <- function(milo, bins=50){
    if (! isTRUE(.valid_neighbourhood(milo))){
        stop("Not a valid Milo object - neighbourhoods are missing. Please run makeNeighbourhoods() first.")
    }
    df <- data.frame(nh_size=sapply(milo@neighbourhoods, function(x) length(x))) 
    ggplot(data=df, aes(nh_size)) + geom_histogram(bins=bins) +
        xlab("Neighbourhood size") +
        theme_classic(base_size = 16)
}


#' @importFrom igraph is_igraph
.valid_neighbourhood <- function(milo){
    # check for a valid neighbourhood slot
    n_neigh <- length(milo@neighbourhoods)
    is_not_empty <- n_neigh > 0
    is_igraph_vx <- class(milo@neighbourhoods[[sample(1:n_neigh, 1)]]) == "igraph.vs" 
    if (isTRUE(is_igraph_vx & is_not_empty)){
        TRUE
    } else {
        FALSE
    }
}


