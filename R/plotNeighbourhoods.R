###########################
### MILO PLOTTING UTILS ###
###########################

### Plotting neighbourhoods stats ###

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

### Plotting DA test results ###

#' Plot Milo test results on reduced dimensiona
#' 
#' Visualize log-FC estimated with differential neighbourhood abundance testing
#' on embedding of original single-cell dataset. 
#'
#' @param x A \code{\linkS4class{Milo}} object 
#' @param milo_results A `data.frame` containing the results of differential neighbourhood abundance testing (output of \code{testNeighborhoods})
#' --> this will need to be changed/removed when output of testNeighbourhoods changes
#' @param reduced_dims a character indicating the name of the \code{reducedDim} slot in the 
#' \code{\linkS4class{Milo}} object to use as (default: 'UMAP').
#' @param filter_alpha the spatialFDR cutoff used as a significance threshold. If not \code{NULL} the logFC will be plotted only for 
#' significantly DA neighbourhoods (default: NULL)
#' @param pt_size size of scatterplot points (default: 1.5)
#' @param components vector of reduced dimensions components to plot (default: c(1,2))
#' 
#' @return a \code{\linkS4class{ggplot}} object
#' 
#' @author Emma Dann
#' 
#' @examples 
#' NULL
#' 
#' @export
#' @rdname plotMiloReducedDim
#' @import ggplot2
#' @importFrom dplyr left_join mutate arrange
plotMiloReducedDim <- function(x, milo_results, reduced_dims="UMAP", filter_alpha=NULL, 
                               pt_size=1.5, components=c(1,2)
){
  ## Join test results and dimensionality reductions
  if (reduced_dims %in% reducedDimNames(x)){
    rdim_df <- data.frame(reducedDim(x, reduced_dims))[,components]
  } else {
    stop(paste(reduced_dims, "is not the name of a dimensionality reduction result in x. Available reductions are:", paste(reducedDimNames(x), collapse = ", ")))
  }
  colnames(rdim_df) <- c('X','Y')
  rdim_df[,"nhIndex"] <- 1:nrow(rdim_df)
  nhIndex <- unlist(neighbourhoodIndex(x))
  milo_results[,"nhIndex"] <- nhIndex
  viz2_df  <- left_join(rdim_df, milo_results, by="nhIndex") 
  
  ## Filter significant DA neighbourhoods 
  if (!is.null(filter_alpha)) {
    viz2_df <- mutate(viz2_df, logFC = ifelse(SpatialFDR > filter_alpha, NA, logFC)) 
  }
  
  ## Plot 
  pl <-
    ggplot(data = arrange(viz2_df, abs(logFC)),
           aes(X, Y)) +
    geom_point(aes(color = ''), size = pt_size / 3, alpha = 0.5) +
    geom_point(
      data = . %>% filter(!is.na(SpatialFDR)),
      aes(fill = logFC),
      size = pt_size,
      stroke = 0.1,
      # colour="black",
      shape = 21
    ) +
    scale_fill_gradient2(
      midpoint = 0,
      high = "red",
      low = "blue",
      name = "log-FC"
    ) +
    xlab(paste(reduced_dims, components[1], sep="_")) +
    ylab(paste(reduced_dims, components[2], sep="_"))
  
  if (!is.null(filter_alpha)) {
    pl <- pl +
      scale_color_manual(values = 'grey', label = paste("SpatialFDR <", round(filter_alpha, 2))) +
      guides(colour = guide_legend(
        '',
        override.aes = list(
          shape = 21,
          colour = "black",
          fill = "grey50",
          size = pt_size,
          alpha = 1,
          stroke = 0.1
        )
      ))
  } else {
    pl <- pl + 
      scale_color_manual(values = 'grey') +
      guides(color="none")
  }
  
  pl <- pl +
    theme_classic(base_size = 16) +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  return(pl)
}

