###########################
### MILO PLOTTING UTILS ###
###########################


#' Plot histogram of neighbourhood sizes
#'
#' This function plots the histogram of the number of cells belonging to
#' each neighbourhood
#'
#' @param milo A \code{\linkS4class{Milo}} object with a non-empty \code{nhoods}
#' slot.
#' @param bins number of bins for \code{geom_histogram}
#'
#' @return A \code{ggplot-class} object
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
    stop("Not a valid Milo object - nhoods are missing. Please run makeNhoods() first.")
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
  if (is_not_empty) {
    is_igraph_vx <- class(milo@nhoods[[sample(1:n_neigh, 1)]]) == "igraph.vs"
    if (isTRUE(is_igraph_vx)){
      TRUE
    } else {
        FALSE
      }
  } else {
    FALSE
  }
}

### Plotting DA test results ###

#' Plot Milo test results on reduced dimensiona
#'
#' Visualize log-FC estimated with differential nhood abundance testing
#' on embedding of original single-cell dataset.
#'
#' @param x A \code{\linkS4class{Milo}} object
#' @param milo_results A `data.frame` containing the results of differential nhood abundance testing (output of \code{testNhoods})
#' --> this will need to be changed/removed when output of testNhoods changes
#' @param nhood_reduced_dims a character indicating the name of the \code{reducedDim} slot in the
#' \code{\linkS4class{Milo}} object to use as (default: 'UMAP').
#' @param filter_alpha the spatialFDR cutoff used as a significance threshold. If not \code{NULL} the logFC will be plotted only for
#' significantly DA nhoods (default: NULL)
#' @param split_by A character indicating the \code{colData} column in \code{x} to use for faceting
#' e.g. useful to visualize results by cell type
#' @param pt_size size of scatterplot points (default: 1.5)
#' @param components vector of reduced dimensions components to plot (default: c(1,2))
#'
#' @return a \code{ggplot-class} object
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
#' @importFrom stats filter
plotMiloReducedDim <- function(x, milo_results, nhood_reduced_dims="UMAP", filter_alpha=NULL, split_by=NULL,
                               pt_size=1.5, components=c(1,2)){
  ## Check for valid nhoodReducedDim object
  # Should have nrows = no. of nhoods + no. of cells
  if (!nhood_reduced_dims %in% names(nhoodReducedDim(x))){
    stop(paste(nhood_reduced_dims, "is not the name of an embedding in nhoodReducedDims(x). Available reductions are:", paste((nhoodReducedDim(x)), collapse = ", ")))
  } else if (nrow(nhoodReducedDim(x, nhood_reduced_dims)) != ncol(x) + length(nhoods(x))) {
    stop(paste(nhood_reduced_dims, "is not a valid nhoodReducedDims(x) object. The number of rows should match the sum of the number of cells and the number of neighbourhoods"))
    }

  ## Join test results and dimensionality reductions
  rdim_df <- data.frame(nhoodReducedDim(x, hood_reduced_dims)[,components])
  colnames(rdim_df) <- c('X','Y')

  n_nhoods <- length(nhoods(x))
  rdim_df[,"Nhood"] <- ifelse(1:nrow(rdim_df) %in% c(1:n_nhoods), c(1:n_nhoods), NA)
  viz_df  <- left_join(rdim_df, milo_results, by="Nhood")
  viz_df[is.na(viz_df["nhIndex"]),'nhIndex'] <- 1:ncol(x) # Add index also to single-cells

  if (!is.null(split_by)){
    split_df <- data.frame(split_by=colData(x)[,split_by])
    split_df[,"nhIndex"] <- 1:nrow(split_df)
    viz_df  <- left_join(viz_df, split_df, by="nhIndex")
  }

  ## Filter significant DA nhoods
  if (!is.null(filter_alpha)) {
    if (filter_alpha > 0) {
      viz_df <- mutate(viz_df, logFC = ifelse(SpatialFDR > filter_alpha, NA, logFC))
    }
  }

  ## Plot
  pl <-
    ggplot(data = arrange(viz_df, abs(logFC)),
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
    xlab(paste(nhood_reduced_dims, components[1], sep="_")) +
    ylab(paste(nhood_reduced_dims, components[2], sep="_"))

  if (!is.null(split_by)) {
    pl <- pl + facet_wrap(split_by~.)
  }
  if (!is.null(filter_alpha)) {
    pl <- pl +
      scale_color_manual(values = 'grey', label = paste("SpatialFDR >", round(filter_alpha, 2))) +
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

