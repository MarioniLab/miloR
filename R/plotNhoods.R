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
#' require(igraph)
#' require(SingleCellExperiment)
#' ux.1 <- matrix(rpois(12000, 5), ncol=400)
#' ux.2 <- matrix(rpois(12000, 4), ncol=400)
#' ux <- rbind(ux.1, ux.2)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#'
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'                             reducedDims=SimpleList(PCA=pca$x))
#' colnames(sce) <- paste0("Cell", seq_len(ncol(sce)))
#' milo <- Milo(sce)
#' milo <- buildGraph(milo, k=20, d=10, transposed=TRUE)
#'
#' milo <- makeNhoods(milo, d=10, prop=0.1)
#' plotNhoodSizeHist(milo)
#'

#' @export
#' @rdname plotNhoodSizeHist
#' @importFrom ggplot2 ggplot geom_histogram xlab theme_classic
#' @importFrom igraph neighbors
#' @importFrom grDevices colorRampPalette
plotNhoodSizeHist <- function(milo, bins=50){
  if (! isTRUE(.valid_nhood(milo))){
    stop("Not a valid Milo object - nhoods are missing. Please run makeNhoods() first.")
  }
  df <- data.frame(nh_size=colSums(nhoods(milo)))

  ggplot(data=df, aes(nh_size)) + geom_histogram(bins=bins) +
    xlab("Neighbourhood size") +
    theme_classic(base_size = 16)
}


#' @importFrom igraph is_igraph
.valid_nhood <- function(milo){
  # check for a valid nhood slot
  n_neigh <- ncol(nhoods(milo))
  is_not_empty <- n_neigh > 0
  if (is_not_empty) {
    # is_graph_vx <- is(milo@nhoods[[sample(seq_len(n_neigh), 1)]], "igraph.vs")
    # is_numeric_vc <- is(milo@nhoods[[sample(seq_len(n_neigh), 1)]], "numeric")
    # if (isTRUE(is_igraph_vx) | isTRUE(is_numeric_vc)){
      TRUE
    # } else {
    #     FALSE
    #   }
  } else {
    FALSE
  }
}


### Plotting neighbourhood graph ###

#' Plot graph of neighbourhood
#'
#' Visualize graph of neighbourhoods
#'
#' @param x A \code{\linkS4class{Milo}} object
#' @param layout this can be (a) a character indicating the name of the \code{reducedDim} slot in the
#' \code{\linkS4class{Milo}} object to use for layout (default: 'UMAP') (b) an igraph layout object
#' @param colour_by this can be a data.frame of milo results or a character corresponding to a column in colData
#' @param subset.nhoods A logical, integer or character vector indicating a subset of nhoods to show in plot
#' (default: NULL, no subsetting)
#' @param size_range a numeric vector indicating the range of node sizes to use for plotting (to avoid overplotting
#' in the graph)
#' @param node_stroke a numeric indicating the desired thickness of the border around each node
#' @param ... arguments to pass to \code{ggraph}
#'
#' @return a \code{ggplot-class} object
#'
#' @author Emma Dann
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotNhoodGraph
#' @import igraph
#' @import ggraph
#' @importFrom SummarizedExperiment colData<-
#' @importFrom RColorBrewer brewer.pal
plotNhoodGraph <- function(x, layout="UMAP", colour_by=NA, subset.nhoods=NULL, size_range=c(0.5,3),
                           node_stroke= 0.3, ... ){
  ## Check for valid nhoodGraph object
  if(!.valid_graph(nhoodGraph(x))){
    stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
  }
  if (is.character(layout)) {
    if (!layout %in% names(reducedDims(x))) {
      stop(layout, "isn't in readucedDim(x) - choose a different layout")
    }
  }
  nh_graph <- nhoodGraph(x)

  ## Subset
  if (!is.null(subset.nhoods)) {
    nh_graph <- igraph::induced_subgraph(nh_graph, vids = which(as.numeric(V(nh_graph)$name) %in% unlist(nhoodIndex(x)[subset.nhoods])))
  }


  ## Order vertex ids by size (so big nhoods are plotted first)
  nh_graph <- permute(nh_graph, order(vertex_attr(nh_graph)$size, decreasing=TRUE))

  ## Define layout
  if (is.character(layout)) {
    redDim <- layout
    layout <- reducedDim(x, redDim)[as.numeric(vertex_attr(nh_graph)$name),]
    # make sure this is a matrix!
    if(!any(class(layout) %in% c("matrix"))){
        warning("Coercing layout to matrix format")
        layout <- as(layout, "matrix")
    }
  }

  ## Define node color
  if (!is.na(colour_by)) {
    if (colour_by %in% colnames(colData(x))) {

      col_vals <- colData(x)[as.numeric(vertex_attr(nh_graph)$name), colour_by]
      if (!is.numeric(col_vals)) {
        col_vals <- as.character(col_vals)
        }
      V(nh_graph)$colour_by <- col_vals
    } else {
      stop(colour_by, "is not a column in colData(x)")
    }
  } else {
    V(nh_graph)$colour_by <- V(nh_graph)$size
    colour_by <- "Nhood size"
  }

  if(colour_by %in% c("logFC")){
    plot.g <- simplify(nh_graph)

    pl <- ggraph(simplify(nh_graph), layout = layout) +
      geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
      geom_node_point(aes(fill = colour_by, size = size), shape=21, stroke=node_stroke) +
      scale_size(range =size_range, name="Nhood size") +
      scale_edge_width(range = c(0.2,3), name="overlap size") +
      theme_classic(base_size=14) +
      theme(axis.line = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), axis.title = element_blank())
    # theme_graph()

  } else{
    pl <- ggraph(simplify(nh_graph), layout = layout) +
      geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
      geom_node_point(aes(fill = colour_by, size = size), shape=21, stroke=node_stroke) +
      scale_size(range = size_range, name="Nhood size") +
      scale_edge_width(range = c(0.2,3), name="overlap size") +
      theme_classic(base_size=14) +
      theme(axis.line = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), axis.title = element_blank())
      # theme_graph()
  }

  if (is.numeric(V(nh_graph)$colour_by)) {
    pl <- pl + scale_fill_gradient2(name=colour_by)
  } else {
    mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(V(nh_graph)$colour_by)))
    pl <- pl + scale_fill_manual(values=mycolors, name=colour_by, na.value="white")
  }
  pl
  }

#' Plot Milo results on graph of neighbourhood
#'
#' Visualize log-FC estimated with differential nhood abundance testing
#' on embedding of original single-cell dataset.
#'
#' @param x A \code{\linkS4class{Milo}} object
#' @param milo_res a data.frame of milo results
#' @param alpha significance level for Spatial FDR (default: 0.05)
#' @param res_column which column of \code{milo_res} object to use for color (default: logFC)
#' @param ... arguments to pass to \code{plotNhoodGraph}
#'
#' @return a \code{ggplot} object
#'
#' @author Emma Dann
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotNhoodGraphDA
#' @import igraph
plotNhoodGraphDA <- function(x, milo_res, alpha=0.05, res_column = "logFC", ... ){
  if(!.valid_graph(nhoodGraph(x))){
    stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
  }
  if (is.character(layout)) {
    if (!layout %in% names(reducedDims(x))) {
      stop(layout, "is not in readucedDim(x) - choose a different layout")
    }
  }

  ## Add milo results to colData
  signif_res <- milo_res
  signif_res[signif_res$SpatialFDR > alpha,res_column] <- 0
  colData(x)[res_column] <- NA
  colData(x)[unlist(nhoodIndex(x)[signif_res$Nhood]),res_column] <- signif_res[,res_column]

  ## Plot logFC
  plotNhoodGraph(x, colour_by = res_column, ... )
}

#' Plot graph of neighbourhoods coloring by nhoodGroups
#'
#' Visualize grouping of neighbourhoods obtained with \code{groupNhoods}
#'
#' @param x A \code{\linkS4class{Milo}} object
#' @param milo_res a data.frame of milo results containing the \code{nhoodGroup} column
#' @param show_groups a character vector indicating which groups to plot
#' all other neighbourhoods will be gray
#' @param ... arguments to pass to \code{plotNhoodGraph}
#'
#' @return a \code{ggplot} object
#'
#' @author Emma Dann
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotNhoodGraphDA
#' @import igraph
plotNhoodGroups <- function(x, milo_res, show_groups=NULL, ... ){
  if(!.valid_graph(nhoodGraph(x))){
    stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
  }
  if (is.character(layout)) {
    if (!layout %in% names(reducedDims(x))) {
      stop(layout, "is not in reducedDim(x) - choose a different layout")
    }
  }

  if (!"NhoodGroup" %in% colnames(milo_res)) {
    stop("'NhoodGroup' columns is missing from milo_res. Please run groupNhoods() or define neighbourhood groupings otherwise.")
  }

  ## Add groups to colData
  # signif_res <- milo_res
  # signif_res[signif_res$SpatialFDR > alpha,res_column] <- 0
  if (!is.null(show_groups)) {
    plot_groups <- show_groups
  } else {
    plot_groups <- unique(milo_res$NhoodGroup)
  }
  colData(x)["NhoodGroup"] <- NA
  groups_res <- milo_res[milo_res$NhoodGroup %in% plot_groups,]
  colData(x)[unlist(nhoodIndex(x)[groups_res$Nhood]),"NhoodGroup"] <- groups_res$NhoodGroup

  ## Plot logFC
  plotNhoodGraph(x, colour_by = "NhoodGroup", ... )
}

#' Visualize gene expression in neighbourhoods
#'
#' Plots the average gene expression in neighbourhoods, sorted by DA fold-change
#'
#' @param x A \code{\linkS4class{Milo}} object
#' @param da.res a data.frame of DA testing results
#' @param features a character vector of features to plot (they must be in rownames(x))
#' @param alpha significance level for Spatial FDR (default: 0.1)
#' @param subset.nhoods A logical, integer or character vector indicating a subset of nhoods to show in plot
#' (default: NULL, no subsetting)
#' @param cluster_features logical indicating whether features should be clustered with hierarchical clustering.
#' If FALSE then the order in \code{features} is maintained (default: FALSE)
#' @param assay A character scalar that describes the assay slot to use for calculating neighbourhood expression.
#' (default: logcounts)
#' Of note: neighbourhood expression will be computed only if the requested features are not in the \code{nhoodExpression} slot
#' of the milo object. If you wish to plot average neighbourhood expression from a different assay, you should run
#' \code{calcNhoodExpression(x)} with the desired assay.
#' @param scale_to_1 A logical scalar to re-scale gene expression values between 0 and 1 for visualisation.
#' @param show_rownames A logical scalar whether to plot rownames or not. Generally useful to set this to
#' \code{show_rownames=FALSE} when plotting many genes.
#' @param highlight_features A character vector of feature names that should be highlighted on the right side of
#' the heatmap. Generally useful in conjunction to \code{show_rownames=FALSE}, if you are interested in only a few
#' features
#' @return a \code{ggplot} object
#'
#' @author Emma Dann
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotNhoodExpressionDA
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr mutate left_join filter percent_rank first group_by summarise
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @importFrom stats hclust
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_replace
plotNhoodExpressionDA <- function(x, da.res, features, alpha=0.1,
                                  subset.nhoods=NULL, cluster_features=FALSE, assay="logcounts",
                                  scale_to_1 = FALSE,
                                  show_rownames=TRUE,
                                  highlight_features = NULL){
  if (length(features) <= 0 | is.null(features)) {
    stop("features is empty")
  }
  ## Check if features are in rownames(x)
  if (!all(features %in% rownames(x))) {
    stop("Some features are not in rownames(x)")
  }
  ## Check if nhood expression exists
  if (dim(nhoodExpression(x))[2] == 1){
    warning("Nothing in nhoodExpression(x): computing for requested features...")
    x <- calcNhoodExpression(x, assay = assay, subset.row = features)
  }
  ## Check if all features are in nhoodExpression
  if (!all(features %in% rownames(nhoodExpression(x)))) {
    warning("Not all features in nhoodExpression(x): recomputing for requested features...")
    x <- calcNhoodExpression(x, assay = assay, subset.row = features)
  }

  expr_mat <- nhoodExpression(x)[features, ]
  colnames(expr_mat) <- seq_len(ncol(nhoods(x)))

  ## Get nhood expression matrix
  if (!is.null(subset.nhoods)) {
      expr_mat <- expr_mat[,subset.nhoods, drop=FALSE]
  }

  if (!isFALSE(scale_to_1)) {
      expr_mat <- t(apply(expr_mat, 1, function(X) (X - min(X))/(max(X)- min(X))))
      # force NAs to 0?
      if(sum(is.na(expr_mat)) > 0){
          warning("NA values found - resetting to 0")
          expr_mat[is.na(expr_mat)] <- 0
      }
  }

  rownames(expr_mat) <- sub(pattern = "-", replacement = ".", rownames(expr_mat)) ## To avoid problems when converting to data.frame

  pl_df <- data.frame(t(expr_mat)) %>%
    rownames_to_column("Nhood") %>%
    mutate(Nhood=as.double(Nhood)) %>%
    left_join(da.res, by="Nhood") %>%
    mutate(logFC_rank=percent_rank(logFC))

  ## Top plot: nhoods ranked by DA log FC
  pl_top <- pl_df %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, paste0("SpatialFDR < ", alpha), NA)) %>%
    ggplot(aes(logFC_rank, logFC)) +
    geom_hline(yintercept = 0, linetype=2) +
    geom_point(size=0.2, color="grey") +
    geom_point(data=.%>% filter(!is.na(is_signif)), aes(color=is_signif), size=1) +
    theme_bw(base_size=16) +
    ylab("DA logFC") +
    scale_color_manual(values="red", name="") +
    scale_x_continuous(expand = c(0.01, 0)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

  ## Bottom plot: gene expression heatmap
  if (isTRUE(cluster_features)) {
      row.order <- hclust(dist(expr_mat))$order # clustering
      ordered_features <- rownames(expr_mat)[row.order]
  } else {
    ordered_features <- rownames(expr_mat)
  }

  # this code assumes that colnames do not begin with numeric values
  # add 'X' to feature names with numeric first characters
  rownames(expr_mat) <- str_replace(rownames(expr_mat), pattern="(^[0-9]+)", replacement="X\\1")

  pl_df <- pl_df %>%
    pivot_longer(cols=rownames(expr_mat), names_to='feature', values_to="avg_expr") %>%
    mutate(feature=factor(feature, levels=ordered_features))

  if (!is.null(highlight_features)) {
    if (!all(highlight_features %in% pl_df$feature)){
      missing <- highlight_features[which(!highlight_features %in% pl_df$feature)]
      warning('Some elements of highlight_features are not in features and will not be highlighted. \nMissing features: ', paste(missing, collapse = ', ') )
    }
    pl_df <- pl_df %>%
      mutate(label=ifelse(feature %in% highlight_features, as.character(feature), NA))
  }

  pl_bottom <- pl_df %>%
    ggplot(aes(logFC_rank, feature, fill=avg_expr)) +
    geom_tile() +
    scale_fill_viridis_c(option="magma", name="Avg.Expr.") +
    xlab("Neighbourhoods") + ylab("Features") +
    scale_x_continuous(expand = c(0.01, 0)) +
    theme_classic(base_size = 16) +
    coord_cartesian(clip="off") +
    theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.spacing = margin(2, 2, 2, 2, "cm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(10,10,10,10)
          )

  if (!is.null(highlight_features)) {
    pl_bottom <- pl_bottom +
      geom_text_repel(data=. %>%
                        filter(!is.na(label)) %>%
                        group_by(label) %>%
                        summarise(logFC_rank=max(logFC_rank), avg_expr=mean(avg_expr), feature=first(feature)),
                      aes(label=label, x=logFC_rank),
                      size=4,
                      xlim = c(max(pl_df$logFC_rank) + 0.01, max(pl_df$logFC_rank) + 0.02),
                      min.segment.length = 0,
                      max.overlaps=Inf,
                      seed=42
      )

  }

  if(isFALSE(show_rownames)){
      pl_bottom <- pl_bottom +
          theme(axis.text.y=element_blank())
  }

  ## Assemble plot
  (pl_top / pl_bottom) +
    plot_layout(heights = c(1,4), guides = "collect") &
    theme(legend.justification=c(0, 1),
          legend.margin = margin(0,0,0,50))
}


#' Visualize gene expression in neighbourhood groups
#'
#' Plots the average gene expression in neighbourhood groups
#'
#' @param x A \code{\linkS4class{Milo}} object
#' @param da.res a data.frame of DA testing results
#' @param features a character vector of features to plot (they must be in rownames(x))
#' @param subset.nhoods A logical, integer or character vector indicating a subset of nhoods to show in plot
#' (default: NULL, no subsetting)
#' @param cluster_features logical indicating whether features should be clustered with hierarchical clustering.
#' If FALSE then the order in \code{features} is maintained (default: FALSE)
#' @param assay A character scalar that describes the assay slot to use for calculating neighbourhood expression.
#' (default: logcounts)
#' Of note: neighbourhood expression will be computed only if the requested features are not in the \code{nhoodExpression} slot
#' of the milo object. If you wish to plot average neighbourhood expression from a different assay, you should run
#' \code{calcNhoodExpression(x)} with the desired assay.
#' @param scale_to_1 A logical scalar to re-scale gene expression values between 0 and 1 for visualisation.
#' @param show_rownames A logical scalar whether to plot rownames or not. Generally useful to set this to
#' \code{show_rownames=FALSE} when plotting many genes.
#' @param highlight_features A character vector of feature names that should be highlighted on the right side of
#' the heatmap. Generally useful in conjunction to \code{show_rownames=FALSE}, if you are interested in only a few
#' features
#' @param grid.space a character setting the \code{space} parameter for \code{facet.grid} (\code{'fixed'} for equally sized facets,
#' \code{'free'} to adapt the size of facent to number of neighbourhoods in group)
#'
#' @return a \code{ggplot} object
#'
#' @author Emma Dann
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotNhoodExpressionDA
#' @import ggplot2
#' @importFrom dplyr mutate left_join filter first group_by summarise ungroup
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @importFrom stats hclust
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_replace
plotNhoodExpressionGroups <- function(x, da.res, features, alpha=0.1,
                                      subset.nhoods=NULL, cluster_features=FALSE, assay="logcounts",
                                      scale_to_1 = FALSE,
                                      show_rownames=TRUE,
                                      highlight_features = NULL,
                                      grid.space="free"){
  if (length(features) <= 0 | is.null(features)) {
    stop("features is empty")
  }
  ## Check if features are in rownames(x)
  if (!all(features %in% rownames(x))) {
    stop("Some features are not in rownames(x)")
  }
  ## Check if nhood expression exists
  if (dim(nhoodExpression(x))[2] == 1){
    warning("Nothing in nhoodExpression(x): computing for requested features...")
    x <- calcNhoodExpression(x, assay = assay, subset.row = features)
  }
  ## Check if all features are in nhoodExpression
  if (!all(features %in% rownames(nhoodExpression(x)))) {
    warning("Not all features in nhoodExpression(x): recomputing for requested features...")
    x <- calcNhoodExpression(x, assay = assay, subset.row = features)
  }

  expr_mat <- nhoodExpression(x)[features, ]
  colnames(expr_mat) <- seq_len(ncol(nhoods(x)))

  ## Get nhood expression matrix
  if (!is.null(subset.nhoods)) {
    expr_mat <- expr_mat[,subset.nhoods, drop=FALSE]
  }

  if (!isFALSE(scale_to_1)) {
    expr_mat <- t(apply(expr_mat, 1, function(X) (X - min(X))/(max(X)- min(X))))
    # force NAs to 0?
    if(sum(is.na(expr_mat)) > 0){
      warning("NA values found - resetting to 0")
      expr_mat[is.na(expr_mat)] <- 0
    }
  }

  rownames(expr_mat) <- sub(pattern = "-", replacement = ".", rownames(expr_mat)) ## To avoid problems when converting to data.frame

  pl_df <- data.frame(t(expr_mat)) %>%
    rownames_to_column("Nhood") %>%
    mutate(Nhood=as.double(Nhood)) %>%
    left_join(da.res, by="Nhood") %>%
    group_by(NhoodGroup) %>%
    mutate(logFC_rank=rank(logFC, ties.method="random")) %>%
    ungroup()

  ## plot: gene expression heatmap
  if (isTRUE(cluster_features)) {
    row.order <- hclust(dist(expr_mat))$order # clustering
    ordered_features <- rownames(expr_mat)[row.order]
  } else {
    ordered_features <- rownames(expr_mat)
  }

  # this code assumes that colnames do not begin with numeric values
  # add 'X' to feature names with numeric first characters
  rownames(expr_mat) <- str_replace(rownames(expr_mat), pattern="(^[0-9]+)", replacement="X\\1")

  pl_df <- pl_df %>%
    pivot_longer(cols=rownames(expr_mat), names_to='feature', values_to="avg_expr") %>%
    mutate(feature=factor(feature, levels=ordered_features))

  if (!is.null(highlight_features)) {
    if (!all(highlight_features %in% pl_df$feature)){
      missing <- highlight_features[which(!highlight_features %in% pl_df$feature)]
      warning('Some elements of highlight_features are not in features and will not be highlighted. \nMissing features: ', paste(missing, collapse = ', ') )
    }
    pl_df <- pl_df %>%
      mutate(label=ifelse(feature %in% highlight_features, as.character(feature), NA))
  }

  pl_bottom <- pl_df %>%
    ggplot(aes(logFC_rank, feature, fill=avg_expr)) +
    geom_tile() +
    scale_fill_viridis_c(option="magma", name="Avg.Expr.") +
    xlab("Neighbourhoods") + ylab("Features") +
    scale_x_continuous(expand = c(0.01, 0)) +
    theme_classic(base_size = 16) +
    coord_cartesian(clip="off") +
    facet_grid(.~NhoodGroup, labeller = "label_both", scales="free", space=grid.space) +
    theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          # panel.spacing = margin(2, 2, 2, 2, "cm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(10,10,10,10)
    )

  if (!is.null(highlight_features)) {
    pl_bottom <- pl_bottom +
      geom_text_repel(data=. %>%
                        filter(!is.na(label)) %>%
                        group_by(label) %>%
                        summarise(logFC_rank=max(logFC_rank), avg_expr=mean(avg_expr), feature=first(feature)),
                      aes(label=label, x=logFC_rank),
                      size=4,
                      xlim = c(max(pl_df$logFC_rank) + 0.01, max(pl_df$logFC_rank) + 0.02),
                      min.segment.length = 0,
                      seed=42
      )

  }

  if(isFALSE(show_rownames)){
    pl_bottom <- pl_bottom +
      theme(axis.text.y=element_blank())
  }

  pl_bottom
}


#' Visualize DA results as a beeswarm plot
#'
#' @param da.res a data.frame of DA testing results
#' @param group.by a character scalar determining which column of \code{da.res} to use for grouping.
#' This can be a column added to the DA testing results using the `annotateNhoods` function.
#' If \code{da.res[,group.by]} is a character or a numeric, the function will coerce it to a factor (see details)
#' (default: NULL, no grouping)
#' @param alpha significance level for Spatial FDR (default: 0.1)
#' @param subset.nhoods A logical, integer or character vector indicating a subset of nhoods to show in plot
#' (default: NULL, no subsetting)
#'
#' @return a \code{ggplot} object
#'
#' @details The group.by variable will be coerced to a factor. If you want the variables in group.by to be
#' in a given order make sure you set the column to a factor with the levels in the right order before running the
#' function.
#'
#' @author Emma Dann
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotDAbeeswarm
#' @import ggplot2
#' @importFrom dplyr mutate filter arrange
#' @importFrom ggbeeswarm geom_quasirandom
plotDAbeeswarm <- function(da.res, group.by=NULL, alpha=0.1, subset.nhoods=NULL){
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", group.by,")?")
    }
    if (is.numeric(da.res[,group.by])) {
      stop(group.by, " is a numeric variable. Please bin to use for grouping.")
    }
    da.res <- mutate(da.res, group_by = da.res[,group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }

  if (!is.factor(da.res[,"group_by"])) {
    message("Converting group.by to factor...")
    da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))
    # anno_vec <- factor(anno_vec, levels=unique(anno_vec))
  }

  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods,]
  }

  da.res %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    ggplot(aes(group_by, logFC, color=logFC_color)) +
    scale_color_gradient2() +
    guides(color="none") +
    xlab(group.by) + ylab("Log Fold Change") +
    geom_quasirandom(alpha=1) +
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))

}


#' Visualize DA results as an MAplot
#'
#' @param da.res A data.frame of DA testing results
#' @param null.mean A numeric scalar determining the expected value of the log fold change under the null
#' hypothesis. \code{default=0}.
#' @param alpha A numeric scalar that represents the Spatial FDR threshold for statistical significance.
#'
#' @return a \code{ggplot} object
#'
#' @details MA plots provide a useful means to evaluate the distribution of log fold changes after differential
#' abundance testing. In particular, they can be used to diagnose global shifts that occur in the presence of
#' confounding between the number of cells acquired and the experimental variable of interest. The expected null
#' value for the log FC distribution (grey dashed line), along with the mean observed log fold change for non-DA
#' neighbourhoods (purple dashed line) are plotted for reference. The deviation between these two lines can give
#' an indication of biases in the results, such as in the presence of a single strong region of DA leading to an
#' increase in false positive DA neighbourhoods in the opposite direction.
#'
#' @author Mike Morgan
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotNhoodMA
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
plotNhoodMA <- function(da.res, alpha=0.05, null.mean=0){
  if(isFALSE(any(class(da.res) %in% c("tibble", "data.frame")))){
    stop("Input `da.res` not recognised - please input the results data.frame from DA testing.")
  }

  if(!is.numeric(alpha)){
    stop("alpha is not numeric - please provide a numeric value for the Spatial FDR threshold")
  }

  if(!is.numeric(null.mean)){
    stop("typeof null.mean not recoginised, please provide a numeric value")
  }

  if(isFALSE(!all(colnames(da.res) %in% c("SpatialFDR", "logCPM", "logFC")))){
    stop("Please provide a results data.frame containing the logCPM, logFC and SpatialFDR results")
  }

  sig.cols <- c("black", "red")
  names(sig.cols) <- c(FALSE, TRUE)
  max.lfc <- max(abs(da.res$logFC))
  max.eps <- max.lfc * 0.1

  emp.null <- mean(da.res$logFC[da.res$SpatialFDR >= alpha])
  min.x <- min(da.res$logCPM)
  minx.eps <- min.x * 0.01
  max.x <- max(da.res$logCPM)
  maxx.eps <- max.x * 0.01

  da.res$Sig <- da.res$SpatialFDR < alpha
  ma.p <- ggplot(da.res, aes(x=logCPM, y=logFC, colour=Sig)) +
    geom_hline(yintercept=null.mean, lty=2, colour='grey80') +
    geom_hline(yintercept=emp.null, lty=2, colour='purple') +
    geom_point() +
    theme_cowplot() +
    scale_colour_manual(values=sig.cols) +
    scale_y_continuous(limits=c(-max.lfc - max.eps, max.lfc + max.eps)) +
    scale_x_continuous(limits=c(min.x-minx.eps, max.x + maxx.eps)) +
    labs(x="Mean log scaled counts", y="Log fold-change") +
    guides(colour=guide_legend(title=paste0("SpatialFDR < ", alpha))) +
    NULL

  return(ma.p)
}


#' Plot the number of cells in a neighbourhood per sample and condition
#'
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{nhoodCounts}
#' slot.
#' @param nhoods A string vector specifying the IDs of the neighbourhoods to plot.
#' These should correspond to row names in \code{nhoodCounts(milo)}
#' @param design.df A \code{data.frame} which matches samples to a condition of interest.
#' The row names should correspond to the samples. You can use the same \code{design.df}
#' that you already used in the \code{testNhoods} function.
#' @param condition String specifying the condition of interest Has to be a column in the \code{design}.
#' @param n_col Number of columns in the output \code{ggplot}.
#'
#' @return A \code{ggplot-class} object
#'
#' @author
#' Nick Hirschmüller
#'
#' @examples
#'
#' require(SingleCellExperiment)
#' ux.1 <- matrix(rpois(12000, 5), ncol=300)
#' ux.2 <- matrix(rpois(12000, 4), ncol=300)
#' ux <- rbind(ux.1, ux.2)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#'
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'                             reducedDims=SimpleList(PCA=pca$x))
#' milo <- Milo(sce)
#' milo <- buildGraph(milo, k=20, d=10, transposed=TRUE)
#' milo <- makeNhoods(milo, k=20, d=10, prop=0.3)
#' milo <- calcNhoodDistance(milo, d=10)
#'
#' cond <- sample(c("A","B","C"),300,replace=T)
#'
#' meta.df <- data.frame(Condition=cond, Replicate=c(rep("R1", 100), rep("R2", 100), rep("R3", 100)))
#' meta.df$SampID <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
#' milo <- countCells(milo, meta.data=meta.df, samples="SampID")
#'
#' design.mtx <- data.frame("Condition"=c(rep("A", 3), rep("B", 3), rep("C",3)),
#'                          "Replicate"=rep(c("R1", "R2", "R3"), 3))
#' design.mtx$SampID <- paste(design.mtx$Condition, design.mtx$Replicate, sep="_")
#' rownames(design.mtx) <- design.mtx$SampID
#'
#' plotNhoodCounts(x = milo,
#'                 nhoods = c(1,2),
#'                 design.df = design.mtx,
#'                 condition = "Condition")
#'
#' @export
#' @rdname plotNhoodCounts
#' @importFrom ggplot2 ggplot geom_point stat_summary facet_wrap ylab
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
plotNhoodCounts <- function(x, nhoods, design.df, condition, n_col=3){
  if (!is(x, "Milo")) {
    stop("Unrecognised input type - must be of class Milo")
  }
  if (ncol(nhoodCounts(x)) == 1 & nrow(nhoodCounts(x)) == 1) {
    stop("Neighbourhood counts missing - please run countCells() first")
  }
  if (ncol(nhoods(x)) == 1 & nrow(nhoods(x)) == 1) {
    stop("No neighbourhoods found. Please run makeNhoods() first.")
  }
  if (!all(nhoods %in% rownames(nhoodCounts(x)))) {
    stop("Specified neighbourhoods do not exist - these should correspond to row names in nhoodCounts(x)")
  }
  if (!is(design.df,"data.frame")){
    stop("The design.df has to be of type data.frame")
  }
  if (!condition %in% colnames(design.df)){
    stop("Condition of interest has to be a column in the design matrix")
  }


  nhood.counts.df <- data.frame(as.matrix(nhoodCounts(x)[nhoods, , drop=F]))
  nhood.counts.df <- rownames_to_column(nhood.counts.df, "nhoods.id")
  nhood.counts.df.long <- pivot_longer(nhood.counts.df, cols=2:ncol(nhood.counts.df),
                                       names_to = "experiment",
                                       values_to = "values")

  tmp.desgin <- rownames_to_column(design.df, "experiment")[,c("experiment",condition)]
  nhood.counts.df.long <- left_join(nhood.counts.df.long,
                                    tmp.desgin,
                                    by="experiment")
  nhood.counts.df.long$nhoods.id <- paste("Nhood:", nhood.counts.df.long$nhoods.id)



  p <- ggplot(nhood.counts.df.long, aes_string(x=condition, y="values"))+
    geom_point()+
    stat_summary(fun="mean", geom="crossbar",
                 mapping=aes(ymin=..y.., ymax=..y..), width=0.22,
                 position=position_dodge(),show.legend = FALSE, color="red")+
    facet_wrap(~nhoods.id, ncol = n_col)+
    ylab("# cells in neighbourhood")

  return(p)
}
