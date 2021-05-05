### Grouping neighbourhoods ###

#' Group neighbourhoods
#'
#' This function groups overlapping and concordantly DA neighbourhoods, using the louvain
#' community detection algorithm.
#'
#' @param x A \code{\linkS4class{Milo}} object containing single-cell gene expression
#' and neighbourhoods.
#' @param da.res A \code{data.frame} containing DA results, as expected from running
#' \code{testNhoods}.
#' @param da.fdr A numeric scalar that determines at what FDR neighbourhoods are declared
#' DA for the purposes of aggregating across concorantly DA neighbourhoods.
#' @param overlap A scalar integer that determines the number of cells that must
#' overlap between adjacent neighbourhoods for merging.
#' @param max.lfc.delta A scalar that determines the absolute difference in log fold change below
#' which neighbourhoods should not be considered adjacent. Default=NULL
#' @param merge.discord A logical scalar that overrides the default behaviour and allows
#' adjacent neighbourhoods to be merged if they have discordant log fold change signs. Using
#' this argument is generally discouraged, but may be useful for constructing an empirical null
#' group of cells, regardless of DA sign.
#' @param subset.nhoods A logical, integer or character vector indicating which neighbourhoods
#' to subset before grouping. All other neighbourhoods will be assigned NA
#' @param compute.new A logical scalar indicating whether to force computing a new neighbourhood
#' adjacency matrix if already present.
#' @param na.function A valid NA action function to apply, should be one of
#' \code{na.fail, na.omit, na.exclude, na.pass} (default='na.pass').
#'
#' @return A \code{data.frame} of model results (as \code{da.res} input) with a new column storing the assigned
#' group label for each neighbourhood (\code{NhoodGroup} column)
#'
#' @details Louvain clustering is applied to the neighbourhood graph. This graph is first modified
#' based on two criteria: 1) neighbourhoods share at least \code{overlap} number of cells,
#' and 2) the DA log fold change sign is concordant.
#' This behaviour can be modulated by setting \code{overlap} to be more or less stringent.
#' Additionally, a threshold on the log fold-changes can be set, such that \code{max.lfc.delta}
#' is required to retain edges between adjacent neighbourhoods. Note: adjacent neighbourhoods will
#' never be merged with opposite signs.

#'
#' @author Emma Dann & Mike Morgan
#' @export
groupNhoods <- function(x, da.res, da.fdr=0.1,
                        overlap=1, max.lfc.delta=NULL,
                        merge.discord=FALSE,
                        subset.nhoods=NULL,
                        compute.new=FALSE,
                        na.function="na.pass"
                        ){
  if(!is(x, "Milo")){
    stop("Unrecognised input type - must be of class Milo")
  }

  if(is.null(na.function)){
    warning("NULL passed to na.function, using na.pass")
    na.func <- get("na.pass")
  } else{
    tryCatch({
      na.func <- get(na.function)
    }, warning=function(warn){
      warning(warn)
    }, error=function(err){
      stop(paste0("NA function ", na.function, " not recognised"))
    }, finally={
    })
  }

  n.da <- sum(na.func(da.res$SpatialFDR < da.fdr))

  if(!is.na(n.da) & n.da == 0){
    stop("No DA neighbourhoods found")
  }

  if(any(is.na(da.res$SpatialFDR))){
    warning("NA values found in SpatialFDR vector")
  }

  message("Found ", n.da, " DA neighbourhoods at FDR ", da.fdr*100, "%")

  ## Check if adjacency matrix exists, if not build
  if((ncol(nhoodAdjacency(x)) == ncol(nhoods(x))) & isFALSE(compute.new)){
    message("nhoodAdjacency found - using for nhood grouping")
  } else {
    x <- buildNhoodGraph(x, overlap = overlap)
  }

  ## Make neighbourhood groups
  nhs_groups <- .group_nhoods_from_adjacency(nhoods(x),
                                             nhood.adj=nhoodAdjacency(x),
                                             da.res=da.res,
                                             is.da=da.res$SpatialFDR < da.fdr,
                                             merge.discord=merge.discord,
                                             max.lfc.delta=max.lfc.delta,
                                             overlap=overlap,
                                             subset.nhoods=subset.nhoods
                                             )

  ## Save in DAres data.frame
  da.res['NhoodGroup'] <- NA
  if (!is.null(subset.nhoods)) {
    da.res[subset.nhoods,"NhoodGroup"] <- as.character(nhs_groups)
  } else {
    da.res['NhoodGroup'] <- as.character(nhs_groups)
  }
  return(da.res)
}

#' @importFrom igraph graph_from_adjacency_matrix components cluster_louvain
.group_nhoods_from_adjacency <- function(nhs, nhood.adj, da.res, is.da,
                                         merge.discord=FALSE,
                                         max.lfc.delta=NULL,
                                         overlap=1,
                                         subset.nhoods=NULL
                                         ){

  if(is.null(colnames(nhs))){
    warning("No names attributed to nhoods. Converting indices to names")
    colnames(nhs) <- as.character(c(1:ncol(nhs)))
  }

  # assume order of nhs is the same as nhood.adj
  if(!is.null(subset.nhoods)){
    if(mode(subset.nhoods) %in% c("character", "logical", "numeric")){
      # force use of logicals for consistency
      if(mode(subset.nhoods) %in% c("character")){
        sub.log <- colnames(nhs) %in% subset.nhoods
      } else if (mode(subset.nhoods) %in% c("numeric")) {
        sub.log <- colnames(nhs) %in% colnames(nhs)[subset.nhoods]
      } else{
        sub.log <- subset.nhoods
      }

      nhood.adj <- nhood.adj[sub.log, sub.log]

      if(length(is.da) == ncol(nhs)){
        nhs <- nhs[sub.log]
        is.da <- is.da[sub.log]
        da.res <- da.res[sub.log, ]
      } else{
        stop("Subsetting `is.da` vector length does not equal nhoods length")
      }
    } else{
      stop(paste0("Incorrect subsetting vector provided:", class(subset.nhoods)))
    }
  } else{
    if(length(is.da) != ncol(nhood.adj)){
      stop("Subsetting `is.da` vector length is not the same dimension as adjacency")
    }
  }

  ## check for concordant signs (only for significant DA) - assume order is the same as nhoods
  if(isFALSE(merge.discord)){
    discord.sign <- sign(da.res[is.da, 'logFC'] %*% t(da.res[is.da, 'logFC'])) < 0
    nhood.adj[is.da, is.da][discord.sign] <- 0
  }

  if(overlap > 1){
    nhood.adj[nhood.adj < overlap] <- 0
  }

  ## Remove edges if the difference is higher than max.lfc.delta
  if(!is.null(max.lfc.delta)){
    lfc.diff <- sapply(da.res[,"logFC"], "-", da.res[,"logFC"])
    nhood.adj[abs(lfc.diff) > max.lfc.delta] <- 0
  }

  # binarise
  nhood.adj <- as.matrix((nhood.adj > 0) + 0)

  n.dim <- ncol(nhood.adj)
  if(!isSymmetric(nhood.adj)){
    stop("Overlap matrix is not symmetric")
  }

  if(nrow(nhood.adj) != ncol(nhood.adj)){
    stop("Non-square distance matrix - check nhood subsetting")
  }

  g <- graph_from_adjacency_matrix(nhood.adj, mode="undirected", diag=FALSE)
  groups <- cluster_louvain(g)$membership
  names(groups) <- colnames(nhood.adj)

  ## The user might actually want to compare DA vs Not DA
  # # only keep the groups that contain >= 1 DA neighbourhoods
  # keep.groups <- intersect(unique(groups[is.da]), unique(groups))
  # return(groups[groups %in% keep.groups])
  return(groups)
}

## Do we still need this??
#' @importFrom igraph graph_from_adjacency_matrix components
.group_nhoods_by_overlap <- function(nhs, da.res, is.da, overlap=1,
                                     max.lfc.delta=NULL, merge.discord=FALSE,
                                     subset.nhoods=NULL, cells=NULL){

  ## Build adjacency matrix for nhoods
  nhood.adj <- .build_nhood_adjacency(nhs)
  groups <- .group_nhoods_from_adjacency(nhs=nhs, nhood.adj=nhood.adj,
                                         is.da=is.da, da.res=da.res,
                                         subset.nhoods=subset.nhoods,
                                         overlap=overlap,
                                         max.lfc.delta=max.lfc.delta,
                                         merge.discord=merge.discord)

  return(groups)
}
