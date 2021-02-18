
#' Identify post-hoc neighbourhood marker genes
#'
#' This function will perform differential gene expression analysis on
#' groups of neighbourhoods. Adjacent and concordantly DA neighbourhoods can be defined using
#' \code{groupNhoods} or by the user. Cells \emph{between} these
#' aggregated groups are compared. For differential gene experession based on an input design
#' \emph{within} DA neighbourhoods see \code{\link{testDiffExp}}.
#'
#' @param x A \code{\linkS4class{Milo}} object containing single-cell gene expression
#' and neighbourhoods.
#' @param da.res A \code{data.frame} containing DA results, as expected from running
#' \code{testNhoods}, as a \code{NhoodGroup} column specifying the grouping of neighbourhoods,
#' as expected from
#' @param da.fdr A numeric scalar that determines at what FDR neighbourhoods are declared
#' DA for the purposes of aggregating across concorantly DA neighbourhoods.
#' @param assay A character scalar determining which \code{assays} slot to extract from the
#' \code{\linkS4class{Milo}} object to use for DGE testing.
#' @param aggregate.samples logical indicating wheather the expression values for cells in the same sample
#' and neighbourhood group should be merged for DGE testing. This allows to perform testing exploiting the replication structure
#' in the experimental design, rather than treating single-cells as independent replicates. The function used for aggregation depends on the
#' selected gene expression assay: if \code{assay="counts"} the expression values are summed, otherwise we take the mean.
#' @param sample_col a character scalar indicating the column in the colData storing sample information
#' (only relevant if \code{aggregate.samples==TRUE})
#' @param overlap A scalar integer that determines the number of cells that must
#' overlap between adjacent neighbourhoods for merging.
#' @param lfc.threshold A scalar that determines the absolute log fold change above
#' which neighbourhoods should be considerd 'DA' for merging. Default=NULL
#' @param merge.discord A logical scalar that overrides the default behaviour and allows
#' adjacent neighbourhoods to be merged if they have discordant log fold change signs. Using
#' this argument is generally discouraged, but may be useful for constructing an empirical null
#' group of cells, regardless of DA sign.
#' @param subset.row A logical, integer or character vector indicating the rows
#' of \code{x} to use for sumamrizing over cells in neighbourhoods.
#' @param gene.offset A logical scalar the determines whether a per-cell offset
#' is provided in the DGE GLM to adjust for the number of detected genes with
#' expression > 0.
#' @param return.groups A logical scalar that returns a \code{\link{data.frame}} of the
#' aggregated groups per single-cell. Cells that are members of non-DA neighbourhoods contain
#' \code{NA} values.
#' @param subset.nhoods A logical, integer or character vector indicating which neighbourhoods
#' to subset before aggregation and DGE testing (default: NULL).
#' @param subset.groups A character vector indicating which groups to test for markers (default: NULL)
#' @param na.function A valid NA action function to apply, should be one of
#' \code{na.fail, na.omit, na.exclude, na.pass}.
#' @param compute.new A logical scalar indicating whether to force computing a new neighbourhood
#' adjacency matrix if already present.
#' 
#' 
#' @details
#' Using a one vs. all approach, each aggregated group of cells is compared to all others
#' using the single-cell log normalized gene expression with a GLM
#' (for details see \code{\link[limma]{limma-package}}), or the single-cell counts using a
#' negative binomial GLM (for details see \code{\link[edgeR]{edgeR-package}}). When using
#' the latter it is recommended to set \code{gene.offset=TRUE} as this behaviour adjusts
#' the model offsets by the number of detected genes in each cell.
#'
#'
#' @return A \code{data.frame} of DGE results containing a log fold change and adjusted
#' p-value for each aggregated group of neighbourhoods. If \code{return.groups} then
#' the return value is a list with the slots \code{groups} and \code{dge} containing the
#' aggregated neighbourhood groups per single-cell and marker gene results, respectively.
#'
#' \emph{Warning}: If all neighbourhoods are grouped together, then it is impossible to
#' run \code{findNhoodMarkers}. In this (hopefully rare) instance, this function will return
#' a warning and return \code{NULL}.
#'
#' @author Emma Dann
#'
#' @examples
#'
#' @name findNhoodGroupMarkers
#' @export
#' @importFrom stats model.matrix as.formula
#' @importFrom Matrix colSums
findNhoodGroupMarkers <- function(x, da.res, assay="logcounts",
                             aggregate.samples=FALSE, sample_col=NULL,
                             subset.row=NULL, gene.offset=TRUE,
                             subset.nhoods=NULL, subset.groups=NULL,
                             na.function="na.pass", compute.new=FALSE){
  
  if(!is(x, "Milo")){
    stop("Unrecognised input type - must be of class Milo")
  } else if(any(!assay %in% assayNames(x))){
    stop(paste0("Unrecognised assay slot: ", assay))
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
  
  if (isTRUE(aggregate.samples) & is.null(sample_col)) {
    stop("if aggregate.samples is TRUE, the column storing sample information must be specified by setting 'sample_col'")
  }
  
  if (!"NhoodGroup" %in% colnames(da.res)) {
    stop("'NhoodGroup' columns is missing from da.res. Please run groupNhoods() or define neighbourhood groupings otherwise.")
  }
  
  nhs.da.gr <- da.res$NhoodGroup
  names(nhs.da.gr) <- da.res$Nhood
  
  if(!is.null(subset.nhoods)){
    nhs.da.gr <- nhs.da.gr[subset.nhoods]
  }
  
  nhood.gr <- unique(nhs.da.gr)
  # perform DGE _within_ each group of cells using the input design matrix
  fake.meta <- data.frame("CellID"=colnames(x), "Nhood.Group"=rep(NA, ncol(x)))
  rownames(fake.meta) <- fake.meta$CellID
  
  # do we want to allow cells to be members of multiple groups? This will create
  # chaos for the LM as there will be a dependency structure comparing 2 different
  # groups that contain overlapping cells.
  # this approach means that the latter group takes precedent.
  # maybe exclude the cells that fall into separate groups?
  
  for(i in seq_along(nhood.gr)){
    nhood.x <- da.res$Nhood[which(nhs.da.gr == nhood.gr[i])]
    
    # get the nhoods
    nhs <- nhoods(x)
    if(!is.null(subset.nhoods)){
      nhs <- nhs[,subset.nhoods]
    }
    
    nhood.gr.cells <- rowSums(nhs[, nhood.x, drop=FALSE]) > 0
    ## set group to NA if a cell was already assigned to a group
    fake.meta[nhood.gr.cells,"Nhood.Group"] <- ifelse(is.na(fake.meta[nhood.gr.cells,"Nhood.Group"]), nhood.gr[i], NA)
  }
  
  # # only compare against the other DA neighbourhoods
  # x <- x[, !is.na(fake.meta$Nhood.Group)]
  # fake.meta <- fake.meta[!is.na(fake.meta$Nhood.Group), ]
  
  if(!is.null(subset.row)){
    x <- x[subset.row, , drop=FALSE]
  }
  
  exprs <- assay(x, assay)
  
  marker.list <- list()
  i.contrast <- c("TestTest - TestRef") # always use contrasts for this
  
  # if there is only 1 group, then need to make sure that all neighbourhoods
  # are not in this group - otherwise can't do any DGE testing
  if(length(nhood.gr) == 1){
    if(sum(fake.meta$Nhood.Group == nhood.gr[1]) == nrow(fake.meta)){
      warning("All graph neighbourhoods are in the same group - cannot perform DGE testing. Returning NULL")
      return(NULL)
    }
  }

  ## Aggregate expression by sample
  # To avoid treating cells as independent replicates
  if (isTRUE(aggregate.samples)) {
    fake.meta[,"sample_id"] <- colData(x)[[sample_col]]
    fake.meta[,'sample_group'] <- paste(fake.meta[,"sample_id"], fake.meta[,"Nhood.Group"], sep="_")
    
    sample_gr_mat <- matrix(0, nrow=nrow(fake.meta), ncol=length(unique(fake.meta$sample_group)))
    colnames(sample_gr_mat) <- unique(fake.meta$sample_group)
    rownames(sample_gr_mat) <- rownames(fake.meta)
    
    for (s in colnames(sample_gr_mat)) {
      sample_gr_mat[which(fake.meta$sample_group == s),s] <- 1
    }
    
    ## Summarise expression by sample
    exprs_smp <- matrix(0, nrow=nrow(exprs), ncol=ncol(sample_gr_mat))
    if (assay=='counts') {
      summFunc <- rowSums
    } else {
      summFunc <- rowMeans
    }
    
    for (i in 1:ncol(sample_gr_mat)){
      if (sum(sample_gr_mat[,i]) > 1) {
        exprs_smp[,i] <- summFunc(exprs[,which(sample_gr_mat[,i] > 0)])
      } else {
        exprs_smp[,i] <- exprs[,which(sample_gr_mat[,i] > 0)]
      }
    }
    rownames(exprs_smp) <- rownames(exprs)
    colnames(exprs_smp) <- colnames(sample_gr_mat)
    
    smp_meta <- unique(fake.meta[,c("sample_group","Nhood.Group")])
    rownames(smp_meta) <- smp_meta[,"sample_group"]
    
    fake.meta <- smp_meta
    exprs <- exprs_smp
  }
  
  ## Test for nhood groups markers
  ## Subset to tests of interest
  if(!is.null(subset.groups)){
    nhood.gr <- subset.groups  
  }
  
  for(i in seq_along(nhood.gr)){
    i.meta <- fake.meta
    i.meta$Test <- "Ref"
    i.meta$Test[fake.meta$Nhood.Group == nhood.gr[i]] <- "Test"
    
    if(ncol(exprs) > 1 & nrow(i.meta) > 1){
      i.design <- as.formula(" ~ 0 + Test")
      i.model <- model.matrix(i.design, data=i.meta)
      rownames(i.model) <- rownames(i.meta)
    }
    
    sink(file="/dev/null")
    gc()
    sink(file=NULL)
    
    if(assay == "logcounts"){
      i.res <- .perform_lognormal_dge(exprs, i.model, model.contrasts=i.contrast,
                                      gene.offset=gene.offset)
    } else if(assay == "counts"){
      i.res <- .perform_counts_dge(exprs, i.model, model.contrasts=i.contrast,
                                   gene.offset=gene.offset)
      colnames(i.res)[ncol(i.res)] <- "adj.P.Val"
    } else{
      warning("Assay type is not counts or logcounts - assuming (log)-normal distribution. Use these results at your peril")
      i.res <- .perform_lognormal_dge(exprs, i.model,
                                      model.contrasts=i.contrast,
                                      gene.offset=gene.offset)
    }
    
    i.res$adj.P.Val[is.na(i.res$adj.P.Val)] <- 1
    i.res$logFC[is.infinite(i.res$logFC)] <- 0
    
    i.res <- i.res[, c("logFC", "adj.P.Val")]
    colnames(i.res) <- paste(colnames(i.res), nhood.gr[i], sep="_")
    marker.list[[paste0(nhood.gr[i])]] <- i.res
    
    sink(file="/dev/null")
    gc()
    sink(file=NULL)
  }
  
  marker.df <- do.call(cbind.data.frame, marker.list)
  colnames(marker.df) <- gsub(colnames(marker.df), pattern="^[0-9]+\\.", replacement="")
  marker.df$GeneID <- rownames(i.res)
  
  return(marker.df)
  }
