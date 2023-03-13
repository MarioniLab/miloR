### Annotate neighbourhood results ###

#' Add annotations from colData to DA testing results
#'
#' This function assigns a categorical label to neighbourhoods in the differential abundance results
#' data.frame (output of \code{testNhoods}), based on the most frequent label among cells in each
#' neighbourhood. This can be useful to stratify DA testing results by cell types or samples.
#' Also the fraction of cells carrying that label is stored.
#'
#' @param x A \code{\linkS4class{Milo}} object containing single-cell gene expression
#' and neighbourhoods.
#' @param da.res A \code{data.frame} containing DA results, as expected from running
#' \code{testNhoods}.
#' @param coldata_col A character scalar determining which column of \code{colData(x)} stores
#' the annotation to be added to the neighbourhoods
#' @param subset.nhoods A character, numeric or logical vector that will subset the annotation to the specific nhoods. If
#' a character vector these should correspond to row names of \code{nhoodCounts}. If a logical vector then
#' these should have the same \code{length} as \code{nrow} of \code{nhoodCounts}. If numeric, then these are assumed
#' to correspond to indices of \code{nhoodCounts} - if the maximal index is greater than \code{nrow(nhoodCounts(x))}
#' an error will be produced. This is necessary if \code{testNhoods} was run using \code{subset.nhoods=...}.
#'
#' @details
#' For each neighbourhood, this calculates the most frequent value of \code{colData(x)[coldata_col]}
#' among cells in the neighbourhood and assigns that value as annotation for the neighbourhood, adding a column in the
#' \code{da.res} data.frame. In addition, a \code{coldata_col_fraction} column will be added, storing the fraction of cells
#' carrying the assigned label. While in practice neighbourhoods are often homogeneous, one might choose to remove an
#' annotation label when the fraction of cells with the label is too low (e.g. below 0.6).
#'
#' @return A \code{data.frame} of model results (as \code{da.res} input) with two new columns: (1) \code{coldata_col} storing
#' the assigned label for each neighbourhood; (2) \code{coldata_col_fraction} storing the fraction of cells in the neighbourhood with
#' the assigned label.
#'
#' @author
#' Emma Dann
#'
#' @examples
#'
#' NULL
#'
#' @export
#' @rdname annotateNhoods
annotateNhoods <- function(x, da.res, coldata_col, subset.nhoods=NULL){
    if(!is(x, "Milo")){
        stop("Unrecognised input type - must be of class Milo")
    }

    if(!coldata_col %in% names(colData(x))){
        stop(coldata_col, " is not a column in colData(x)")
    }

    if(is.null(subset.nhoods)){
        if(ncol(nhoods(x)) != nrow(da.res)){
            stop("the number of rows in da.res does not match the number of neighbourhoods in nhoods(x). Are you sure da.res is the output of testNhoods(x) or did you use subset.nhoods?")
        }
        keep.nh <- rep(TRUE, ncols(nhoods(x)))
    } else{
        keep.nh <- subset.nhoods
    }

  anno_vec <- colData(x)[[coldata_col]]
  if (!is.factor(anno_vec)) {
    message("Converting ", coldata_col, " to factor...")
    anno_vec <- factor(anno_vec, levels=unique(anno_vec))
  }

  ## Count occurrence of labels in each nhood
  n.levels <- length(levels(anno_vec))
  nhood_counts <- vapply(seq_len(ncol(nhoods(x)[, keep.nh, drop=FALSE])), FUN=function(n) table(anno_vec[which(nhoods(x)[ ,n]==1)]),
                         FUN.VALUE=numeric(n.levels))
  nhood_counts <- t(nhood_counts)
  rownames(nhood_counts) <- seq_len(ncol(nhoods(x)[, keep.nh, drop=FALSE]))

  ## Fetch the most frequent label
  max_val <- apply(nhood_counts, 1, function(X) colnames(nhood_counts)[which.max(X)])
  max_frac <- apply(nhood_counts, 1, function(X) max(X)/sum(X))

  ## Add to da.res table
  da.res[coldata_col] <- max_val
  da.res[paste0(coldata_col, "_fraction")] <- max_frac

  return(da.res)
  }



