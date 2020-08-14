#' Perform differential neighbourhood abundance testing
#'
#' This will perform differential neighbourhood abundance testing after cell
#' counting.
#' @param x A \code{\linkS4class{Milo}} object with a non-empty
#' \code{nhoodCounts} slot.
#' @param design A \code{formula} or \code{model.matrix} object describing the
#' experimental design for differential abundance testing. The last component
#' of the formula or last column of the model matrix are by default the test
#' variable. This behaviour can be overridden by setting the \code{model.contrasts}
#' argument
#' @param design.df A \code{data.frame} containing meta-data to which \code{design}
#' refers to
#' @param min.mean A scalar used to threshold neighbourhoods on the minimum
#' average cell counts across samples.
#' @param model.contrasts A string vector that defines the contrasts used to perform
#' DA testing.
#' @param fdr.weighting The spatial FDR weighting scheme to use. Choice from edge,
#' vertex, neighbour-distance or k-distance (default). If \code{none} is passed no
#' spatial FDR correction is performed and returns a vector of NAs.
#' @param seed Seed number used for pseudorandom number generators.
#'
#'
#' @details
#' This function wraps up several steps of differential abundance testing using
#' the \code{edgeR} functions. These could be performed separately for users
#' who want to exercise more contol over their DA testing. By default this
#' function sets the \code{lib.sizes} to the log10(colSums(x)), and uses the
#' Quasi-Likelihood F-test in \code{glmQLFTest} for DA testing. FDR correction
#' is performed separately as the default multiple-testing correction is
#' inappropriate for neighbourhoods with overlapping cells.
#'
#' @return A \code{\linkS4class{Milo}}
#'
#' @author Mike Morgan
#'
#' @examples
#' NULL
#'
#' @name testNhoods
NULL


#' @export
#' @importFrom stats model.matrix
#' @importFrom Matrix colSums rowMeans
#' @importFrom stats dist median
#' @importFrom limma makeContrasts
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags
testNhoods <- function(x, design, design.df,
                               fdr.weighting=c("k-distance", "neighbour-distance", "edge", "vertex", "none"),
                               min.mean=0, model.contrasts=NULL, seed=42){
    set.seed(seed)
    if(class(design) == "formula"){
        model <- model.matrix(design, data=design.df)
        rownames(model) <- rownames(design.df)
    } else if(class(design) == "matrix"){
        model <- design
        if(any(rownames(model) != rownames(design.df))){
            warning("Design matrix and design matrix dimnames are not the same")
        }
    }

    if(class(x) != "Milo"){
        stop("Unrecognised input type - must be of class Milo")
    } else if(.is_empty(x, "nhoodCounts")){
        stop("Neighbourhood counts missing - please run countCells first")
    }

    if(ncol(nhoodCounts(x)) != nrow(model)){
        stop(paste0("Design matrix (", nrow(model), ") and nhood counts (",
                    ncol(nhoodCounts(x)), ") are not the same dimension"))
    }

    # assume nhoodCounts and model are in the same order
    # cast as DGEList doesn't accept sparse matrices
    # what is the cost of cast a matrix that is already dense vs. testing it's class
    if(min.mean > 0){
        keep.nh <- rowMeans(nhoodCounts(x)) >= min.mean
    } else{
        keep.nh <- rep(TRUE, nrow(nhoodCounts(x)))
    }

    dge <- DGEList(counts=nhoodCounts(x)[keep.nh, ],
                   lib.size=log(colSums(nhoodCounts(x))))

    dge <- estimateDisp(dge, model)
    fit <- glmQLFit(dge, model, robust=TRUE)
    if(!is.null(model.contrasts)){
        mod.constrast <- makeContrasts(contrasts=model.contrasts, levels=model)
        res <- as.data.frame(topTags(glmQLFTest(fit, contrast=mod.constrast),
                                     sort.by='none', n=Inf))
    } else{
        n.coef <- ncol(model)
        res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))
    }

    res$Nhood <- as.numeric(rownames(res))
    message(paste0("Performing spatial FDR correction with", fdr.weighting, " weighting"))
    mod.spatialfdr <- graphSpatialFDR(x.nhoods=nhoods(x),
                                      graph=graph(x),
                                      weighting=fdr.weighting,
                                      pvalues=res[order(res$Nhood), ]$PValue,
                                      indices=nhoodIndex(x),
                                      distances=nhoodDistances(x),
                                      reduced.dimensions=reducedDim(x, "PCA"))

    res$SpatialFDR[order(res$Nhood)] <- mod.spatialfdr
    res
}


#' @importFrom methods slot
#' @importFrom Matrix rowSums
.is_empty <- function(x, attribute){
    # check if a Milo object slot is empty or not
    x.slot <- slot(x, attribute)

    if(class(x.slot) == "list" & names(slot(x, "graph")) == "graph"){
        return(length(x.slot[[1]]) > 0)
    } else if(class(x.slot) == "list" & is.null(names(x.slot))){
        return(length(x.slot))
    } else if(any(class(x.slot) %in% c("dgCMatrix", "dsCMatrix"))){
        return(sum(rowSums(x.slot)) == 0)
    }
}