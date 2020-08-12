#' Perform differential neighbourhood abundance testing
#'
#' This will perform differential neighbourhood abundance testing after cell
#' counting.
#' @param x A \code{\linkS4class{Milo}} object with a non-empty
#' \code{neighbourhoodCounts} slot.
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
#' vertex, neighbour-distance or k-distance (default).
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
#' @name testNeighbourhoods
NULL


#' @export
#' @importFrom stats model.matrix
#' @importFrom Matrix colSums
#' @importFrom stats dist median
#' @importFrom limma makeContrasts
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags
testNeighbourhoods <- function(x, design, design.df,
                               fdr.weighting=c("k-distance", "neighbour-distance", "edge", "vertex"),
                               min.mean=0, model.contrasts=NULL){

    if(class(design) == "formula"){
        model <- model.matrix(design, data=design.df)
        rownames(model) <- rownames(design.df)
    } else if(class(design) == "matrix"){
        model <- design
    }

    if(class(x) != "Milo"){
        stop("Unrecognised input type - must be of class Milo")
    } else if(.is_empty(x, "neighbourhoodCounts")){
        stop("Neighbourhood counts missing - please run countCells first")
    }

    # need to assume rownames of model are already set?
    # could check and warn at least if they aren't equal
    if(rownames(model) != rownames(design.df)){
        warning("Design matrix and design.df matrix dimnanes are not the same")
    }

    # assume neighbourhoodCounts and model are in the same order
    # cast as DGEList doesn't accept sparse matrices
    # what is the cost of cast a matrix that is already dense vs. testing it's class
    dge <- DGEList(counts=neighbourhoodCounts(x),
                   lib.size=log(colSums(neighbourhoodCounts(x))))

    dge <- estimateDisp(dge, model)
    fit <- glmQLFit(dge, model, robust=TRUE)
    if(!is.null(model.contrasts)){
        mod.constrast <- makeContrasts(model.contrasts, levels=model)
        res <- as.data.frame(topTags(glmQLFTest(fit, contrast=model.contrasts),
                                     sort.by='none', n=Inf))
    } else{
        n.coef <- ncol(model)
        res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))
    }

    res$Neighbourhood <- as.numeric(rownames(res))
    message("Performing spatial FDR correction")
    mod.spatialfdr <- graphSpatialFDR(nhoods=neighbourhoods(x),
                                      graph=graph(x),
                                      weighting=fdr.weighting,
                                      pvalues=res[order(res$Neighbourhood), ]$PValue,
                                      indices=neighbourhoodIndex(x),
                                      distances=neighbourDistances(x))

    res$SpatialFDR[order(res$Neighbourhood)] <- mod.spatialfdr
    res
}


#' @importFrom methods slot
#' @importFrom Matrix rowSums
.is_empty <- function(x, attribute){
    # check if a Milo object slot is empty or not
    x.slot <- slot(x, attribute)

    if(class(x.slot) == "list" & names(slot(sim1.mylo, "graph")) == "graph"){
        return(length(x.slot[[1]]) > 0)
    } else if(class(x.slot) == "list" & is.null(names(x.slot))){
        return(length(x.slot))
    } else if(any(class(x.slot) %in% c("dgCMatrix", "dsCMatrix"))){
        return(sum(rowSums(x.slot)) == 0)
    }
}
