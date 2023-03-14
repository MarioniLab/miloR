#' Check for separation of count distributions by variables
#'
#' Check the count distributions for each nhood according to a test
#' variable of interest. This is important for checking if there is separation
#' in the GLMM to inform either nhood subsetting or re-computation of the
#' NN-graph and refined nhoods.
#' @param x \code{\linkS4class{Milo}} object with a non-empty
#' \code{nhoodCounts} slot.
#' @param design.df A \code{data.frame} containing meta-data in which \code{condition}
#' is a column variable. The rownames must be the same as, or a subset of, the
#' colnames of \code{nhoodCounts(x)}.
#' @param condition A character scalar of the test variable contained in \code{design.df}.
#' This should be a factor variable if it is numeric or character it will be cast to a
#' factor variable.
#' @param min.val A numeric scalar that sets the minimum number of counts across condition level
#' samples, below which separation is defined.
#' @param factor.check A logical scalar that sets the factor variable level checking. See \emph{details}
#' for more information.
#'
#' @details
#' This function checks across nhoods for separation based on the separate levels
#' of an input factor variable. It checks if \emph{condition} is a factor variable,
#' and if not it will cast it to a factor. Note that the function first checks for the
#' number of unique values - if this exceeds > 50% of the number of elements an
#' error is generated. Users can override this behaviour with \code{factor.check=FALSE}.
#'
#' @return
#' A logical vector of the same length as \code{ncol(nhoodCounts(x))} where \emph{TRUE}
#' values represent nhoods where separation is detected. The output of this function
#' can be used to subset nhood-based analyses
#' e.g. \code{testNhoods(..., subset.nhoods=checkSepartion(x, ...))}.
#'
#' @author Mike Morgan
#'
#' @examples
#' library(SingleCellExperiment)
#' ux.1 <- matrix(rpois(12000, 5), ncol=400)
#' ux.2 <- matrix(rpois(12000, 4), ncol=400)
#' ux <- rbind(ux.1, ux.2)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#'
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'                             reducedDims=SimpleList(PCA=pca$x))
#'
#' milo <- Milo(sce)
#' milo <- buildGraph(milo, k=20, d=10, transposed=TRUE)
#' milo <- makeNhoods(milo, k=20, d=10, prop=0.3)
#' milo <- calcNhoodDistance(milo, d=10)
#'
#' cond <- rep("A", ncol(milo))
#' cond.a <- sample(1:ncol(milo), size=floor(ncol(milo)*0.25))
#' cond.b <- setdiff(1:ncol(milo), cond.a)
#' cond[cond.b] <- "B"
#' meta.df <- data.frame(Condition=cond, Replicate=c(rep("R1", 132), rep("R2", 132), rep("R3", 136)))
#' meta.df$SampID <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
#' milo <- countCells(milo, meta.data=meta.df, samples="SampID")
#'
#' test.meta <- data.frame("Condition"=c(rep("A", 3), rep("B", 3)), "Replicate"=rep(c("R1", "R2", "R3"), 2))
#' test.meta$Sample <- paste(test.meta$Condition, test.meta$Replicate, sep="_")
#' rownames(test.meta) <- test.meta$Sample
#'
#' check.sep <- checkSeparation(milo, design.df=test.meta, condition='Condition')
#' sum(check.sep)
#'
#' @name checkSeparation
NULL

#' @export
#'
checkSeparation <- function(x, design.df, condition, factor.check=TRUE){
    if(!any(colnames(design.df) %in% condition)){
        stop(condition, " is not a variable in design.df")
    }

    if(is.null(nhoodCounts(x))){
        stop("nhoodCounts not found - please run countCells() first")
    }

    if(is.null(rownames(design.df))){
        stop("Please add rownames to design.df that are the same as the colnames of nhoodCounts(x)")
    } else{
        if(sum(rownames(design.df) %in% colnames(nhoodCounts(x))) < 1){
            stop("rownames of design.df are not a subset of nhoodCounts colnames")
        }
    }

    if(isTRUE(factor.check)){
        if(!is.factor(design.df[, condition])){
            n.unique <- length(unique(as.vector(design.df[, condition])))
            if(n.unique >= floor(nrow(design.df)/2.0)){
                stop("Too many levels in ", condition, ". This is not a suitable variable")
            } else{
                cond.vec <- as.factor(design.df[, condition])
            }
        } else{
            cond.vec <- design.df[, condition]
        }
    } else{
        cond.vec <- as.function(design.df[, condition])
    }

    names(cond.vec) <- rownames(design.df)

    any_separate <- apply(nhoodCounts(x)[, names(cond.vec)],
                          FUN=function(NR, conditions, minval=1) {
                              cond.levels <- levels(conditions)
                              nr.tab <- unlist(by(NR, INDICES=conditions, FUN=sum, simplify=FALSE))
                              any(nr.tab < minval)
                          }, MARGIN=1, conditions=cond.vec, minval=min.val)
    return(any_separate)
}



