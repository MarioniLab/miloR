#' Count cells in neighbourhoods
#'
#' This function quantifies the number of cells in each neighbourhood according
#' to an input experimental design. This forms the basis for the differential
#' neighbourhood abundance testing.
#' @param x A \code{\linkS4class{Milo}} object with non-empty \code{graph} and
#' \code{nhoods} slots.
#' @param meta.data A cell X variable \code{data.frame}  containing study meta-data
#' including experimental sample IDs. Assumed to be in the same order as the
#' cells in the input \code{\linkS4class{Milo}} object.
#' @param samples Either a string specifying which column of \code{data}
#' should be used to identify the experimental samples for counting, or a
#' named vector of sample ids mapping each single cell to it's respective sample.
#' @details
#' This function generates a counts matrix of \code{nhoods} X samples,
#' and populates the \code{nhoodCounts} slot of the input
#' \code{\linkS4class{Milo}} object. This matrix is used down-stream for
#' differential abundance testing.
#'
#' @return A \code{\linkS4class{Milo}} object containing a counts matrix in the
#' \code{nhoodCounts} slot.
#'
#' @author Mike Morgan, Emma Dann
#'
#' @examples
#'
#' library(igraph)
#' m <- matrix(rnorm(100000), ncol=100)
#' milo <- buildGraph(t(m), k=20, d=10)
#' milo <- makeNhoods(milo, k=20, d=10, prop=0.3)
#'
#' cond <- rep("A", nrow(m))
#' cond.a <- sample(1:nrow(m), size=floor(nrow(m)*0.25))
#' cond.b <- setdiff(1:nrow(m), cond.a)
#' cond[cond.b] <- "B"
#' meta.df <- data.frame(Condition=cond, Replicate=c(rep("R1", 330), rep("R2", 330), rep("R3", 340)))
#' meta.df$SampID <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
#' milo <- countCells(milo, meta.data=meta.df, samples="SampID")
#' milo
#'
#' @name countCells
NULL

#' @export
#' @rdname countCells
#' @importFrom Matrix Matrix
#' @importClassesFrom S4Vectors DataFrame
countCells <- function(x, samples, meta.data=NULL){

    # cast dplyr objects to data.frame
    if(!is.data.frame(meta.data) & !is.null(meta.data)){
        meta.data <- as.data.frame(meta.data)
    }

    if(length(samples) > 1 & !is.null(meta.data)){
        stop("Multiple sample columns provided, please specify a unique column name")
    } else if(is.null(meta.data) & length(samples) != ncol(x)){
        stop(paste0("Length of vector does not match dimensions of object. Length:",
                    length(samples), " Dimensions: ", ncol(x)))
    }

    # check the nhoods slot is populated
    if(ncol(nhoods(x)) == 1 & nrow(nhoods(x)) == 1){
        stop("No neighbourhoods found. Please run makeNhoods() first.")
    }

    message("Checking meta.data validity")
    if(!is.null(meta.data)){
        if (is.factor(meta.data[, samples])){
            samp.ids <- levels(meta.data[, samples])
        } else {
            samp.ids <- unique(as.character(meta.data[, samples]))
        }
    } else {
        if (is.factor(samples)){
            samp.ids <- levels(samples)
        } else {
            samp.ids <- unique(as.character(samples))
        }
    }

    num.hoods <- ncol(nhoods(x))

    ## Convert meta data to binary dummies in sparse matrix
    dummy.meta.data <- Matrix(data=0, nrow=nrow(meta.data), ncol = length(samp.ids), sparse = TRUE)
    colnames(dummy.meta.data) <- samp.ids
    rownames(dummy.meta.data) <- rownames(meta.data)
    for (s in seq_along(samp.ids)){
        i.s <- samp.ids[s]
        s.ixs <- which(meta.data[samples]==i.s)
        dummy.meta.data[s.ixs, as.character(i.s)] <- 1
    }

    message("Counting cells in neighbourhoods")
    count.matrix <- Matrix::t(nhoods(x)) %*% dummy.meta.data

    # add to the object
    rownames(count.matrix) <- c(1:num.hoods)
    nhoodCounts(x) <- count.matrix

    return(x)
}
