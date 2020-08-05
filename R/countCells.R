#' Count cells in neighbourhoods
#'
#' This function quantifies the number of cells in each neighbourhood according
#' to an input experimental design. This forms the basis for the differential
#' neighbourhood abundance testing.
#' @param x A \code{\linkS4class{Milo}} object with non-empty \code{graph} and
#' \code{neighbourhoods} slots.
#' @param data A cell X variable \code{data.frame}  containing study meta-data
#' including experimental sample IDs. Assumed to be in the same order as the
#' cells in the input \code{\linkS4class{Milo}} object.
#' @param samples Either a string specifying which column of \code{data}
#' should be used to identify the experimental samples for counting, or a
#' named vector of sample ids mapping each single cell to it's respective sample.
#' @details
#' This function generates a counts matrix of \code{neighbourhoods} X samples,
#' and populates the \code{neighbourhoodCounts} slot of the input
#' \code{\linkS4class{Milo}} object. This matrix is used down-stream for
#' differential abundance testing.
#'
#' @return A \code{\linkS4class{Milo}} object containing a counts matrix in the
#' \code{neighbourhoodCounts} slot.
#'
#' @author Mike Morgan
#'
#' @examples
#'
#' requires(igraph)
#' m <- matrix(rnorm(10000), ncol=10)
#' milo <- buildGraph(m, d=10)
#' milo <- makeNeighbourhoods(milo, prop=0.3)
#'
#' cond <- rep("A", nrow(m))
#' cond.a <- sample(1:nrow(m), size=floor(nrow(m)*0.25))
#' cond.b <- setdiff(1:nrow(m), cond.a)
#' cond[cond.b] <- "B"
#' meta.df <- data.frame(Condition=cond,
#'                       Replicate=c(rep("R1", floor(nrow(m)*0.33)),
#'                                   rep("R2", floor(nrow(m)*0.33)),
#'                                   rep("R3", nrow(m)-(2\*floor(nrow(m)\*0.33))))
#'
#' meta.df$SampID <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
#' milo <- countCells(milo, data=meta.df, sample.column="SampID")
#' milo
#'
#' @name countCells
NULL

#' @export
#' @importFrom Matrix Matrix
#' @importClassesFrom S4Vectors DataFrame
countCells <- function(x, samples, data=NULL){

    if(length(samples) > 1 & !is.null(data)){
        stop("Multiple sample columns provided, please specify a unique column name")
    } else if(is.null(data) & length(samples) != ncol(x)){
        stop(paste0("Length of vector does not match dimensions of object. Length:",
                    length(samples), " Dimensions: ", ncol(x)))
    }

    # check the neighbourhoods slot is populated
    if(length(neighbourhoods(x)) == 0){
        stop("No neighbourhoods found. Please run makeNeighbourhoods() first.")
    }

    message("Checking data validity")
    if(!is.null(data)){
        samp.ids <- unique(data[, samples])
    } else{
        samp.ids <- unique(samples)
    }

    n.hoods <- length(neighbourhoods(x))
    message(paste0("Setting up matrix with ", n.hoods, " neighbourhoods"))
    count.matrix <- Matrix(0L, ncol=length(samp.ids), nrow=n.hoods, sparse=TRUE)
    colnames(count.matrix) <- samp.ids

    message("Counting cells in neighbourhoods")
    for(i in seq_along(1:n.hoods)){
        v.i <- neighbourhoods(x)[[i]]
        for(j in seq_along(1:length(samp.ids))){
            j.s <- samp.ids[j]

            if(is.null(data)){
                # samples is a vector of N cells
                j.s.vertices <- intersect(v.i, names(samples[samples == j.s]))
            } else{
                j.s.vertices <- intersect(v.i, which(data[, samples] == j.s))
            }
            count.matrix[i, j] <- length(j.s.vertices)
        }
    }

    # add to the object
    rownames(count.matrix) <- c(1:n.hoods)
    neighbourhoodCounts(x) <- count.matrix
    return(x)
}
