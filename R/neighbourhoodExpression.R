#' Average expression within neighbourhoods
#'
#' This function calculates the mean expression of each feature in the
#' Milo object stored in the assays slot. Neighbourhood expression data are
#' stored in a new slot \code{neighbourhoodExpression}.
#'
#' @param x A \code{Milo} object with neighbourhoods slot populated, alternatively a list
#' containing a vector of cell indices; one per neighbourhood.
#' @param subset.row A logical, integer or character vector indicating the rows
#' of \code{x} to use for sumamrizing over cells in neighbourhoods.
#' @param exprs If \code{x} is a list of neighbourhoods, \code{exprs} is a matrix of genes X
#' cells to use for calculating neighbourhood expression.
#'
#' @details
#' This function computes the mean expression of each gene, subset by \code{subset.rows}
#' where present, across the cells contained within each neighbourhood.
#'
#' @return A \code{\linkS4class{Milo}} object with the neighbourhoodExpression slot populated.
#'
#' @author
#' Mike Morgan
#'
#' @examples
#' require(SingleCellExperiment)
#' m <- matrix(rnorm(100000), ncol=100)
#' milo <- Milo(SingleCellExperiment(assays(logcounts=m)))
#' milo <- buildGraph(m, d=30, transposed=TRUE)
#' milo <- makeNeighbourhoods(milo)
#' milo <- calcNeighbourhoodExpression(milo)
#' dim(neighbourhoodExpression(milo))
#'
#' @name calcNeighbourhoodExpression
NULL

#' @export
#' @rdname buildGraph
calcNeighbourhoodExpression <- function(x, assay="logcounts", subset.row=NULL, exprs=NULL){

    if(class(x) == "Milo"){
        # are neighbourhoods calculated?
        if(length(neighbourhoods(x)) == 0){
            stop("No neighbourhoods found - run makeNeighbourhoods first")
        }

        if(!is.null(assay(x, assay))){
            n.exprs <- .calc_expression(neighbourhoods=neighbourhoods(x),
                                        data.set=assay(x, assay),
                                        subset.row=subset.row)
            neighbourhoodExpression(x) <- n.exprs
            return(x)
        }
    } else if(class(x) == "list"){
        if(is.null(exprs)){
            stop("No expression data found. Please specific a gene expression matrix to exprs")
        } else{
            n.exprs <- .calc_expression(neighbourhoods=x,
                                        data.set=exprs,
                                        subset.row=subset.row)
            x.milo <- Milo(SingleCellExperiment(assay=list(logcounts=exprs)))
            neighbourhoodExpression(x.milo) <- n.exprs
            return(x.milo)
        }
    }
}


#' @importFrom Matrix colSums
.calc_expression <- function(neighbourhoods, data.set, subset.row=NULL){
    neighbour.model <- matrix(0L, ncol=length(neighbourhoods), nrow=ncol(data.set))

    for(x in seq_along(1:length(neighbourhoods))){
        neighbour.model[neighbourhoods[[x]], x] <- 1
    }

    if(!is.null(subset.row)){
        neigh.exprs <- t(neighbour.model) %*% t(data.set[subset.row, ])
    } else{
        neigh.exprs <- t(neighbour.model) %*% t(data.set)
    }
    neigh.exprs <- t(apply(neigh.exprs, 2, FUN=function(XP) XP/colSums(neighbour.model)))

    if(is.null(subset.row)){
        rownames(neigh.exprs) <- rownames(data.set)
    } else{
        rownames(neigh.exprs) <- rownames(data.set)[subset.row]
    }

    return(neigh.exprs)
}


