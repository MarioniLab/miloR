## These are utility functions not meant to be exposed to the user

#' @importFrom methods slot
#' @importFrom Matrix rowSums
.check_empty <- function(x, attribute){
    # check if a Milo object slot is empty or not
    x.slot <- slot(x, attribute)

    if(is.list(x.slot) & names(slot(x, "graph")) == "graph"){
        return(length(x.slot[[1]]) > 0)
    } else if(is.list(x.slot) & is.null(names(x.slot))){
        return(length(x.slot))
    } else if(any(class(x.slot) %in% c("dgCMatrix", "dsCMatrix", "ddiMatrix", "matrix"))){
        return(sum(rowSums(x.slot)) == 0)
    }
}


.check_binary <- function(x){
    # check if a matrix is binary or not
    sum.zeros <- sum(x == 0)
    sum.ones <- sum(x == 1)
    n.comps <- nrow(x) * ncol(x)

    return(sum(c(sum.zeros, sum.ones)) == n.comps)
}


#' @importFrom igraph make_graph simplify
.neighborsToKNNGraph <- function(nn, directed=FALSE) {
    start <- as.vector(row(nn))
    end <- as.vector(nn)
    interleaved <- as.vector(rbind(start, end))

    if (directed) {
        g <- make_graph(interleaved, directed=TRUE)

    } else {
        g <- make_graph(interleaved, directed=FALSE)
        g <- simplify(g, edge.attr.comb = "first")
    }
    g
}

# setting internals for replacement methods that require multiple arguments - borrowed from SingleCellExperiment
#' @importFrom methods slot
.set_reduced_dims <- function(x, value, slot.x=NULL, rdim=NULL){
    x <- updateObject(x)
    content <- slot(x, slot.x)

    if(slot.x == "nhoodReducedDim"){

        if(!is.null(rdim)){
            content[[rdim]] <- value
            x@nhoodReducedDim <- content
        } else{
            stop("No reduced dimensionality slot provided")
        }
    }else{
        stop("replacement method not implemented for ", slot)
    }

    x
}


# parse design formula
#' @export
.parse_formula <- function(in.form, design.df, vtype=c("re", "fe")){
    ## parse the formula and return the X and Z matrices
    # need to decide on how to handle intercept terms - i.e. FE or RE
    sp.form <- unlist(strsplit(as.character(in.form),
                               split="+", fixed=TRUE))

    if(vtype %in% c("re")){
        v.terms <- unlist(lapply(sp.form, FUN=function(sp) {
            return(ifelse(grepl(trimws(sp), pattern="\\|"), .rEParse(trimws(sp)), NA))
            }))
        v.terms <- v.terms[!is.na(v.terms)]
        d.mat <- as.matrix(design.df[, trimws(v.terms)])
        if (is.character(d.mat)) {
            d.mat <- matrix(unlist(lapply(data.frame(d.mat)[, , drop = FALSE],
                                          function(x) as.integer(factor(x)))), ncol = length(v.terms))
        }

        # add the residual variance term
        d.mat <- cbind(d.mat, matrix(data=1L, nrow=nrow(d.mat), ncol=1))
        colnames(d.mat) <- c(trimws(v.terms), "residual")

    } else if(vtype %in% c("fe")){
        v.terms <- trimws(unlist(sp.form[!grepl(trimws(sp.form), pattern="~|\\|")]))
        if(length(v.terms) > 1){
            v.terms <- paste(v.terms, collapse=" + ")
        }

        d.mat <- model.matrix(as.formula(paste("~ 1 +", v.terms)), data = design.df)
        d.mat <- d.mat[ ,!grepl("1*\\|", colnames(d.mat))] # drop the intercept term
    } else{
        stop("vtype ", vtype, " not recognised")
    }

    return(d.mat)
}


#' @export
.rEParse <- function(re.form) {

    .x <- gsub(unlist(strsplit(re.form, split="|", fixed=TRUE)),
               pattern="\\)", replacement="")

    return(.x[length(.x)])
}


######################################
## neighbourhood grouping functions
######################################

#### nhood adjacency matrix function
# Build adjacency matrix of overlap between neighbourhoods
#' @importFrom gtools permutations
#' @importFrom Matrix crossprod
.build_nhood_adjacency <- function(nhoods, overlap=1){
    nh_intersect_mat <- Matrix::crossprod(nhoods)
    nh_intersect_mat[nh_intersect_mat < overlap] <- 0

    rownames(nh_intersect_mat) <- colnames(nhoods)
    colnames(nh_intersect_mat) <- colnames(nhoods)
    return(nh_intersect_mat)
}

