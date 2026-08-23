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
.parse_formula <- function(in.form, design.df, vtype=c("re", "fe"), add.int=FALSE){
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

        # add the residual variance term if appropriate
        if(isTRUE(add.int)){
            d.mat <- cbind(d.mat, matrix(data=1L, nrow=nrow(d.mat), ncol=1))
            colnames(d.mat) <- c(trimws(v.terms), "residual")
        } else{
            colnames(d.mat) <- trimws(v.terms)
        }

    } else if(vtype %in% c("fe")){
        v.terms <- trimws(unlist(sp.form[!grepl(trimws(sp.form), pattern="~|\\|")]))
        if(length(v.terms) > 1){
            v.terms <- paste(v.terms, collapse=" + ")
        }

        # the intercept is a fixed effect in this model
        d.mat <- model.matrix(as.formula(paste("~ 1 +", v.terms)), data = design.df)
        d.mat <- d.mat[ ,!grepl("1*\\|", colnames(d.mat))] # drop the random terms

        if(isFALSE(add.int)){
            # drop the intercept term if required.
            d.mat[, !grepl("Intercept", colnames(d.mat)), drop=FALSE]
        }

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

#' Check a mixed model design for degenerate variance components
#'
#' A variance component is only identifiable if its design contributes
#' something the fixed effects do not. If the columns of Z_j lie in the column
#' space of X then sigma_j carries no information: the fit either fails on a
#' singular Hessian or converges to an arbitrary value with the Satterthwaite
#' degrees of freedom collapsing to zero and every p-value equal to 1. Neither
#' is a useful answer, and the second is easy to miss, so the design is
#' rejected up front.
#'
#' rank([X Z_j]) is always at least rank(X); equality means Z_j adds nothing.
#'
#' @param x.model matrix of fixed effects
#' @param z.model matrix mapping observations to random effect levels, one
#' column per random effect variable
#' @param geno.only logical, whether the only random effect is the genetic one,
#' whose design is the identity and so is never contained in the span of X
#'
#' @return invisible NULL, called for the side effect of erroring
#'
#' @importFrom stats model.matrix
.checkDesignRank <- function(x.model, z.model, geno.only=FALSE){
    x.rank <- qr(x.model)$rank

    if(x.rank < ncol(x.model)){
        stop("Fixed effect design matrix is rank deficient: rank ", x.rank,
             " with ", ncol(x.model), " columns. Two or more fixed effects are ",
             "collinear - drop one, or combine them into a single variable.")
    }

    # the genetic design is the identity, which is never contained in span(X)
    if(isTRUE(geno.only)){
        return(invisible(NULL))
    }

    for(j in seq_len(ncol(z.model))){
        zj <- model.matrix(~ 0 + factor(z.model[, j]))
        if(qr(cbind(x.model, zj))$rank <= x.rank){
            stop("Random effect '", colnames(z.model)[j], "' is collinear with the ",
                 "fixed effects: its levels are fully determined by the fixed effect ",
                 "design, so its variance component is not identifiable. Remove it ",
                 "from the random effects, or remove the corresponding fixed effect.")
        }
    }

    invisible(NULL)
}
