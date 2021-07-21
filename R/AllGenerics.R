
#' @export
setGeneric("graph", function(x) standardGeneric("graph"))

#' @export
setGeneric("graph<-", function(x, value) standardGeneric("graph<-"))


#' @export
setGeneric("nhoodDistances", function(x) standardGeneric("nhoodDistances"))

#' @export
setGeneric("nhoodDistances<-", function(x, value) standardGeneric("nhoodDistances<-"))


#' @export
setGeneric("nhoods", function(x) standardGeneric("nhoods"))

#' @export
setGeneric("nhoods<-", function(x, value) standardGeneric("nhoods<-"))


#' @export
setGeneric("nhoodCounts", function(x) standardGeneric("nhoodCounts"))

#' @export
setGeneric("nhoodCounts<-", function(x, value) standardGeneric("nhoodCounts<-"))


#' @export
setGeneric("nhoodIndex", function(x) standardGeneric("nhoodIndex"))

#' @export
setGeneric("nhoodIndex<-", function(x, value) standardGeneric("nhoodIndex<-"))


#' @export
setGeneric("nhoodExpression", function(x) standardGeneric("nhoodExpression"))

#' @export
setGeneric("nhoodExpression<-", function(x, value) standardGeneric("nhoodExpression<-"))

#' @export
setGeneric("nhoodReducedDim", function(x, value="PCA") standardGeneric("nhoodReducedDim"))

#' @export
setGeneric("nhoodReducedDim<-", function(x, rdim="PCA", ..., value) standardGeneric("nhoodReducedDim<-"))

#' @export
setGeneric("nhoodGraph", function(x) standardGeneric("nhoodGraph"))

#' @export
setGeneric("nhoodGraph<-", function(x, value) standardGeneric("nhoodGraph<-"))

#' @export
setGeneric("nhoodAdjacency", function(x) standardGeneric("nhoodAdjacency"))

#' @export
setGeneric("nhoodAdjacency<-", function(x, value) standardGeneric("nhoodAdjacency<-"))
