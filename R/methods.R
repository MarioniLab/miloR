######## Methods ########
#' @export
#' @aliases Milo
#' @describeIn Milo get graph
setMethod("graph", "Milo", function(x) {
    if(length(x@graph)){
        x@graph[[1]]
    } else{
        warning("Graph not set")
        list()
        }
    })

#' @export
#' @aliases Milo
#' @describeIn Milo set graph
setMethod("graph<-", "Milo", function(x, value){
    x@graph <- list("graph"=value)
    validObject(x)
    x
    })


#' @export
#' @aliases Milo
#' @describeIn Milo get nhoodDistances
setMethod("nhoodDistances", "Milo", function(x) x@nhoodDistances)

#' @export
#' @aliases Milo
#' @describeIn Milo set nhoodDistances
setMethod("nhoodDistances<-", "Milo", function(x, value){
    if(!any(class(value) %in% c("dgCMatrix"))){
        x@nhoodDistances <- as(value, "dgCMatrix")
    } else{
        x@nhoodDistances <- value
    }

    validObject(x)
    x
})


#' @export
#' @aliases Milo
#' @describeIn Milo get nhoods
setMethod("nhoods", "Milo", function(x) x@nhoods)

#' @export
#' @aliases Milo
#' @describeIn Milo set nhoods
setMethod("nhoods<-", "Milo", function(x, value){
    x@nhoods <- value
    validObject(x)
    x
})


#' @export
#' @aliases Milo
#' @describeIn Milo get nhoodCounts
setMethod("nhoodCounts", "Milo", function(x) x@nhoodCounts)

#' @export
#' @aliases Milo
#' @describeIn Milo set nhoodCounts
setMethod("nhoodCounts<-", "Milo", function(x, value){
    x@nhoodCounts <- value
    validObject(x)
    x
})


#' @export
#' @aliases Milo
#' @describeIn Milo get nhoodIndex
setMethod("nhoodIndex", "Milo", function(x) x@nhoodIndex)

#' @export
#' @aliases Milo
#' @describeIn Milo set nhoodIndex
setMethod("nhoodIndex<-", "Milo", function(x, value){
    x@nhoodIndex <- value
    validObject(x)
    x
})


#' @export
#' @aliases Milo
#' @describeIn Milo get nhoodExpression
setMethod("nhoodExpression", "Milo", function(x) x@nhoodExpression)

#' @export
#' @aliases Milo
#' @describeIn Milo set nhoodExpression
setMethod("nhoodExpression<-", "Milo", function(x, value){
    x@nhoodExpression <- value
    validObject(x)
    x
})


#' @export
#' @aliases Milo
#' @describeIn Milo get nhoodReducedDim
setMethod("nhoodReducedDim", "Milo", function(x, value="PCA") {
    x@nhoodReducedDim[[value]]
    })

#' @export
#' @aliases Milo
#' @describeIn Milo set nhoodReducedDim
setMethod("nhoodReducedDim<-", "Milo", function(x, value, rdim="PCA"){
    x@nhoodReducedDim[[rdim]] <- value
    validObject(x)
    x
})


#' @importFrom S4Vectors coolcat
#' @importFrom methods callNextMethod
.milo_show <- function(object) {
    callNextMethod()
    coolcat("nhoods dimensions(%d): %s\n", length(object@nhoods))
    coolcat("nhoodCounts dimensions(%d): %s\n", dim(object@nhoodCounts))
    coolcat("nhoodDistances dimensions(%d): %s\n", dim(object@nhoodDistances))
    coolcat("graph names(%d): %s\n", names(object@graph))
    coolcat("nhoodIndex names(%d): %s\n", length(object@nhoodIndex))
    coolcat("nhoodExpression dimension(%d): %s\n", dim(object@nhoodExpression))
    coolcat("nhoodReducedDim names(%d): %s\n", names(object@nhoodReducedDim))

    sink(file="/dev/null")
    gc()
    sink(file=NULL)
}

#' @export
#' @aliases Milo
#' @describeIn Milo show method
#' @import methods
setMethod("show", "Milo", .milo_show)
