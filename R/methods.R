######## Methods ########
#' @export
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
#' @describeIn Milo set graph
setMethod("graph<-", "Milo", function(x, value){
    x@graph <- list("graph"=value)
    validObject(x)
    x
    })


#' @export
#' @describeIn Milo get nhoodDistances
setMethod("nhoodDistances", "Milo", function(x) x@nhoodDistances)

#' @export
#' @describeIn Milo set nhoodDistances
setMethod("nhoodDistances<-", "Milo", function(x, value){
    x@nhoodDistances <- value
    validObject(x)
    x
})


#' @export
#' @describeIn Milo get nhoods
setMethod("nhoods", "Milo", function(x) x@nhoods)

#' @export
#' @describeIn Milo set nhoods
setMethod("nhoods<-", "Milo", function(x, value){
    x@nhoods <- value
    validObject(x)
    x
})


#' @export
#' @describeIn Milo get nhoodCounts
setMethod("nhoodCounts", "Milo", function(x) x@nhoodCounts)

#' @export
#' @describeIn Milo set nhoodCounts
setMethod("nhoodCounts<-", "Milo", function(x, value){
    x@nhoodCounts <- value
    validObject(x)
    x
})


#' @export
#' @describeIn Milo get nhoodIndex
setMethod("nhoodIndex", "Milo", function(x) x@nhoodIndex)

#' @export
#' @describeIn Milo set nhoodIndex
setMethod("nhoodIndex<-", "Milo", function(x, value){
    x@nhoodIndex <- value
    validObject(x)
    x
})


#' @export
#' @describeIn Milo get nhoodExpression
setMethod("nhoodExpression", "Milo", function(x) x@nhoodExpression)

#' @export
#' @describeIn Milo set nhoodExpression
setMethod("nhoodExpression<-", "Milo", function(x, value){
    x@nhoodExpression <- value
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

    sink(file="/dev/null")
    gc()
    sink(file=NULL)
}

#' @export
#' @describeIn Milo show method
#' @import methods
setMethod("show", "Milo", .milo_show)
