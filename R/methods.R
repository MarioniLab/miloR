######## Methods ########
#' @export
#' @describeIn Milo get graph
setMethod("graph", "Milo", function(x) x@graph[[1]])

#' @export
#' @describeIn Milo set graph
setMethod("graph<-", "Milo", function(x, value){
    x@graph <- list("graph"=value)
    validObject(x)
    x
    })


#' @export
#' @describeIn Milo get neighbourDistances
setMethod("neighbourDistances", "Milo", function(x) x@neighbourDistances)

#' @export
#' @describeIn Milo set neighbourDistances
setMethod("neighbourDistances<-", "Milo", function(x, value){
    x@neighbourDistances <- value
    validObject(x)
    x
})


#' @export
#' @describeIn Milo get neighbourhoods
setMethod("neighbourhoods", "Milo", function(x) x@neighbourhoods)

#' @export
#' @describeIn Milo set neighbourhoods
setMethod("neighbourhoods<-", "Milo", function(x, value){
    x@neighbourhoods <- value
    validObject(x)
    x
})


#' @export
#' @describeIn Milo get neighbourhoodCounts
setMethod("neighbourhoodCounts", "Milo", function(x) x@neighbourhoodCounts)

#' @export
#' @describeIn Milo set neighbourhoodCounts
setMethod("neighbourhoodCounts<-", "Milo", function(x, value){
    x@neighbourhoodCounts <- value
    validObject(x)
    x
})


#' @importFrom S4Vectors coolcat
#' @importFrom methods callNextMethod
.milo_show <- function(object) {
    callNextMethod()
    coolcat("neighbourhoods dimensions(%d): %s\n", length(object@neighbourhoods))
    coolcat("neighbourhoodCounts dimensions(%d): %s\n", dim(object@neighbourhoodCounts))
    coolcat("neighbourDistances dimensions(%d): %s\n", dim(object@neighbourDistances))
    coolcat("graph names(%d): %s\n", names(object@graph))

    sink(file="/dev/null")
    gc()
    sink(file=NULL)
}

#' @export
#' @describeIn Milo show method
#' @import methods
setMethod("show", "Milo", .milo_show)
