library(methods)
library(SingleCellExperiment)

## allow matrices to be dense or sparse on instantiation
## cast matrices to be sparse where possible





######## Methods ########

### graph slot
setGeneric("graph", function(x) standardGeneric("graph"))
setGeneric("graph<-", function(x, value) standardGeneric("graph<-"))

setMethod("graph", "Milo", function(x) x@graph)
setMethod("graph<-", "Milo", function(x, value){
    x@graph <- value
    validObject(x)
    x
    })

### adjacency matrix slot
setGeneric("adjacency", function(x) standardGeneric("adjacency"))
setGeneric("adjacency<-", function(x, value) standardGeneric("adjacency<-"))

setMethod("adjacency", "Milo", function(x) x@adjacency)
setMethod("adjacency<-", "Milo", function(x, value){
    x@adjacency <- value
    validObject(x)
    x
})

### distance matrix slot
setGeneric("distance", function(x) standardGeneric("distance"))
setGeneric("distance<-", function(x, value) standardGeneric("distance<-"))

setMethod("distance", "Milo", function(x) x@distance)
setMethod("distance<-", "Milo", function(x, value){
    x@adjacency <- value
    validObject(x)
    x
})


### neighbourhoodCounts matrix slot
setGeneric("neighbourhoodCounts", function(x) standardGeneric("neighbourhoodCounts"))
setGeneric("neighbourhoodCounts<-", function(x, value) standardGeneric("neighbourhoodCounts<-"))

setMethod("neighbourhoodCounts", "Milo", function(x) x@neighbourhoodCounts)
setMethod("neighbourhoodCounts<-", "Milo", function(x, value){
    x@neighbourhoodCounts <- value
    validObject(x)
    x
})


#' @importFrom S4Vectors coolcat
#' @importFrom methods callNextMethod
#'
.milo_show <- function(object) {
    methods::callNextMethod()
    S4Vectors::coolcat("neighbourhoodCounts names(%d): %s\n", names(object@neighbourhoodCounts))
    S4Vectors::coolcat("distance names(%d): %s\n", names(object@distance))
    S4Vectors::coolcat("adjacency names(%d): %s\n", names(object@adjacency))
    S4Vectors::coolcat("graph names(%d): %s\n", names(object@graph))

}

setMethod("show", "Milo", .milo_show)
