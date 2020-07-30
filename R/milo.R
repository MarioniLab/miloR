library(methods)

# Do I really want the SCE as the backbone?

setClass("Milo",
         contains = "SingleCellExperiment",
         slots=c(
             graph = "ANY", # this should be NA or an igraph object
             adjacency = "ANY", # this should be NA or a matrix
             distance = "ANY", # this should be NA or a matrix
             counts = "ANY" # this should be NA or a matrix
             ),
         prototype = list(
             graph = NA_real_,
             adjacency = NA_real_,
             distance = NA_real_,
             counts = NA_real_
             )
         )

## class helper
Milo <- function(graph=NA_real_, adjacency=NA_real_, distance=NA_real_, counts=NA_real_){
    new("Milo", graph=graph, adjacency=adjacency, distance=distance, counts=counts)
}

## class validator
setValidity("Milo", function(object){
    if (class(object@counts) != "matrix"){
        "@counts must be a matrix format"
    } else{
        TRUE
    }

    if (class(object@graph) != "numeric"){
        "@graph must be of type numeric"
    } else{
        TRUE
    }
})





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


### counts matrix slot
setGeneric("counts", function(x) standardGeneric("counts"))
setGeneric("counts<-", function(x, value) standardGeneric("counts<-"))

setMethod("counts", "Milo", function(x) x@counts)
setMethod("counts<-", "Milo", function(x, value){
    x@counts <- value
    validObject(x)
    x
})


# ## define show method
# setMethod("show", "Milo", function(object){
#     cat(is(object)[[1]], "\n",
#     " counts: ", class(object@counts), "\n",
#     " graph: ", class(object@graph), "\n",
#     sep=""
#     )
# })






### tests ###
mylo <- Milo(graph=5)

graph(mylo) <- 341
graph(mylo)
counts(mylo)
distance(mylo)
adjacency(mylo)

# create an invalid object
wrong <- Milo(graph="blah")




