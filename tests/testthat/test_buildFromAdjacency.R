context("Test buildFromAdjacency function")
library(miloR)

### Set up a mock data set using simulated data
library(SingleCellExperiment)
library(scran)
library(scater)
library(irlba)
library(MASS)
library(mvtnorm)
library(miloR)

set.seed(42)
r.n <- 1000
n.dim <- 50
block1.cells <- 500
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block1.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block1.eigens <- block1.eigens[order(block1.eigens)]
block1.p <- qr.Q(qr(matrix(rnorm(block1.cells^2, mean=4, sd=0.01), block1.cells)))
block1.sigma <- crossprod(block1.p, block1.p*block1.eigens)
block1.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block1.cells, mean=2, sd=0.01), sigma=block1.sigma))


block2.cells <- 400
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block2.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block2.eigens <- block2.eigens[order(block2.eigens)]
block2.p <- qr.Q(qr(matrix(rnorm(block2.cells^2, mean=4, sd=0.01), block2.cells)))
block2.sigma <- crossprod(block2.p, block2.p*block2.eigens)
block2.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block2.cells, mean=4, sd=0.01), sigma=block2.sigma))


block3.cells <- 200
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block3.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block3.eigens <- block3.eigens[order(block3.eigens)]
block3.p <- qr.Q(qr(matrix(rnorm(block3.cells^2, mean=4, sd=0.01), block3.cells)))
block3.sigma <- crossprod(block3.p, block3.p*block3.eigens)
block3.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block3.cells, mean=5, sd=0.01), sigma=block3.sigma))


sim1.gex <- do.call(cbind, list("b1"=block1.gex, "b2"=block2.gex, "b3"=block3.gex))
colnames(sim1.gex) <- paste0("Cell", 1:ncol(sim1.gex))
rownames(sim1.gex) <- paste0("Gene", 1:nrow(sim1.gex))
sim1.pca <- prcomp_irlba(t(sim1.gex), n=50, scale.=TRUE, center=TRUE)

sim1.sce <- SingleCellExperiment(assays=list(logcounts=sim1.gex),
                                 reducedDims=list("PCA"=sim1.pca$x))
attr(reducedDim(sim1.sce, "PCA"), "rotation") <- sim1.pca$rotation

sim1.mylo <- Milo(sim1.sce)
sim1.mylo <- buildGraph(sim1.mylo, k=21, d=30)
sim1.graph <- miloR::graph(sim1.mylo)
sim1.adj <- as_adjacency_matrix(sim1.graph)
test.graph <- buildFromAdjacency(as(sim1.adj, "matrix"), k=21)

test_that("Incorrect input generates expected errors", {
    expect_error(buildFromAdjacency(list()), "Input 'x' is not a recognisable matrix format")
})

test_that("Casting to sparse matrix generates message", {
    expect_message(buildFromAdjacency(as(sim1.adj, "matrix"), k=21),
                   "Casting to sparse matrix format")
})

test_that("Inferring k works on row-wise adjacency matrix", {
    r <- 1000
    c <- 1000
    k <- 35
    m <- floor(matrix(runif(r*c), r, c))
    for(i in seq_along(1:r)){
        m[i, sample(1:c, size=k)] <- 1
    }

    expect_message(buildFromAdjacency(as(m, "dgTMatrix")),
                   "Inferring k from matrix")
})

test_that("Inferring k works on column-wise adjacency matrix", {
    r <- 1000
    c <- 1000
    k <- 35
    m <- floor(matrix(runif(r*c), r, c))
    for(i in seq_along(1:c)){
        m[sample(1:r, size=k), i] <- 1
    }

    expect_warning(buildFromAdjacency(as(m, "dgTMatrix")),
                   "Row sums are not all equal")
})

test_that("Error produced when matrix is undirected and no k is provided", {
    expect_error(suppressWarnings(buildFromAdjacency(sim1.adj)),
                 "Cannot infer k from matrix")
})

test_that("Inputting a non-square binary matrix result generates the expected error", {
    r <- 1000
    c <- 900
    k <- 35
    m <- floor(matrix(runif(r*c), r, c))
    for(i in seq_along(1:c)){
        m[sample(1:r, size=k), i] <- 1
    }

    expect_error(buildFromAdjacency(as(m, "dgTMatrix"), k=35),
                   "Input matrix is binary but not square")
})

test_that("Providing a distance matrix generates a message", {
    sim1.dist <- nhoodDistances(sim1.mylo)
    expect_message(suppressWarnings(buildFromAdjacency(sim1.dist, k=21)),
                   "Adding nhoodDistances to Milo object")
})

test_that("buildFrom adjacency is reproducible", {
    # is the inferred graph the same as the original?
    # using identical_graph doesn't seem to work, despite these graphs being identical
    # check characteristics of the graph to be sure
    adj.graph <- miloR::graph(buildFromAdjacency(sim1.adj, k=21))
    adj.degree <- degree(adj.graph)
    expect_equal(degree(sim1.graph), adj.degree)

    adj.n.vertices <- length(V(adj.graph))
    expect_equal(length(V(sim1.graph)), adj.n.vertices)

    adj.n.edges <- length(E(adj.graph))
    expect_equal(length(E(sim1.graph)), adj.n.edges)
})

test_that("A valid Milo object is created", {
    adj.mylo <- suppressMessages(buildFromAdjacency(sim1.adj, k=21))
    expect_equal(class(graph(adj.mylo)), "igraph")

    sim1.dist <- nhoodDistances(sim1.mylo)
    adj.mylo <- suppressMessages(buildFromAdjacency(sim1.dist, k=21))
    expect_equal(class(graph(adj.mylo)), "igraph")
    expect_true(is(nhoodDistances(adj.mylo), "sparseMatrix"))
    expect_equal(ncol(nhoodDistances(adj.mylo)), length(V(graph(adj.mylo))))
    expect_equal(nrow(nhoodDistances(adj.mylo)), length(V(graph(adj.mylo))))
})





