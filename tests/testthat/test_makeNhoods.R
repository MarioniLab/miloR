context("Testing makeNhoods function")
library(miloR)

### Set up a mock data set using simulated data
library(SingleCellExperiment)
library(scran)
library(scater)
library(irlba)
library(MASS)
library(mvtnorm)

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
sim1.pca <- prcomp_irlba(t(sim1.gex), n=50+1, scale.=TRUE, center=TRUE)

set.seed(42)
block1.cond <- rep("A", block1.cells)
block1.a <- sample(1:block1.cells, size=floor(block1.cells*0.9))
block1.b <- setdiff(1:block1.cells, block1.a)
block1.cond[block1.b] <- "B"

block2.cond <- rep("A", block2.cells)
block2.a <- sample(1:block2.cells, size=floor(block2.cells*0.05))
block2.b <- setdiff(1:block2.cells, block2.a)
block2.cond[block2.b] <- "B"

block3.cond <- rep("A", block3.cells)
block3.a <- sample(1:block3.cells, size=floor(block3.cells*0.5))
block3.b <- setdiff(1:block3.cells, block3.a)
block3.cond[block3.b] <- "B"

meta.df <- data.frame("Block"=c(rep("B1", block1.cells), rep("B2", block2.cells), rep("B3", block3.cells)),
                      "Condition"=c(block1.cond, block2.cond, block3.cond),
                      "Replicate"=c(rep("R1", floor(block1.cells*0.33)), rep("R2", floor(block1.cells*0.33)),
                                    rep("R3", block1.cells-(2*floor(block1.cells*0.33))),
                                    rep("R1", floor(block2.cells*0.33)), rep("R2", floor(block2.cells*0.33)),
                                    rep("R3", block2.cells-(2*floor(block2.cells*0.33))),
                                    rep("R1", floor(block3.cells*0.33)), rep("R2", floor(block3.cells*0.33)),
                                    rep("R3", block3.cells-(2*floor(block3.cells*0.33)))))
colnames(meta.df) <- c("Block", "Condition", "Replicate")
# define a "sample" as teh combination of condition and replicate
meta.df$Sample <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
meta.df$Vertex <- c(1:nrow(meta.df))

sim1.sce <- SingleCellExperiment(assays=list(logcounts=sim1.gex),
                                 reducedDims=list("PCA"=sim1.pca$x))

sim1.mylo <- Milo(sim1.sce)

sim1.graph <- graph(buildGraph(sim1.mylo, k=21, d=30))

test_that("Incorrect input gives informative error", {
    # missing graph
    graph(sim1.mylo) <- list()
    expect_error(makeNhoods(sim1.mylo, d=30), "Not a valid Milo object - graph is missing")

    # pass a graph but no reduced dimensions
    expect_error(makeNhoods(sim1.graph, d=30), "No reduced dimensions matrix provided")

    # unexpected input type
    expect_error(makeNhoods(list(), d=30), paste0("Data format: ", class(list()), " not recognised."))

    # wrong refinement scheme given
    expect_error(makeNhoods(sim1.graph, d = 30, refined = TRUE, refinement_scheme = "NA"), "When refined == TRUE, refinement_scheme must be one of \"reduced_dim\" or \"graph\".")
})

test_that("No errors arise with any of the possible inputs (igraph or milo object, named or unnamed)", {

    sim1.mylo <- buildGraph(sim1.mylo, k=21)
    sim1.graph.vertex <- set_vertex_attr(sim1.graph, "name", index = V(sim1.graph), as.character(colnames(sim1.mylo)))
    sim1.mylo.nonames <- sim1.mylo
    colnames(sim1.mylo.nonames) <- NULL
    reduced_dimensions <- sim1.pca$x
    rownames(reduced_dimensions) <- colnames(sim1.sce)

    ### input is igraph

    # pass a graph with no vertex attributes and no reduced dimensions - should not throw an error for refinement_scheme = graph
    expect_error(makeNhoods(sim1.graph, d=30, refined = TRUE, refinement_scheme = "graph"), regexp = NA)

    # pass a graph with vertex attributes and no reduced dimensions - should not throw an error for refinement_scheme = graph
    expect_error(makeNhoods(sim1.graph.vertex, d=30, refined = TRUE, refinement_scheme = "graph"), regexp = NA)

    # pass a graph with no vertex attributes and reduced dimensions matrix - should not throw an error for refinement_scheme = "reduced_dim"
    expect_error(makeNhoods(sim1.graph, d=30, refined = TRUE, reduced_dims = reduced_dimensions, refinement_scheme = "reduced_dim"), regexp = NA)

    # pass a graph with vertex attributes and reduced dimensions matrix - should not throw an error for refinement_scheme = "reduced_dim"
    expect_error(makeNhoods(sim1.graph.vertex, d=30, refined = TRUE,  reduced_dims = reduced_dimensions, refinement_scheme = "reduced_dim"), regexp = NA)

    # pass a graph with vertex attributes and reduced dimensions matrix - should throw a warning for refinement_scheme = "graph"
    expect_warning(makeNhoods(sim1.graph.vertex, d=30, refined = TRUE,  reduced_dims = reduced_dimensions, refinement_scheme = "graph"),
                   "Ignoring reduced dimensions matrix because refinement_scheme = graph was selected.")

    ### input is milo object

    # pass a milo object with no colnames - should not throw an error with refinement_scheme = "reduced_dim"
    expect_warning(makeNhoods(sim1.mylo.nonames, d=30, refined = TRUE, refinement_scheme = "reduced_dim"), regexp = "Rownames not set")

    # pass a milo object with colnames - should not throw an error with refinement_scheme = "reduced_dim"
    expect_error(makeNhoods(sim1.mylo, d=30, refined = TRUE, refinement_scheme = "reduced_dim"), regexp = NA)

    # pass a milo object with no colnames - should not throw an error with refinement_scheme = "graph"
    expect_error(makeNhoods(sim1.mylo.nonames, d=30, refined = TRUE, refinement_scheme = "graph"), regexp = NA)

    # pass a milo object with colnames - should not throw an error with refinement_scheme = "graph"
    expect_error(makeNhoods(sim1.mylo, d=30, refined = TRUE, refinement_scheme = "graph"), regexp = NA)

})


test_that("Passing d parameter higher than ncols of reducedDims gives warning", {
    sim1.mylo <- buildGraph(sim1.mylo, k=21)
    expect_warning(makeNhoods(sim1.mylo, d=ncol(reducedDim(sim1.mylo, "PCA")) + 1),
                   "Specified d is higher than the total number of dimensions in")
})


test_that("Random sampling returns the same vertex list with a set seed", {
    sim1.mylo <- buildGraph(sim1.mylo, k=21)
    set.seed(101)
    random.vertices <- nhoods(makeNhoods(sim1.mylo, refined=FALSE))
    set.seed(101)
    expect_equal(nhoods(makeNhoods(sim1.mylo, refined=FALSE)), random.vertices)

    # indices should also be the same
    set.seed(101)
    random.indx <- nhoodIndex(makeNhoods(sim1.mylo, refined=FALSE))
    set.seed(101)
    expect_equal(nhoodIndex(makeNhoods(sim1.mylo, refined=FALSE)), random.indx)
})


test_that("Refined sampling strictly returns equal to or fewer neighbourhoods than random sampling", {
    sim1.mylo <- buildGraph(sim1.mylo, k=21)
    set.seed(101)
    random.vertices <- nhoods(makeNhoods(sim1.mylo, refined=FALSE))
    expect_true(ncol(nhoods(makeNhoods(sim1.mylo, refined=TRUE))) <= ncol(random.vertices))
    expect_true(ncol(nhoods(makeNhoods(sim1.mylo, refined = TRUE, refinement_scheme = "graph"))) <= ncol(random.vertices))
})





