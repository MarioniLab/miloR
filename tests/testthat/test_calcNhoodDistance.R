context("Test calcNhoodDistance function")
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


block2.cells <- 500
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block2.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block2.eigens <- block2.eigens[order(block2.eigens)]
block2.p <- qr.Q(qr(matrix(rnorm(block2.cells^2, mean=4, sd=0.01), block2.cells)))
block2.sigma <- crossprod(block2.p, block2.p*block2.eigens)
block2.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block2.cells, mean=4, sd=0.01), sigma=block2.sigma))


block3.cells <- 650
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

# build a graph - this can take a while for large graphs - will need to play
# around with the parallelisation options
sim1.mylo <- buildGraph(sim1.mylo, k=21, d=30)

# define neighbourhoods - this is slow for large data sets
# how can this be sped up? There are probably some parallelisable steps
sim1.mylo <- makeNhoods(sim1.mylo, k=21, prop=0.1, refined=TRUE,
                        d=30,
                        reduced_dims="PCA")

test_that("Incorrect input gives expected errors", {
    expect_error(calcNhoodDistance("blah", d=30), "Input is not a valid Milo object")

    reducedDim(sim1.mylo, "PCA") <- NULL
    expect_error(suppressWarnings(calcNhoodDistance(sim1.mylo, d=30, reduced.dim="blah")),
                                  "not found in reducedDim slot")
})

test_that("The populated slot is a list of sparse matrices", {
    return.obj <- nhoodDistances(calcNhoodDistance(sim1.mylo, d=30))
    expect_identical(class(return.obj), "list")

    expect_true(all(unlist(lapply(return.obj, function(X) class(X) %in% c("dgCMatrix")))))
})

test_that("The list of matrices correspond to the same sizes as the neighbourhoods", {
    return.cols <- lapply(nhoodDistances(calcNhoodDistance(sim1.mylo, d=30)),
                          function(X) ncol(X))
    return.rows <- lapply(nhoodDistances(calcNhoodDistance(sim1.mylo, d=30)),
                          function(X) nrow(X))
    expect_identical(return.cols, return.rows)

    nhood.sizes <- colSums(nhoods(sim1.mylo))
    expect_equal(nhood.sizes, unlist(return.cols))
    expect_equal(nhood.sizes, unlist(return.rows))
})

test_that("calcNhoodDistance produces identical output", {
    .test_calc_distance <- function(in.x){

        dist.list <- list()
        for(i in seq_along(1:nrow(in.x))){
            i.diff <- t(apply(in.x, 1, FUN=function(P) P - in.x[i, ]))
            i.dist <- sqrt(rowSums(i.diff**2))
            dist.list[[paste0(i)]] <- list("rowIndex"=rep(i, nrow(in.x)), "colIndex"=c(1:length(i.dist)),
                                           "dist"=i.dist)
        }

        dist.df <- do.call(rbind.data.frame, dist.list)
        out.dist <- sparseMatrix(i=dist.df$rowIndex, j=dist.df$colIndex, x=dist.df$dist,
                                 dimnames=list(rownames(in.x), rownames(in.x)),
                                 repr="T")
        return(out.dist)
    }

    nhood.dists <- sapply(colnames(nhoods(sim1.mylo)),
                          function(X) .test_calc_distance(reducedDim(sim1.mylo, "PCA")[nhoods(sim1.mylo)[, X] > 0,
                                                                                       c(1:30),drop=FALSE]))
    names(nhood.dists) <- nhoodIndex(sim1.mylo)
    return.obj  <- nhoodDistances(calcNhoodDistance(sim1.mylo, d=30))
    expect_identical(length(nhood.dists), length(return.obj))
    expect_identical(names(nhood.dists), names(return.obj))
    expect_identical(lapply(nhood.dists, dim), lapply(return.obj, dim))
})

