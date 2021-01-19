context("Test spatialFDR function")
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
sim1.mylo <- calcNhoodDistance(sim1.mylo, d=30)

sim1.meta <- data.frame("Condition"=c(rep("A", 3), rep("B", 3)),
                        "Replicate"=rep(c("R1", "R2", "R3"), 2))
sim1.meta$Sample <- paste(sim1.meta$Condition, sim1.meta$Replicate, sep="_")
rownames(sim1.meta) <- sim1.meta$Sample

sim1.mylo <- countCells(sim1.mylo, samples="Sample", meta.data=meta.df)

test_that("Incorrect parameter values produce errors", {
    # call outside of testNhoods function
    da.ref <- testNhoods(sim1.mylo, design=~Condition,
                         fdr.weighting="none",
                         design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])

    expect_error(graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                                 graph=miloR::graph(sim1.mylo),
                                 weighting="neighbour-distance",
                                 pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                                 indices=nhoodIndex(sim1.mylo),
                                 distances=nhoodDistances(sim1.mylo),
                                 reduced.dimensions=NULL),
                 "A matrix of reduced dimensions")

    expect_error(graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                                 graph=graph(sim1.mylo),
                                 weighting="k-distance",
                                 pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                                 indices=NULL,
                                 distances=nhoodDistances(sim1.mylo),
                                 reduced.dimensions=NULL),
                 "No neighbourhood indices found")

    expect_error(graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                                 graph=graph(sim1.mylo),
                                 weighting="k-distance",
                                 pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                                 indices=nhoodIndex(sim1.mylo),
                                 distances=NULL,
                                 reduced.dimensions=NULL),
                 "k-distance weighting requires")

    expect_error(graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                                 graph=graph(sim1.mylo),
                                 weighting="blah",
                                 pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                                 indices=nhoodIndex(sim1.mylo),
                                 distances=nhoodDistances(sim1.mylo),
                                 reduced.dimensions=reducedDim(sim1.mylo, "PCA")),
                 "Weighting option not recognised")
})


test_that("Input 'none' produces NAs", {
    da.ref <- testNhoods(sim1.mylo, design=~Condition,
                                 fdr.weighting="none",
                                 design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])

    out.p <- graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                             graph=graph(sim1.mylo),
                             weighting="none",
                             pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                             indices=nhoodIndex(sim1.mylo),
                             distances=nhoodDistances(sim1.mylo),
                             reduced.dimensions=reducedDim(sim1.mylo, "PCA"))
    expect_true(sum(is.na(out.p)) == length(out.p))
})


test_that("graphSpatialFDR produces reproducible results for neighbour-distance weighting", {
    # calculate spatial FDR here using the actual code
    da.ref <- testNhoods(sim1.mylo, design=~Condition,
                                 fdr.weighting="none",
                                 design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])

    nhoods <- nhoods(sim1.mylo)
    graph <- graph(sim1.mylo)
    reduced.dimensions <- reducedDim(sim1.mylo, "PCA")

    pvalues <- da.ref[order(da.ref$Nhood), ]$PValue
    haspval <- !is.na(pvalues)
    if (!all(haspval)) {
        coords <- coords[haspval, , drop=FALSE]
        pvalues <- pvalues[haspval]
    }
    t.connect <- sapply(colnames(nhoods)[haspval],
                        FUN=function(PG) {
                            x.pcs <- reduced.dimensions[nhoods[, PG] > 0, ]
                            x.euclid <- as.matrix(dist(x.pcs))
                            x.distdens <- mean(x.euclid[lower.tri(x.euclid, diag=FALSE)])
                            return(x.distdens)})
    w <- 1/unlist(t.connect)
    w[is.infinite(w)] <- 0

    # Computing a density-weighted q-value.
    o <- order(pvalues)
    pvalues <- pvalues[o]
    w <- w[o]

    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
    adjp <- pmin(adjp, 1)

    if (!all(haspval)) {
        refp <- rep(NA_real_, length(haspval))
        refp[haspval] <- adjp
        adjp <- refp
    }

    func.fdr <- graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                                graph=graph(sim1.mylo),
                                weighting="neighbour-distance",
                                pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                                indices=nhoodIndex(sim1.mylo),
                                distances=nhoodDistances(sim1.mylo),
                                reduced.dimensions=reducedDim(sim1.mylo, "PCA"))

    expect_equal(func.fdr, adjp)
})

test_that("graphSpatialFDR produces reproducible results for k-distance weighting", {
    # calculate spatial FDR here using the actual code - first with PCA, then distances
    da.ref <- testNhoods(sim1.mylo, design=~Condition,
                                 fdr.weighting="none",
                                 design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])

    nhoods <- nhoods(sim1.mylo)
    graph <- graph(sim1.mylo)
    indices <- nhoodIndex(sim1.mylo)
    distances <- nhoodDistances(sim1.mylo)
    reduced.dimensions <- reducedDim(sim1.mylo, "PCA")

    pvalues <- da.ref[order(da.ref$Nhood), ]$PValue
    haspval <- !is.na(pvalues)
    if (!all(haspval)) {
        coords <- coords[haspval, , drop=FALSE]
        pvalues <- pvalues[haspval]
    }
    t.connect <- unlist(lapply(indices, FUN=function(X) max(distances[[as.character(X)]])))
    w <- 1/unlist(t.connect)
    w[is.infinite(w)] <- 0

    # Computing a density-weighted q-value.
    o <- order(pvalues)
    pvalues <- pvalues[o]
    w <- w[o]

    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
    adjp <- pmin(adjp, 1)

    if (!all(haspval)) {
        refp <- rep(NA_real_, length(haspval))
        refp[haspval] <- adjp
        adjp <- refp
    }

    func.fdr <- graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                                graph=graph(sim1.mylo),
                                weighting="k-distance",
                                pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                                indices=nhoodIndex(sim1.mylo),
                                distances=nhoodDistances(sim1.mylo),
                                reduced.dimensions=reducedDim(sim1.mylo, "PCA"))

    expect_equal(func.fdr, adjp)

    # with reduced dims
    suppressWarnings(require(BiocNeighbors))
    pvalues <- da.ref[order(da.ref$Nhood), ]$PValue
    haspval <- !is.na(pvalues)
    if (!all(haspval)) {
        coords <- coords[haspval, , drop=FALSE]
        pvalues <- pvalues[haspval]
    }

    t.connect <- unlist(lapply(indices,
                               FUN=function(X) max(findKNN(reduced.dimensions,
                                                           get.distance=TRUE,
                                                           subset=X, k=21)[["distance"]])))
    w <- 1/unlist(t.connect)
    w[is.infinite(w)] <- 0

    # Computing a density-weighted q-value.
    o <- order(pvalues)
    pvalues <- pvalues[o]
    w <- w[o]

    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
    adjp <- pmin(adjp, 1)

    if (!all(haspval)) {
        refp <- rep(NA_real_, length(haspval))
        refp[haspval] <- adjp
        adjp <- refp
    }

    func.fdr <- graphSpatialFDR(x.nhoods=nhoods(sim1.mylo),
                                graph=graph(sim1.mylo),
                                weighting="k-distance",
                                pvalues=da.ref[order(da.ref$Nhood), ]$PValue,
                                indices=nhoodIndex(sim1.mylo),
                                distances=NULL,
                                reduced.dimensions=reducedDim(sim1.mylo, "PCA"))

    expect_equal(func.fdr, adjp)
})

