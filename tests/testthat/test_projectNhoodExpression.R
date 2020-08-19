context("Testing projectNhoodExpression function")
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
sim1.mylo <- buildGraph(sim1.mylo, k=21, d=30)

test_that("Incorrect input gives informative error", {
    # incorrect input
    expect_error(projectNhoodExpression(graph(sim1.mylo)), "Unrecognised input type")

    # missing nhoods slot
    expect_error(projectNhoodExpression(sim1.mylo), "Not a valid Milo object - nhoods are missing")

    # missing rotation information
    sim1.mylo <- makeNhoods(sim1.mylo, k=21, d=30)
    expect_error(projectNhoodExpression(sim1.mylo), "loading matrix for 'reduced_dim' needs to be stored as attribute")
})

test_that("Warnings produced when expected", {

    # d is too high
    attr(reducedDim(sim1.mylo, "PCA"), "rotation") <- sim1.pca$rotation
    sim1.mylo <- makeNhoods(sim1.mylo, k=21, d=30)
    expect_warning(projectNhoodExpression(sim1.mylo, d=60), "d is higher than the total number of dimensions")

    # nhoodExpression is missing
    expect_warning(projectNhoodExpression(sim1.mylo, d=30), "nhoodExpression slot is empty")
})


test_that("Nhood projection is reproducible", {
    # recreate the expected projection
    attr(reducedDim(sim1.mylo, "PCA"), "rotation") <- sim1.pca$rotation
    sim1.mylo <- makeNhoods(sim1.mylo, k=21, d=30)
    sim1.mylo <- calcNhoodExpression(sim1.mylo)

    # need to check that combining the nhoods and cells doesn't cause any corruption of original values
    loadings <- attr(reducedDim(sim1.mylo, "PCA"), "rotation")[, c(1:30)]
    nhood.project <- t(scale(nhoodExpression(sim1.mylo), scale=TRUE, center=TRUE)) %*% loadings
    rownames(nhood.project) <- paste0("nh_", 1:nrow(nhood.project))

    func.projct <- nhoodReducedDim(projectNhoodExpression(sim1.mylo, d=30))
    nh.projct <- func.projct[grepl(rownames(func.projct), pattern="nh"), ]
    cx.projct <- func.projct[!grepl(rownames(func.projct), pattern="nh"), ]

    expect_identical(nhood.project, nh.projct)
    expect_identical(cx.projct, reducedDim(sim1.mylo, "PCA")[, 1:30])
})
