context("Testing nhood grouping function")
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

# need to through in random 0's for realism
block1.drop <- matrix(sapply(1:1000, FUN=function(X) rbinom(n=500, size=1, prob=0.05) == 1), nrow=nrow(block1.gex), byrow=TRUE)
block1.gex[block1.drop] <- 0

block2.cells <- 500
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block2.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block2.eigens <- block2.eigens[order(block2.eigens)]
block2.p <- qr.Q(qr(matrix(rnorm(block2.cells^2, mean=3, sd=0.01), block2.cells)))
block2.sigma <- crossprod(block2.p, block2.p*block2.eigens)
block2.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block2.cells, mean=4, sd=0.01), sigma=block2.sigma))

# need to through in random 0's for realism
block2.drop <- matrix(sapply(1:1000, FUN=function(X) rbinom(n=500, size=1, prob=0.05) == 1), nrow=nrow(block2.gex), byrow=TRUE)
block2.gex[block2.drop] <- 0

block3.cells <- 650
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block3.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block3.eigens <- block3.eigens[order(block3.eigens)]
block3.p <- qr.Q(qr(matrix(rnorm(block3.cells^2, mean=4.5, sd=0.01), block3.cells)))
block3.sigma <- crossprod(block3.p, block3.p*block3.eigens)
block3.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block3.cells, mean=5, sd=0.01), sigma=block3.sigma))

# need to through in random 0's for realism
block3.drop <- matrix(sapply(1:1000, FUN=function(X) rbinom(n=500, size=1, prob=0.05) == 1), nrow=nrow(block3.gex), byrow=TRUE)
block3.gex[block3.drop] <- 0

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
rownames(meta.df) <- colnames(sim1.gex)

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

sim1.meta <- data.frame("Condition"=c(rep("A", 3), rep("B", 3)),
                        "Replicate"=rep(c("R1", "R2", "R3"), 2))
sim1.meta$Sample <- paste(sim1.meta$Condition, sim1.meta$Replicate, sep="_")
rownames(sim1.meta) <- sim1.meta$Sample
sim1.mylo <- countCells(sim1.mylo, samples="Sample", meta.data=meta.df)
sim1.res <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                       design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])

test_that("Incorrect input gives the expected errors", {
    expect_error(groupNhoods(matrix(0L, nrow=nrow(sim1.mylo), ncol=ncol(sim1.mylo))),
                 "Unrecognised input type")

    fake.res <- sim1.res
    fake.res$SpatialFDR <- 1
    expect_error(groupNhoods(sim1.mylo, fake.res),
                 "No DA neighbourhoods found")
})

test_that("Output is correct type", {
    
    expect_length(groupNhoods(sim1.mylo, sim1.res)$NhoodGroup, nrow(sim1.res))
    
    })


test_that("subsetting gives the expected result", {
    full.out <- suppressWarnings(groupNhoods(sim1.mylo, sim1.res,
                                                  compute.new=TRUE, subset.nhoods = c(1:10)))
    expect_true(any(is.na(full.out$NhoodGroup)))
    expect_equal(sum(!is.na(full.out$NhoodGroup)), 10)
    full.out <- suppressWarnings(groupNhoods(sim1.mylo, sim1.res,
                                             compute.new=TRUE, subset.nhoods = sim1.res$logFC < 0.1))
    expect_equal(sum(!is.na(full.out$NhoodGroup)), sum(sim1.res$logFC < 0.1))
    full.out <- suppressWarnings(groupNhoods(sim1.mylo, sim1.res,
                                             compute.new=TRUE, subset.nhoods = colnames(nhoods(sim1.mylo))[1:10]))
    expect_equal(sum(!is.na(full.out$NhoodGroup)), 10)
    
})


# test_that("Every nhood is assigned a group", {
#     full.out <- suppressWarnings(groupNhoods(sim1.mylo, sim1.res,
#                                                   compute.new=TRUE))
#     
#     
#     expect_false(any(is.na(full.out$NhoodGroup)))
# })

test_that("Setting max.lfc.delta increases the number of groups", {
    delta.out <- suppressWarnings(groupNhoods(sim1.mylo, sim1.res, max.lfc.delta = 0.3,
                                             compute.new=TRUE))
    full.out <- suppressWarnings(groupNhoods(sim1.mylo, sim1.res,
                                             compute.new=TRUE))
      
    expect_lt(length(unique(full.out$NhoodGroup)), length(unique(delta.out$NhoodGroup)))
})



