context("Testing testNhoods function")
library(miloR)

### Set up a mock data set using simulated data
library(SingleCellExperiment)
library(scran)
library(scater)
library(irlba)
library(MASS)
library(mvtnorm)
library(BiocParallel)

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

test_that("Wrong input gives errors", {
    # missing neighbourhood counts
    expect_error(testNhoods(sim1.mylo, design=~Condition,
                                    design.df=sim1.meta),
                 "Neighbourhood counts missing - please run countCells first")

    # count cells
    sim1.mylo <- countCells(sim1.mylo, samples="Sample", meta.data=meta.df)
    expect_error(testNhoods(nhoodCounts(sim1.mylo), design=~Condition,
                            design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]),
                 "Unrecognised input type - must be of class Milo")

    # missing reduced dimension slot
    expect_error(testNhoods(sim1.mylo, design=~Condition,
                            design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ],
                            reduced.dim="blah"),
                 "is not found in reducedDimNames")

})

sim1.mylo <- countCells(sim1.mylo, samples="Sample", meta.data=meta.df)

test_that("Discordant dimension names gives an error", {
    design.matrix <- model.matrix(~Condition, data=sim1.meta)
    rownames(design.matrix) <- c(1:nrow(design.matrix))
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                             design.df=sim1.meta)),
                                  "Design matrix and model matrix rownames are not a subset")

    # warning if only a subset are present, or order is wrong
    design.matrix <- model.matrix(~Condition, data=sim1.meta)
    rownames(design.matrix) <- rownames(sim1.meta)

    set.seed(42)
    design.matrix <- design.matrix[sample(rownames(design.matrix)), ]
    expect_warning(testNhoods(sim1.mylo, design=design.matrix,
                              design.df=sim1.meta),
                   "Sample names in design matrix and nhood counts are not matched. Reordering")

})

test_that("Discordant dimensions between input and design gives an error", {
    add.meta <- sim1.meta[c(1:5), ]
    rownames(add.meta) <- paste0(rownames(add.meta), "_add")
    big.meta <- rbind.data.frame(sim1.meta, add.meta)
    design.matrix <- model.matrix(~Condition, data=big.meta)

    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                                     design.df=sim1.meta)),
                 "Design matrix and model matrix are not the same dimensionality")
})


test_that("Concordant dimensions between input and output", {
    in.rows <- nrow(nhoodCounts(sim1.mylo))
    out.rows <- nrow(testNhoods(sim1.mylo, design=~Condition,
                                        design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    expect_identical(in.rows, out.rows)
})


test_that("Identical results are produced with identical input", {
    # run for each weighting scheme
    kd.ref1 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    kd.ref2 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(kd.ref1, kd.ref2)

    # mean nhood distance
    nd.ref1 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="neighbour-distance",
                          design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    nd.ref2 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="neighbour-distance",
                          design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(nd.ref1, nd.ref2)
})

test_that("testNhoods produces reproducible results with equivalent input", {
    # same input, but different format should give the same results
    design.matrix <- model.matrix(~Condition, data=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    mod.ref <- suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                                   design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    form.ref <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition,
                                                    design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    # test for each column of output
    expect_equal(mod.ref$Pvalue, form.ref$Pvalue)
    expect_equal(mod.ref$SpatialFDR, form.ref$SpatialFDR)
    expect_equal(mod.ref$logFC, form.ref$logFC)
    expect_equal(mod.ref$`F`, form.ref$`F`)
    expect_equal(mod.ref$FDR, form.ref$FDR)
})

test_that("Filtering nhoods provides reproducible results", {
    require(Matrix)
    exp.nh <- sum(Matrix::rowMeans(nhoodCounts(sim1.mylo)) >= 5)
    out.da <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                 min.mean=5,
                                 design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(nrow(out.da), exp.nh)

    kd.ref1 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  min.mean=5,
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    kd.ref2 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  min.mean=5,
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(kd.ref1, kd.ref2)
})

test_that("Model contrasts provide expected results", {
    design.matrix <- model.matrix(~0 + Condition, data=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    cont.ref <- suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                                    model.contrasts=c("ConditionB-ConditionA"),
                                                    design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    form.ref <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition,
                                                    design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))

    # test for each column of output
    expect_equal(cont.ref$Pvalue, form.ref$Pvalue)
    expect_equal(cont.ref$SpatialFDR, form.ref$SpatialFDR)
    expect_equal(cont.ref$logFC, form.ref$logFC)
    expect_equal(cont.ref$`F`, form.ref$`F`)
    expect_equal(cont.ref$FDR, form.ref$FDR)
})

test_that("Providing a subset model.matrix is reproducible", {
    require(Matrix)
    set.seed(42)
    subset.samples <- sample(rownames(sim1.meta))
    exp.nh <- sum(Matrix::rowMeans(nhoodCounts(sim1.mylo)[, subset.samples]) >= 1)
    out.da <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                          min.mean=1,
                                          design.df=sim1.meta[subset.samples, ]))
    expect_identical(nrow(out.da), exp.nh)

    kd.ref1 <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                           min.mean=1,
                                           design.df=sim1.meta[subset.samples, ]))
    kd.ref2 <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                           min.mean=1,
                                           design.df=sim1.meta[subset.samples, ]))
    expect_identical(kd.ref1, kd.ref2)
})

sim1.meta$Condition_num <- paste0("Condition_num", c(1, 1, 1, 0, 0, 0))
sim1.meta$Replicate_num <- paste0("Replicate_num", c(1, 2, 3, 1, 2, 3))
sim1.meta$Replicate2 <- paste0("Replicate2", c(1, 2, 1, 2, 1, 2))

test_that("Singular Hessians are detectable and fail appropriately", {
    set.seed(42)
    # having a singular Hessian depends on some of the staring values <- this test needs to
    # be reproducible and not depend on setting a specific seed. The easiest way might be to have
    # a variance component that is effectively 0.

    # collinear fixed and random effects
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=~Condition + (1|Condition),
                            design.df=sim1.meta, glmm.solver="Fisher", fail.on.error=TRUE)),
                 "Coefficients Hessian is computationally singular")
})

test_that("Invalid formulae give expected errors", {
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=~Condition + (50|Condition),
                                             design.df=sim1.meta, glmm.solver="Fisher")),
                 "is an invalid formula for random effects")
})

test_that("NA or Inf cell sizes causes the expected errors", {
    cell.sizes.na <- colSums(nhoodCounts(sim1.mylo))
    cell.sizes.na[1] <- NA
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=~Condition,
                                             design.df=sim1.meta,
                                             cell.sizes=cell.sizes.na)),
                 "NA or Infinite values found in cell\\.sizes")

    cell.sizes.inf <- colSums(nhoodCounts(sim1.mylo))
    cell.sizes.inf[1] <- Inf
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=~Condition,
                                             design.df=sim1.meta,
                                             cell.sizes=cell.sizes.inf)),
                                  "NA or Infinite values found in cell\\.sizes")
})


