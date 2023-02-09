context("Testing plotting functions")
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

#test for DA
sim1.mylo <- countCells(sim1.mylo, samples="Sample", meta.data=meta.df)
sim1.da.res <- testNhoods(sim1.mylo, design = ~ Condition, design.df = sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])

## Tests for plotNhoodExpressionDA ##

test_that("Incorrect input features produce proper errors", {
    # input features not in milo object
    expect_error(plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = "GeneA"),
                 "Some features are not in rownames(x)", fixed=TRUE)
    expect_error(plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = NA),
               "Some features are not in rownames(x)", fixed=TRUE)
    sim1.mylo <- calcNhoodExpression(sim1.mylo)
    expect_error(plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = c()),
                 "features is empty", fixed=TRUE)
    expect_error(plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = c("blah")),
                 "Some features are not in rownames(x)", fixed=TRUE)
})

test_that("calcNhoodExpression is run within the function only if needed", {
  feats <- paste0("Gene", 1:100)
  sim1.mylo.2 <- calcNhoodExpression(sim1.mylo, subset.row = feats)
  expect_warning(plotNhoodExpressionDA(sim1.mylo.2, sim1.da.res, features = c("Gene101", "Gene102", "Gene103")))
  expect_warning(plotNhoodExpressionDA(sim1.mylo.2, sim1.da.res, features = c("Gene1", "Gene2", "Gene103")))
  expect_silent(plotNhoodExpressionDA(sim1.mylo.2, sim1.da.res, features = c("Gene1", "Gene2", "Gene3")))
  })

sim1.mylo <- calcNhoodExpression(sim1.mylo)

test_that("Subsetting produces the expected number of neighbourhoods", {
  max.length <- nrow(sim1.da.res)
  subset.length <- max.length - 10
  p <- plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = c("Gene101", "Gene102"),
                             subset.nhoods = c(1:subset.length))
  expect_equal(length(unique(p$data$Nhood)), subset.length)
  })

test_that("Different input types produce the same subsetting", {
  subset_numeric <- 1:10
  subset_logical <- c(rep(TRUE, 10), rep(FALSE, ncol(nhoods(sim1.mylo))-10))
  p <- suppressWarnings(plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = c("Gene101", "Gene102"),
                                              subset.nhoods = subset_numeric))
  p1 <- suppressWarnings(plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = c("Gene101", "Gene102"),
                                               subset.nhoods = subset_logical))

  # extract the components of the 2 plots
  expect_equal(str(p), str(p1))
})

test_that("The order of features is maintained if cluster_features=FALSE", {
  p_feats <- paste0("Gene", 400:300)
  p <- plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = p_feats,
                             cluster_features = FALSE)
  expect_true(all(levels(p$data[["feature"]]) == p_feats))
  p1 <- plotNhoodExpressionDA(sim1.mylo, sim1.da.res, features = p_feats,
                             cluster_features = TRUE)
  expect_false(all(levels(p1$data[["feature"]]) == p_feats))
  })


test_that("Incorrect input produce expected error in plotNhoodCounts", {
  expect_error(plotNhoodCounts(x=sim1.sce,
                               subset.nhoods=c("1", "2"),
                               design.df=sim1.meta,
                               condition="Condition"),
               "Unrecognised input type - must be of class Milo",
               fixed=TRUE)

  tmp.milo = Milo(sim1.sce)
  expect_error(plotNhoodCounts(x=tmp.milo,
                               subset.nhoods=c("1", "2"),
                               design.df=sim1.meta,
                               condition="Condition"),
               "No neighbourhoods found. Please run makeNhoods() first.",
               fixed=TRUE)

  tmp.milo = buildGraph(tmp.milo, k = 30, d = 3)
  tmp.milo = makeNhoods(tmp.milo)
  expect_error(plotNhoodCounts(x=tmp.milo,
                               subset.nhoods=c("1", "2"),
                               design.df=sim1.meta,
                               condition="Condition"),
               "Neighbourhood counts missing - please run countCells() first",
               fixed=TRUE)

  tmp.mdata <- sim1.meta
  rownames(tmp.mdata)<-NULL
  expect_error(plotNhoodCounts(x=sim1.mylo,
                               subset.nhoods=c("1", "2"),
                               design.df=tmp.mdata,
                               condition="Condition"),
               "The design.df has to be of type data.frame with rownames that correspond to the samples.",
               fixed=TRUE)

  expect_error(plotNhoodCounts(x=sim1.mylo,
                               subset.nhoods=c("1", "2"),
                               design.df=sim1.meta,
                               condition="Batch"),
               "Condition of interest has to be a column in the design matrix",
               fixed=TRUE)

  expect_error(plotNhoodCounts(x=sim1.mylo,
                               subset.nhoods=c("1","2","a34"),
                               design.df=sim1.meta,
                               condition="Condition"),
               paste0("Specified subset.nhoods do not exist - ",
                      "these should either be an integer or character vector corresponding to row names in nhoodCounts(x) ",
                      "or a logical vector with length nrow(nhoodCounts(x))."),
               fixed=TRUE)

  expect_error(plotNhoodCounts(x=sim1.mylo,
                               subset.nhoods=c(TRUE, FALSE, FALSE, TRUE),
                               design.df=sim1.meta,
                               condition="Condition"),
               "Length of the logical vector has to match number of rows in nhoodCounts(x)",
               fixed=TRUE)


})

test_that("Data is correctly reshaped and plotted in plotNhoodCounts",{
  nhoods_of_interest = c("1", "2")

  p <- plotNhoodCounts(x=sim1.mylo,
                  subset.nhoods=nhoods_of_interest,
                  design.df=sim1.meta,
                  condition="Condition")

  # check if we have the expected number of rows in our ggplot
  expect_equal(nrow(p$data), length(nhoods_of_interest)*length(unique(sim1.meta$Sample)))
})

test_that("Same result regardless of the type of nhood vector in plotNhoodCounts",{
  nhoods_chr_vector = c("1", "2")
  p_chr <- plotNhoodCounts(x=sim1.mylo,
                       subset.nhoods=nhoods_chr_vector,
                       design.df=sim1.meta,
                       condition="Condition")

  nhoods_num_vector <- c(1,2)
  p_num <- plotNhoodCounts(x=sim1.mylo,
                           subset.nhoods=nhoods_num_vector,
                           design.df=sim1.meta,
                           condition="Condition")

  nhoods_logi_vector <- c(TRUE, TRUE, rep(FALSE, nrow(nhoodCounts(sim1.mylo))-2))

  p_logi <- plotNhoodCounts(x=sim1.mylo,
                           subset.nhoods=nhoods_logi_vector,
                           design.df=sim1.meta,
                           condition="Condition")

  # all ggplot objects should contain the exact same data regardless of vector type.
  expect_identical(p_chr$data, p_num$data)
  expect_identical(p_num$data, p_logi$data)
})





