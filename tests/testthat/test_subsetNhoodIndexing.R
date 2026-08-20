context("Testing nhood indexing consistency under subset.nhoods")
library(miloR)
library(SingleCellExperiment)
library(BiocParallel)

### -------------------------------------------------------------------------
### The GLM and GLMM branches of testNhoods() must label their results with the
### same neighbourhood identifiers.
###
### The GLM branch inherits rownames from the DGEList, which is built on
### nhoodCounts(x)[keep.nh, ], so it returns the original indices with gaps
### where neighbourhoods were dropped. The GLMM branch previously overwrote
### these with 1..n_kept. With subset.nhoods in play the two disagreed, and
### merge(glm.res, glmm.res, by="Nhood") then paired unrelated neighbourhoods -
### which produced an apparent near-zero correlation between GLM and GLMM log
### fold changes on real data where the true correlation is above 0.99.
###
### logCPM is the diagnostic used here because it is a property of the
### neighbourhood and not of the model, so it must agree across the two
### branches for the same neighbourhood regardless of how each one fits.
### -------------------------------------------------------------------------

.makeTestMilo <- function(seed=42, n.cells=500, n.samples=10){
    set.seed(seed)
    ux <- matrix(rpois(20 * n.cells, 5), ncol=n.cells)
    vx <- log2(ux + 1)
    pca <- prcomp(t(vx))
    sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
                                reducedDims=SimpleList(PCA=pca$x))
    milo <- Milo(sce)
    milo <- buildGraph(milo, k=15, d=10, transposed=TRUE)
    milo <- makeNhoods(milo, k=15, d=10, prop=0.3)
    milo <- calcNhoodDistance(milo, d=10)

    per <- n.cells / n.samples
    meta <- data.frame(Sample=rep(paste0("S", seq_len(n.samples)), each=per),
                       Batch=rep(c("B1", "B2"), each=n.cells / 2))
    milo <- countCells(milo, meta.data=meta, samples="Sample")
    milo
}

.designDF <- function(n.samples=10){
    dd <- data.frame(Cond=rep(c(0, 1), n.samples / 2),
                     Batch=rep(c("B1", "B2"), each=n.samples / 2))
    rownames(dd) <- paste0("S", seq_len(n.samples))
    dd
}


test_that("GLM and GLMM return the same nhood identifiers when subsetting", {
    milo <- .makeTestMilo()
    dd <- .designDF()

    n.nh <- nrow(nhoodCounts(milo))
    skip_if(n.nh < 10, "too few neighbourhoods for a meaningful subset test")

    # drop a scattered set so the kept indices are non-consecutive - a
    # contiguous subset would not expose renumbering
    keep <- rep(TRUE, n.nh)
    keep[c(2, 5, 6)] <- FALSE

    glm.res <- testNhoods(milo, design=~Cond, design.df=dd,
                          fdr.weighting="graph-overlap", subset.nhoods=keep,
                          BPPARAM=SerialParam())
    glmm.res <- testNhoods(milo, design=~Cond + (1|Batch), design.df=dd,
                           fdr.weighting="graph-overlap", subset.nhoods=keep,
                           glmm.solver="HE-NNLS", max.iter=30,
                           # force=TRUE: the N>=60 guard is about statistical
                           # advisability, not correctness. This test is about
                           # index bookkeeping, so a small fast object is right.
                           force=TRUE,
                           BPPARAM=SerialParam())

    expect_equal(nrow(glm.res), sum(keep))
    expect_equal(nrow(glmm.res), sum(keep))

    # both branches must report the true indices of the kept neighbourhoods
    expect_equal(sort(glm.res$Nhood), which(keep))
    expect_equal(sort(glmm.res$Nhood), which(keep))
    expect_equal(sort(glm.res$Nhood), sort(glmm.res$Nhood))

    # and neither may renumber to a consecutive run, since keep has gaps
    expect_false(identical(sort(glmm.res$Nhood), seq_len(sum(keep))))
})


test_that("logCPM agrees across GLM and GLMM for the same nhood when subsetting", {
    # This is the invariant that exposed the original defect: joining on Nhood
    # must line up rows describing the same neighbourhood, and logCPM does not
    # depend on which model was fitted.
    milo <- .makeTestMilo()
    dd <- .designDF()

    n.nh <- nrow(nhoodCounts(milo))
    skip_if(n.nh < 10, "too few neighbourhoods for a meaningful subset test")
    keep <- rep(TRUE, n.nh)
    keep[c(1, 4, 7)] <- FALSE

    glm.res <- testNhoods(milo, design=~Cond, design.df=dd,
                          fdr.weighting="graph-overlap", subset.nhoods=keep,
                          BPPARAM=SerialParam())
    glmm.res <- testNhoods(milo, design=~Cond + (1|Batch), design.df=dd,
                           fdr.weighting="graph-overlap", subset.nhoods=keep,
                           glmm.solver="HE-NNLS", max.iter=30,
                           # force=TRUE: the N>=60 guard is about statistical
                           # advisability, not correctness. This test is about
                           # index bookkeeping, so a small fast object is right.
                           force=TRUE,
                           BPPARAM=SerialParam())

    cmp <- merge(glm.res, glmm.res, by="Nhood")
    expect_equal(nrow(cmp), sum(keep))

    ok <- is.finite(cmp$logCPM.x) & is.finite(cmp$logCPM.y)
    skip_if(sum(ok) < 5, "too few finite logCPM values to compare")

    # Under the defect this correlation collapsed towards zero because the rows
    # being compared described different neighbourhoods.
    #
    # Correlation rather than absolute agreement is the right invariant here:
    # the two branches compute logCPM by different routes - the GLMM branch
    # directly from the count matrix, the GLM branch via edgeR's aveLogCPM with
    # its prior count - so the values carry a small offset even when correctly
    # paired. What cannot survive a misaligned join is the correlation.
    paired.cor <- cor(cmp$logCPM.x[ok], cmp$logCPM.y[ok])
    expect_gt(paired.cor, 0.95)

    # Demonstrate the test can actually detect the regression: re-create the
    # old behaviour by renumbering the GLMM rows 1..n_kept and joining on that.
    # If this did not collapse the correlation, the check above would pass
    # whether or not the bug were present.
    renumbered <- glmm.res
    renumbered$Nhood <- rank(renumbered$Nhood)
    bad <- merge(glm.res, renumbered, by="Nhood")
    bok <- is.finite(bad$logCPM.x) & is.finite(bad$logCPM.y)
    skip_if(sum(bok) < 5, "too few rows to evaluate the misaligned join")
    expect_lt(cor(bad$logCPM.x[bok], bad$logCPM.y[bok]), paired.cor)
})


test_that("spatial FDR weights are taken from the subset, not the full set", {
    # graphSpatialFDR builds one weight per column of x.nhoods and then applies
    # w[order(pvalues)]. Passing all neighbourhoods alongside subset p-values
    # indexes that weight vector with subset ranks, so the weights land on the
    # wrong neighbourhoods. Calling it with a consistently subset x.nhoods must
    # give the same answer as testNhoods produces internally.
    set.seed(11)
    n.nh <- 60
    nh <- matrix(rbinom(400 * n.nh, 1, 0.08), nrow=400)
    colnames(nh) <- as.character(seq_len(n.nh))
    nh <- Matrix::Matrix(nh, sparse=TRUE)

    pvals <- runif(n.nh)
    keep <- rep(TRUE, n.nh)
    keep[sample(seq_len(n.nh), 12)] <- FALSE

    consistent <- graphSpatialFDR(x.nhoods=nh[, keep, drop=FALSE], graph=NULL,
                                  weighting="graph-overlap", pvalues=pvals[keep])
    mismatched <- graphSpatialFDR(x.nhoods=nh, graph=NULL,
                                  weighting="graph-overlap", pvalues=pvals[keep])

    expect_equal(length(consistent), sum(keep))
    # the two must differ, otherwise this test cannot detect a regression
    expect_false(isTRUE(all.equal(consistent, mismatched)))
    expect_true(all(consistent >= 0 & consistent <= 1))
})
