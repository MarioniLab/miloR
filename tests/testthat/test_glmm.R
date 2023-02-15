context("Testing fitGLMM function")
### Set up a mock data set using simulated data
suppressWarnings({
    library(miloR)
    library(SingleCellExperiment)
    library(scran)
    library(scater)
    library(irlba)
})

##### ------- Simulate data ------- #####
data(sim_family)
sim.df <- sim_family$DF

set.seed(42)
random.levels <- list("Fam"=paste0("Fam", unique(as.numeric(as.factor(sim.df$Fam)))))
X <- as.matrix(data.frame("Intercept"=rep(1, nrow(sim.df)), "FE2"=as.numeric(sim.df$FE2)))
Z <- as.matrix(data.frame("Fam"=as.numeric(as.factor(sim.df$Fam))))
y <- sim.df$Mean.Count
dispersion <- 0.5
mmcontrol <- glmmControl.defaults()
mmcontrol$solver <- "Fisher"

test_that("Discordant input matrices give errors", {
    # truncated y
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=Z, y=y[seq_len(nrow(X)-1)], offsets=rep(0, nrow(X)), random.levels=random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol), "Dimensions of")

    # trucated X
    set.seed(42)
    expect_error(fitGLMM(X=X[seq_len(nrow(X)-1), ], Z=Z, y=y, offsets=rep(0, nrow(X)), random.levels=random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol), "Dimensions of")

    # trucated Z
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=Z[seq_len(nrow(Z)-1), , drop=FALSE], y=y, offsets=rep(0, nrow(X)), random.levels=random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol), "Dimensions of")

    # non-square covariance matrix
    kin <- sim_family$IBD
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=Z, y=y, offsets=rep(0, nrow(X)), Kin=kin[seq_len(nrow(kin)-1), , drop=FALSE],
                         random.levels=random.levels, REML = TRUE, dispersion=dispersion, glmm.control=mmcontrol),
                 "Input covariance matrix is not square")

    # non-square covariance matrix - covariance only model
    kin <- sim_family$IBD
    g.Z <- diag(nrow(Z))
    colnames(g.Z) <- paste0("Genetic", seq_len(ncol(g.Z)))
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=g.Z, y=y, offsets=rep(0, nrow(X)), Kin=kin[seq_len(nrow(kin)-1), , drop=FALSE], geno.only=TRUE,
                         random.levels=random.levels, REML = TRUE, dispersion=dispersion, glmm.control=mmcontrol),
                 "Input covariance matrix is not square")

    # discordant covariance and Z dimensions - RE and covariance
    set.seed(42)
    expect_error(fitGLMM(X=X[seq_len(nrow(Z)-1), , drop=FALSE], Z=Z[seq_len(nrow(Z)-1), , drop=FALSE], y=y[seq_len(nrow(Z)-1)],
                         offsets=rep(0, nrow(X)), Kin=kin, random.levels=random.levels, REML = TRUE, dispersion=dispersion,
                         glmm.control=mmcontrol),
                         "Input covariance matrix and Z design matrix are discordant")

    # random levels and Z matrix are discordant
    wrong.Z <- Z
    wrong.Z[wrong.Z == 5]  <- 1 # arbitrarily re-set family IDs
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=wrong.Z,
                         y=y, offsets=rep(0, nrow(X)),
                         random.levels=random.levels, REML = TRUE, dispersion=dispersion, glmm.control=mmcontrol),
                 "Columns of Z are discordant with input random effect levels")

    # invalid column names in Z
    inv.Z <- Z
    colnames(inv.Z) <- NULL
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=inv.Z, y=y, offsets=rep(0, nrow(X)), random.levels=random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol), "Columns of Z must have valid names")

    # non unique column names in Z
    fail.random.levels <- list("RE1"=paste("RE1", unique(as.numeric(as.factor(sim.df$RE1))), sep="_"),
                               "RE2"=paste("RE1", as.numeric(unique(sim.df$RE2)), sep="_"))
    fail.X <- as.matrix(data.frame("Intercept"=rep(1, nrow(sim.df)), "FE2"=as.numeric(sim.df$FE2)))
    fail.Z <- as.matrix(data.frame("RE1"=as.numeric(as.factor(sim.df$RE1)), "RE2"=as.numeric(sim.df$RE2)))
    fail.y <- sim.df$Mean.Count

    set.seed(42)
    expect_error(fitGLMM(X=fail.X, Z=fail.Z, y=fail.y, offsets=rep(0, nrow(X)), random.levels=fail.random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol),
                 "Columns of Z are discordant with input random effect levels")

})

test_that("Infinite and NA values fail as expected", {
    inf.offsets <- rep(0, nrow(X))
    inf.offsets[sample(length(inf.offsets), size=1)] <- Inf

    set.seed(42)
    expect_error(fitGLMM(X=X, Z=Z, y=y, offsets=inf.offsets, random.levels=random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol),
                 "Infinite values in initial estimates")

    na.offsets <- rep(0, nrow(X))
    na.offsets[sample(length(na.offsets), size=1)] <- NA
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=Z, y=y, offsets=na.offsets, random.levels=random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol), "NA values in offsets")

    na.X <- X
    na.X[sample(seq_len(nrow(X)), size=1), 2] <- NA
    set.seed(42)
    expect_error(fitGLMM(X=na.X, Z=Z, y=y, offsets=rep(0, nrow(X)), random.levels=random.levels, REML = TRUE,
                           dispersion=dispersion, glmm.control=mmcontrol), "NAs values in initial estimates")

    # force infinite values with large offsets
    set.seed(42)
    expect_error(fitGLMM(X=X, Z=Z, y=y, offsets=rep(10000, nrow(X)), random.levels=random.levels, REML = TRUE,
                         dispersion=dispersion, glmm.control=mmcontrol), "Infinite values in initial estimates - reconsider model")

})



