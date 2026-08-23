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
                         dispersion=dispersion, glmm.control=mmcontrol), "Infinite values in initial estimates")

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




### -------------------------------------------------------------------------
### The genetic random effect is carried in the whitened parameterisation: the
### genetic block of Z is the Cholesky factor L of the relatedness matrix, so
### Z_g Z_g' = L L' = K, and the corresponding block of G is sigma_g I.
###
### G used to hold sigma_g K in that block instead, which - against a Z_g that
### is already L - put sigma_g L K L' into the pseudo-variance and applied K
### twice, while every partial derivative of V* with respect to sigma_g stayed
### K. The two halves described different models, so the score and information
### were not derivatives of the objective being fitted.
###
### A kinship with family block structure is used deliberately. A near-identity
### kinship leaves sigma_g confounded with the negative binomial dispersion and
### makes the comparison uninformative.
### -------------------------------------------------------------------------
.mkBlockKin <- function(n, fam=4){
    K <- diag(n)
    for(f in seq_len(n / fam)){
        i <- ((f - 1) * fam + 1):(f * fam)
        K[i, i] <- 0.5
        diag(K)[i] <- 1
    }
    K
}

.simGeneticFit <- function(seed, n=60, maxit=50){
    K <- .mkBlockKin(n)
    L <- t(chol(K))
    set.seed(seed)
    X <- cbind(1, rbinom(n, 1, 0.5))
    # varying offsets - a constant offset lies in the span of the intercept and
    # would hide offset handling entirely
    offs <- log(rnbinom(n, mu=2000, size=80))
    u.gen <- as.numeric(L %*% rnorm(n, 0, sqrt(0.25)))
    y <- as.numeric(rnbinom(n, mu=exp(offs + log(0.02) + 0.3 * X[, 2] + u.gen), size=5) + 1)
    list(n=n, K=K, L=L, X=X, y=y, offsets=offs, maxit=maxit,
         random.levels=list(Genetic=paste0("Genetic", seq_len(n))))
}

.fitKin <- function(d){
    suppressWarnings(fitGLMM(
        X=d$X, Z=matrix(1, d$n, 1, dimnames=list(NULL, "Genetic")), y=d$y,
        offsets=d$offsets, Kin=d$K, geno.only=TRUE, random.levels=d$random.levels,
        REML=TRUE, dispersion=4,
        glmm.control=list(theta.tol=1e-6, max.iter=d$maxit, solver="Fisher",
                          init.u=rep(0, d$n), init.sigma=NULL, init.beta=NULL)))
}


test_that("the genetic random effect uses the whitened parameterisation", {
    d <- .simGeneticFit(7)
    fit <- .fitKin(d)
    sg <- as.numeric(fit$Sigma[1])
    Ginv <- fit$Ginv

    expect_true(is.finite(sg) && sg > 0)
    expect_equal(dim(Ginv), c(d$n, d$n))

    # G_g = sigma_g I, so its inverse is diagonal. Under the old parameterisation
    # this block was Kinv / sigma_g, which is dense.
    expect_equal(max(abs(Ginv[upper.tri(Ginv)])), 0)
    expect_equal(sg * Ginv, diag(d$n), tolerance=1e-10)

    # and therefore the genetic block of the pseudo-variance is exactly sigma_g K
    G <- solve(Ginv)
    expect_equal(d$L %*% G %*% t(d$L), sg * d$K, tolerance=1e-8)
})


test_that("fitGLMM with a kinship matrix matches the reference genetic fitter", {
    # fitGeneticNullGlmm builds V* = W + sigma_g K directly rather than through
    # Z and G, so it is an independent implementation of the same model.
    d <- .simGeneticFit(5)
    fit <- .fitKin(d)
    b0 <- as.numeric(solve(crossprod(d$X), crossprod(d$X, log(d$y + 1) - d$offsets)))
    ref <- suppressWarnings(miloR:::fitGeneticNullGlmm(
        Z=matrix(0, d$n, 0), X=d$X, K=d$K, muvec=rep(mean(d$y), d$n), offsets=d$offsets,
        curr_beta=b0, curr_u=rep(0, d$n), curr_sigma=0.5, y=d$y, u_indices=list(),
        theta_conv=1e-6, curr_disp=4, REML=TRUE, maxit=d$maxit, disp_as_vc=FALSE))

    # both fitters use the same legacy dispersion estimator here; if they part
    # company on the dispersion the variance components are not comparable, so
    # check that first to keep a failure self-explanatory
    expect_equal(fit$Dispersion, ref$Dispersion, tolerance=1e-4)
    expect_equal(as.numeric(fit$Sigma[1]), as.numeric(ref$Sigma[1]), tolerance=1e-5)
    expect_equal(as.numeric(fit$FE), as.numeric(ref$FE), tolerance=1e-4)
})
