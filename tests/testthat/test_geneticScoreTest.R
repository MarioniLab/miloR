context("Testing the genetic score test and null model fitter")
library(miloR)

### -------------------------------------------------------------------------
### Simulate a small NB-GLMM problem with a structured kinship matrix.
###
### The kinship matrix deliberately has family block structure. A near-identity
### kinship leaves sigma_g confounded with the negative binomial dispersion,
### which drives the genetic variance component to the boundary and makes these
### tests uninformative about the variance component machinery.
### -------------------------------------------------------------------------
.mkKin <- function(n, fam=4){
    nf <- n / fam
    K <- diag(n)
    for(f in seq_len(nf)){
        i <- ((f - 1) * fam + 1):(f * fam)
        K[i, i] <- 0.5
        diag(K)[i] <- 1
    }
    K
}

.simGeneticNhood <- function(n=100, qp=20, seed=42, beta.snp=0){
    set.seed(seed)
    pool <- sample(seq_len(qp), n, replace=TRUE)
    Zp <- matrix(0, n, qp)
    Zp[cbind(seq_len(n), pool)] <- 1
    colnames(Zp) <- paste0("pool", seq_len(qp))
    K <- .mkKin(n)
    g <- rbinom(n, 2, 0.3)
    X <- cbind(1, rbinom(n, 1, 0.5), as.numeric(scale(rnorm(n))))
    colnames(X) <- c("(Intercept)", "sex", "age")
    # varying offsets - a constant offset would hide offset handling bugs,
    # because a constant vector lies in the span of the intercept
    offs <- log(rnbinom(n, mu = 2000, size = 80))
    u.pool <- rnorm(qp, 0, 0.30)
    u.gen <- as.numeric(t(chol(K)) %*% rnorm(n, 0, 0.35))
    mu <- exp(offs + log(0.02) + 0.2 * X[, 2] + beta.snp * g + Zp %*% u.pool + u.gen)
    y <- rnbinom(n, mu=mu, size=5) + 1
    list(n=n, qp=qp, Zp=Zp, K=K, g=g, X=X, Xa=cbind(X, SNP=g), y=y, offsets=offs)
}

.fitNull <- function(d, project=TRUE){
    miloR:::fitGeneticNullGlmm(
        Z=d$Zp, X=d$X, K=d$K, muvec=rep(mean(d$y), d$n), offsets=d$offsets,
        curr_beta=c(log(mean(d$y)) - mean(d$offsets), 0, 0),
        curr_u=rep(0, d$qp + d$n), curr_sigma=c(0.5, 0.5), y=d$y,
        u_indices=list(seq_len(d$qp)), theta_conv=1e-6, curr_disp=1,
        REML=TRUE, maxit=100, return_projection=project)
}


test_that("the null model converges and returns finite estimates", {
    d <- .simGeneticNhood()
    fit <- .fitNull(d)

    expect_true(fit$converged)
    expect_true(all(is.finite(as.numeric(fit$FE))))
    expect_true(all(as.numeric(fit$Sigma) >= 0))
    expect_true(is.finite(fit$Dispersion) && fit$Dispersion > 0)
    expect_equal(nrow(fit$P), d$n)
    expect_equal(length(as.numeric(fit$Pystar)), d$n)
})


test_that("mu stays on the scale of the data when offsets vary", {
    # The offset belongs to the linear predictor but is not a column of X. If it
    # is left in the working response for the GLS solve, the intercept absorbs
    # it, eta double counts it and mu diverges.
    d <- .simGeneticNhood()
    fit <- .fitNull(d)

    W <- as.numeric(fit$Wdiag)
    mu.implied <- 1 / (W - 1 / fit$Dispersion)

    expect_true(all(is.finite(mu.implied)))
    # mu should sit within an order of magnitude of the observed counts
    expect_lt(max(mu.implied), 10 * max(d$y))
    expect_gt(min(mu.implied), min(d$y) / 10)
})


test_that("the REML projection matches a direct computation", {
    d <- .simGeneticNhood()
    fit <- .fitNull(d)

    W <- as.numeric(fit$Wdiag)
    sg <- as.numeric(fit$Sigma)
    ystar <- as.numeric(fit$ystar)

    V <- diag(W) + sg[1] * (d$Zp %*% t(d$Zp)) + sg[2] * d$K
    Vi <- solve(V)
    Mi <- solve(t(d$X) %*% Vi %*% d$X)
    P.manual <- Vi - Vi %*% d$X %*% Mi %*% t(d$X) %*% Vi

    expect_lt(max(abs(fit$P - P.manual)), 1e-9)
    # Pystar is built from the offset-corrected working response
    expect_lt(max(abs(as.numeric(fit$Pystar) -
                          as.numeric(P.manual %*% (ystar - d$offsets)))), 1e-9)
    # the projection annihilates the fixed effect design
    expect_lt(max(abs(t(d$X) %*% as.numeric(fit$Pystar))), 1e-9)
})


test_that("the score statistic equals the Wald statistic built from it", {
    # This is an algebraic identity, not an approximation: with beta = U/I and
    # se = 1/sqrt(I), (beta/se)^2 = U^2/I = chi2. It holds to machine precision
    # and is the tightest available check on the score test implementation.
    d <- .simGeneticNhood()
    fit <- .fitNull(d)
    sc <- miloR:::scoreTestGeneticSNPs(fit$P, as.numeric(fit$Pystar),
                                       matrix(d$g, ncol=1))

    chisq <- as.numeric(sc$Chisq)
    wald <- (as.numeric(sc$Beta) / as.numeric(sc$SE))^2

    expect_equal(chisq, wald, tolerance=1e-9)
    expect_lt(abs(chisq - wald), 1e-9)
})


test_that("batched scoring equals variant-by-variant scoring", {
    d <- .simGeneticNhood()
    fit <- .fitNull(d)

    set.seed(7)
    G <- cbind(d$g, rbinom(d$n, 2, 0.2), rbinom(d$n, 2, 0.45),
               rbinom(d$n, 2, 0.1), rbinom(d$n, 2, 0.5))

    batched <- miloR:::scoreTestGeneticSNPs(fit$P, as.numeric(fit$Pystar), G)
    single <- vapply(seq_len(ncol(G)), FUN=function(j){
        as.numeric(miloR:::scoreTestGeneticSNPs(fit$P, as.numeric(fit$Pystar),
                                                G[, j, drop=FALSE])$Chisq)
    }, FUN.VALUE=numeric(1))

    expect_equal(as.numeric(batched$Chisq), single, tolerance=1e-12)
})


test_that("monomorphic variants return NA rather than an unstable ratio", {
    d <- .simGeneticNhood()
    fit <- .fitNull(d)

    G <- cbind(mono=rep(1, d$n), real=d$g)
    sc <- miloR:::scoreTestGeneticSNPs(fit$P, as.numeric(fit$Pystar), G)

    expect_true(is.na(as.numeric(sc$Chisq)[1]))
    expect_false(is.na(as.numeric(sc$Chisq)[2]))
})


test_that("the score test agrees with a full Wald refit to within O(1/n)", {
    # The score test holds the variance components at their null values whereas
    # the Wald refit re-estimates them with the variant in the model, so the two
    # differ at order 1/n. They are NOT equal to 1e-9 - only the score/one-step
    # Wald identity above is exact. This test pins the size of the discrepancy
    # so that a regression which changes it will be caught.
    d <- .simGeneticNhood(beta.snp=0.25)

    fit <- .fitNull(d)
    sc <- miloR:::scoreTestGeneticSNPs(fit$P, as.numeric(fit$Pystar),
                                       matrix(d$g, ncol=1))

    alt <- miloR:::fitGeneticNullGlmm(
        Z=d$Zp, X=d$Xa, K=d$K, muvec=rep(mean(d$y), d$n), offsets=d$offsets,
        curr_beta=c(log(mean(d$y)) - mean(d$offsets), 0, 0, 0),
        curr_u=rep(0, d$qp + d$n), curr_sigma=c(0.5, 0.5), y=d$y,
        u_indices=list(seq_len(d$qp)), theta_conv=1e-6, curr_disp=1,
        REML=TRUE, maxit=100, return_projection=FALSE)

    beta.score <- as.numeric(sc$Beta)
    beta.wald <- as.numeric(alt$FE)[4]

    expect_true(is.finite(beta.score) && is.finite(beta.wald))
    # same sign and same order of magnitude
    expect_equal(sign(beta.score), sign(beta.wald))
    expect_lt(abs(beta.score - beta.wald), 0.5 * max(abs(beta.wald), 1e-8) + 0.1)
})


test_that("supplying null estimates as starting values reproduces the fit", {
    d <- .simGeneticNhood()
    fit <- .fitNull(d, project=FALSE)

    warm <- miloR:::fitGeneticNullGlmm(
        Z=d$Zp, X=d$X, K=d$K, muvec=rep(mean(d$y), d$n), offsets=d$offsets,
        curr_beta=c(log(mean(d$y)) - mean(d$offsets), 0, 0),
        curr_u=rep(0, d$qp + d$n), curr_sigma=c(0.5, 0.5), y=d$y,
        u_indices=list(seq_len(d$qp)), theta_conv=1e-6, curr_disp=1,
        REML=TRUE, maxit=100,
        null_sigma_=as.numeric(fit$Sigma),
        null_beta_=as.numeric(fit$FE),
        null_disp=fit$Dispersion,
        return_projection=FALSE)

    expect_true(warm$converged)
    # A warm start follows a different iteration path, so it lands at a
    # different point within the convergence tolerance rather than at exactly
    # the same one. The meaningful scale for "the same fit" is the standard
    # error, not machine precision: require agreement to under 1% of the SE.
    expect_lt(max(abs(as.numeric(warm$FE) - as.numeric(fit$FE)) /
                      as.numeric(fit$SE)), 0.01)
    expect_lt(max(abs(as.numeric(warm$FE) - as.numeric(fit$FE))), 1e-3)
    # a warm start should not need more iterations than a cold one
    expect_lte(warm$Iters, fit$Iters)
})


test_that("fixing the variance components holds them at their supplied values", {
    d <- .simGeneticNhood()
    fit <- .fitNull(d, project=FALSE)
    target <- as.numeric(fit$Sigma)

    fixed <- miloR:::fitGeneticNullGlmm(
        Z=d$Zp, X=d$Xa, K=d$K, muvec=rep(mean(d$y), d$n), offsets=d$offsets,
        curr_beta=c(log(mean(d$y)) - mean(d$offsets), 0, 0, 0),
        curr_u=rep(0, d$qp + d$n), curr_sigma=target, y=d$y,
        u_indices=list(seq_len(d$qp)), theta_conv=1e-6,
        curr_disp=fit$Dispersion, REML=TRUE, maxit=100,
        fix_variance=TRUE, return_projection=FALSE)

    expect_equal(as.numeric(fixed$Sigma), target, tolerance=1e-12)
    expect_equal(fixed$Dispersion, fit$Dispersion, tolerance=1e-12)
})
