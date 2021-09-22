### Components list for the GLMM

# dispersion estimation - approximate profile likelihood a la edgeR?
# estimable/predictable functions - LR, Wald-type or score test?
# random effects definitions/covariance definition - need to be able to provide a strict covariance matrix
# define likelihood functions
# implement and test both LMM and NB-LMM


# Newton-Raphson to estimate theta (beta and u), and tau (variance components)
# need theta_hat and tau_hat, will need to iterate between the two to acheive convergence?


#' Perform differential abundance testing using a (generalised) linear mixed model
#'
#' This function will perform DA testing on all nhoods using a linear mixed model using either a
#' normal (LMM) or negative binomial (NB-LMM) error model
#' @param x A \code{\linkS4class{Milo}} object with a non-empty
#' \code{nhoodCounts} slot.
#' @param error.model A string vector dictating the type of error model to use for the LMM - either
#' 'normal' (default) or 'negbinom'. The former should only be used for approximately normally distributed
#' input variables, and the latter for overdispersed counts.
#'
#'
#' @importFrom MASS ginv
#' @export
runLMM <- function(X, Z, y, error.model=c('normal', 'negbinom')){
    # model components
    # X - fixed effects model matrix
    # Z - random effects model matrix
    # A - genetic relationship matrix
    # y - observed phenotype
    out.list <- list()

    # initial parameter values
    init.params <- list("alpha"=exp(runif(1)), "betas"=exp(runif(ncol(X))), "us"=runif(ncol(Z)), "sigmas"=exp(runif(ncol(Z))))
    init.params$zi <- (X %*% init.params$betas) + (Z %*%  init.params$us)

    total.var <- var(y)
    exp.y <- init.params$zi
    var.y <- exp.y + exp.y/init.params$alpha

    init.taus <- c(init.params$alpha, init.params$sigmas)
    init.theta <- c(init.params$betas, init.params$us)

    # eta is the log-linked predictor
    eta <- exp((X %*% init.params$betas) + (Z %*% init.params$us))

    # compute the derivatives for NR
    v <- computeV(y, init.params$alpha, init.params$zi)

    # fixed effect betas
    beta.deriv <- logLikBeta(X, v)
    # invert G
    G <- computeG(init.params$sigmas)
    Ginv <- solve(G)

    # random effects u
    u.deriv <- logLikU(Z, v, Ginv, init.params$us)

    # compute R (co-)variance matrix
    R <- computeR(init.params$alpha, y, init.params$zi)

    # compute the hessian
    out.list[["hessian"]] <- computeHessian(X, R, Z, Ginv)

    # return results
    out.list[["estimates"]] <- do.call(rbind, list(beta.deriv, u.deriv))

    y_bar <- computeYbar(X, init.params$betas, Z, init.params$us, R, v)

    # test solving for new thetas
    new.thetas <- solveMME(X, Z, out.list[["hessian"]], y_bar, R)

    # update the zis
    new.zis <- (X %*% new.thetas[1:ncol(X)]) + (Z %*%  new.thetas[(ncol(X)+1):length(new.thetas)])

    tau.hessian <- as.matrix(out.list[["hessian"]][((ncol(X)+1):length(new.thetas)), (ncol(X)+1):length(new.thetas)])

    tau.likelihood <- laplaceApprox(alpha=init.params$alpha, y=y, zi=new.zis, G, Ginv,
                                    u=new.thetas[(ncol(X)+1):length(new.thetas)], tau.hessian=tau.hessian)

}


singleNR <- function(score_vec, hess_mat, theta_hat, lambda=1e-5, det.tol=1e-10, cond.tol=1e-15){
    # sequentially update the parameter using the Newton-Raphson algorithm
    # theta ~= theta_hat - hess^-1 * score

    # need to check that hess is positive definite - if not then use a modification simialar to the
    # Levenberg-Marquardt method to make it slightly pd
    # a pd matrix has all positive eigen values
    # if eigen values have opposite signs then we have a saddle point

    hess.eigen <- eigen(hess_mat)
    hess.condition <- 1/kappa(hess_mat)

    # check if complex eigen values
    complex.eigens <- any(is.complex(hess.eigen$values))
    if(isTRUE(complex.eigens)){
        pos.eigens <- FALSE
    } else{
        pos.eigens <- all(hess.eigen$values > 0)
    }

    if(pos.eigens | complex.eigens | hess.condition < det.tol){
        if(det(hess_mat) < det.tol | hess.condition < cond.tol){
            theta_new <- theta_hat - ginv(hess_mat) %*% score_vec
        } else{
            theta_new <- theta_hat - solve(hess_mat) %*% score_vec
        }
    } else{
        # warning("Hessian is not positive definite - attempting modified Levenberg-Marquardt adjustment")
        # use the Levenberg-Marquardt method
        # replace hess with hess + lambda * I, such that lambda is chosen large enough to make hessian pd
        is.pd <- FALSE
        while(isFALSE(is.pd)){
            try_hess <- hess_mat + (diag(ncol(hess_mat)) * lambda)
            is.pd <- all(eigen(try_hess)$values > 0)

            if(isFALSE(is.pd)){
                lambda <- lambda + (2*lambda) # double lambda each time, this is very crude
            } else{
                new_hess <- try_hess
            }
        }

        if(det(new_hess) < det.tol){
            theta_new <- theta_hat - ginv(new_hess) %*% score_vec
        } else{
            theta_new <- theta_hat - solve(new_hess) %*% score_vec
        }
        hess_mat <- new_hess
    }

    return(list("theta"=theta_new, "hessian"=hess_mat))
}

computeDinv <- function(mu){
    # D is a matrix of 1/mu, so Dinv is diag(mu_i)
    Dinv <- matrix(0L, ncol=length(mu), nrow=length(mu))
    diag(Dinv) <- mu
    return(Dinv)
}

computeW <- function(mu, r){
    # diagonal matrix containing elements of dV^-1/dmu
    w <- (((mu**2)*(1-(2*mu)))/(r)) - ((2*(mu**3))/(r**2)) - mu
    W <- diag(length(mu))
    diag(W) <- w
    return(W)
}

computeB <- function(y, r, mu){
    # diagonal matrix containing elements of d(y - db(theta)/dtheta)/dmu
    n <- length(y)
    b <- y - (n*mu) - (n*r)/(1 - (r*(mu**-1)))
    B <- diag(n)
    diag(B) <- b
    return(B)
}

computeQ <- function(r, mu){
    # diagonal matrix containing elements of d(y - db(theta)/dtheta)/du
    n <- length(mu)
    q <- -n*(1 + (r**2)/(mu**2*(1-r*(mu**-1))))
    Q <- diag(n)
    diag(Q) <- q
    return(Q)
}

## setup functions
computeG <- function(u_hats, cluster_levels, curr_G, diag=FALSE){
    # compute G from the estimates of u=u_hat
    if(length(u_hats) > 1){
        c <- length(u_hats) - 1
    } else{
        c <- 1
    }

    # cluster_levels should be a list with the component clusters for each random effect
    # I also need to consider the covariances between random effect cluster levels
    u_bars <- do.call(rbind, lapply(cluster_levels,
                                    FUN=function(UX) {
                                        mean(u_hats[UX, ])
                                        }))
    rownames(u_bars) <- names(cluster_levels)

    G_score <- varScore(G_inv=curr_G, u_hat=u_bars)
    G_hess <- varHess(G_inv=curr_G, u_hat=u_bars)

    G <- singleNR(score_vec=G_score, hess_mat=G_hess, theta_hat=curr_G)$theta

    return(G)
}

computeVinv <- function(mu, y, r){
    # compute diagonal matrix of inverse variances
    v.vec <- (mu**2 * (1/r)) - mu
    V <- matrix(0L, ncol=length(mu), nrow=length(mu))
    diag(V) <- v.vec

    # need to check that V is not singular
    v.det <- det(V)
    v.kappa <- tryCatch({
        1/kappa(V)
        },
        error=function(err){
            warning("V is computationally singular - cannot estimate condition number")
            1e-20 # arbitrarily small number
            },
        warning=function(warn) {
            1e-20
        },
        finally={
        })

    is.singular <- (v.det == 0) | is.infinite(v.det) | v.kappa < 1e-15

    if(isFALSE(is.singular)){
        Vinv <- solve(V)
    } else if(is.infinite(det(V))){
        # warning("V has an infinite determinate, using a generalized inverse")
        Vinv <- ginv(V)
    } else {
        # warning("V is computationally singular, using a generalized inverse")
        Vinv <- ginv(V)
    }

    return(Vinv)
}

betaScore <- function(X, D_inv, V_inv, mu, y, r){
    # score vector for fixed effects
    # kx1 vector
    n <- length(y)
    rhs <- y - (n*mu) - (n*r)/(1 - (r*(mu**-1)))
    return(t(X) %*% D_inv %*% V_inv %*% (y - rhs))
}

betaHess <- function(X, D_inv, V_inv, B, W, Q){
    # compute hessian for beta's
    # k x k matrix
    part.1 <- t(X) %*% D_inv %*% V_inv %*% B %*% X
    part.2 <- t(X) %*% D_inv %*% W %*% D_inv %*% B %*% X
    part.3 <- t(X) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% X
    hess <- part.1 + part.2 + part.3

    return(hess)
}

betaUHess <- function(X, Z, D_inv, V_inv, mu, y, r){
    # compute the d S(beta)/du
    part.1 <- t(X) %*% D_inv %*% V_inv %*% B %*% Z
    part.2 <- t(X) %*% D_inv %*% W %*% D_inv %*% B %*% Z
    part.3 <- t(X) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% Z

    hess <- part.1 + part.2 + part.3
    return(hess)
}


uBetaHess <- function(X, Z, D_inv, V_inv, mu, y, r){
    # compute the d S(u)/dbeta
    part.1 <- t(Z) %*% D_inv %*% V_inv %*% B %*% X
    part.2 <- t(Z) %*% D_inv %*% W %*% D_inv %*% B %*% X
    part.3 <- t(Z) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% X

    hess <- part.1 + part.2 + part.3
    return(hess)
}

jointHess <- function(X, Z, D_inv, V_inv, G_inv, B, W, Q){
    # construct the full hessian
    beta.beta <- betaHess(X, D_inv, V_inv, B, W, Q)
    rand.rand <- randHess(Z, D_inv, V_inv, G_inv, B, W, Q)

    beta.rand <- betaUHess(X, Z, D_inv, V_inv, B, W, Q)
    rand.beta <- uBetaHess(X, Z, D_inv, V_inv, B, W, Q)

    top.hess <- cbind(beta.beta, beta.rand)
    bottom.hess <- cbind(rand.beta, rand.rand)

    full.hess <- rbind(top.hess, bottom.hess)
    return(full.hess)
}

randScore <- function(Z, D_inv, V_inv, G_inv, mu, y, r, u_hat){
    # score function for random effects
    # c X 1 vector
    n <- length(y)
    y_diff <- y - (n*mu) - (n*r)/(1 - (r*(mu**-1)))
    LHS <- t(Z) %*% D_inv %*% V_inv %*% y_diff
    RHS <- 0.5 * (G_inv %*% u_hat) - (G_inv %*% u_hat)
    return(LHS - RHS)
}

randHess <- function(Z, D_inv, V_inv, G_inv, B, W, Q){
    # compute Hessian for the random effects
    # c x c matrix
    part.1 <- t(Z) %*% D_inv %*% V_inv %*% B %*% Z
    part.2 <- t(Z) %*% D_inv %*% W %*% D_inv %*% B %*% Z
    part.3 <- t(Z) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% Z
    part.4 <- (0.5 * G_inv) - G_inv

    hess <- part.1 + part.2 + part.3 - part.4

    return(hess)
}


varScore <- function(G_inv, u_hat){
    # compute the score function for the sigmas
    (-0.5 * G_inv) + ((u_hat %*% t(u_hat)) %*% G_inv %*% G_inv)
}

varHess <- function(G_inv, u_hat){
    # compute the Hessian for the sigmas
    (0.5 * (G_inv %*% G_inv)) - (2*(u_hat %*% t(u_hat)) %*% G_inv)
}


laplaceApprox <- function(mu, y, r, G, G_inv, curr_u, hessian){
    ## compute the Laplace approximation to the full extended likelihood
    nb.liklihood <- nbLogLikelihood(mu=mu, r=r, y=y)
    norm.liklihood <- normLogLikelihood(G, G_inv, curr_u)

    det.hess <- det(hessian)
    nb.liklihood + norm.liklihood - ((2*pi) * det.hess)
}


nbLogLikelihood <- function(mu, r, y){
    ## compute the negative binomial log likelihood over our variables and observations
    n <- length(y)
    sum((n * y * log(1 - (r/mu))) - (n * r * log(1 - (r/mu))) + r * log(r) + r*log(mu) + (log(gamma(y+1))/log(gamma(r))))
}


normLogLikelihood <- function(G, Ginv, u){
    ## compute the normal log likelihood over our variables and observations
    det.G <- det(G)
    if(det.G > 0){
        log.detG <- log(det.G)
    } else{
        log.detG <- 0
    }
    u.mat <- as.matrix(u, ncol=1)
    c <- length(u)

    sum(-((c/2) * log(2*pi)) -0.5 * log.detG - (0.5 * (t(u) %*% Ginv %*% u)))
}















######## old  code



newtonRaphson <- function(X, init.beta, Z, init.u, y, y_bar, alphas, sigmas, theta_new, max.iter=10, tol=1e-6){
    # use NR to iteratively estimate beta's and u's
    theta <- c(init.beta, init.u)
    theta.diff <- abs(theta_new - theta)

    var.y <- var(y)
    resid.var <- var.y - sum(sigmas) # variance NOT explained by model parameters

    betas <- theta[1:length(init.beta)]
    us <- theta[(length(init.beta)+1):length(theta)]
    zis <- (X %*% betas) + (Z %*%  us)

    # compute an initial estimate of G
    alphas <- alphas # how to update the alphas? this is from the Laplace/derivative-free method?
    sigma.e <- trigamma(alphas) # as per Templeman & Gianola (1999) - eq [8]
    sigmas <- sigmas # are these the diagonal elements of G?

    iter.t <- 1
    iter.list <- list()
    out.list <- list()
    while(any(theta.diff > tol)){
        # compute the derivatives for NR
        v <- computeV(y, alphas, zis)

        # fixed effect betas
        beta.deriv <- logLikBeta(X, v)

        # invert G
        G <- computeG(sigmas)
        Ginv <- solve(G)

        # random effects u
        u.deriv <- logLikU(Z, v, Ginv, us)

        first.derivs <- do.call(rbind, list(beta.deriv, u.deriv))

        # compute R (co-)variance matrix
        R <- computeR(alphas, y, zis)

        # compute the hessian
        hessian <- computeHessian(X, R, Z, Ginv)
        y_bar <- computeYbar(X, betas, Z, us, R, v)

        # test solving for new thetas
        theta_new <- solveMME(X, Z, hessian, y_bar, R)
        theta.diff <- abs(theta_new - theta)

        betas <- theta_new[1:length(init.beta), , drop=FALSE]
        us <- theta_new[(length(init.beta)+1):length(theta), , drop=FALSE]
        zis <- (X %*% betas) + (Z %*%  us)
        sigmas <- diag(G)

        # need to update the alphas
        message("Theta difference: ", paste(theta.diff, collapse=", "), " on iteration ", iter.t)
        iter.t <- iter.t + 1

        if(iter.t > max.iter){
            warning("Maximum iterations met: ", max.iter, " - model may not have converged. Theta difference: ", paste(theta.diff, collapse=", "))
            theta <- theta_new
            break
        }

        theta <- theta_new

        # compute new marginal log likelihood
        tau.hessian <- as.matrix(hessian[((ncol(X)+1):length(theta_new)), (ncol(X)+1):length(theta_new)])
        mml <- laplaceApprox(alpha=alphas, y=y, zi=zis, G, Ginv, u=us, tau.hessian=tau.hessian)
        message("Log-likelihood: ", mml)
        iter.list[[iter.t-1]] <- list("theta"=theta, "MML"=mml)
    }
    out.list$Iterations <- iter.list
    out.list$Final.Theta <- theta
    out.list$Hessian <- hessian

    out.list
}


computeR <- function(alpha, y, zi){
    # remaining variance components - R is a _diagonal_ matrix
    I.n <- diag(length(y))
    diag(I.n) <- (((alpha + y) * alpha * zi)/((zi + alpha) **2))
    I.n
}


# ### functions to compute derivatives
# computeV <- function(y, alpha, zi){
#     # the components of v are the product of the partial derivs w.r.t. zis and the zis w.r.t. eta
#     y - (((alpha + y)*zi)/(zi + alpha))
# }

logLikBeta <- function(X, v){
    # partial derivatives for fixed effect estimation with NR
    t(X) %*% v
}

logLikU <- function(Z, v, Ginv, u){
    # partial derivatives for random effect estimation with NR
    (t(Z) %*% v) - (Ginv %*% u)
}

### functions to compute likelihoods

computeHessian <- function(X, R, Z, Ginv, eigen.tol=1e-50){
    ## compute the Hessian matrix as the partial second derivatives
    h1 <- t(X) %*% R %*% X
    h2 <- t(X) %*% R %*% Z
    h3 <- t(Z) %*% R %*% X
    h4 <- t(Z) %*% R %*% Z + Ginv

    hess <- cbind(h1, h2)
    hess <- rbind(hess, cbind(h3, h4))

    # check that the hessian is positive-definite, i.e. are all of the eigen values > 0 within some tolerance?
    hess.eigen <- eigen(hess)$values
    eigen.diffs <- abs(0 - abs(hess.eigen))
    if(any(eigen.diffs < eigen.tol)){
        warning("Hessian may not be positive (semi)-definite. Eigen values close to zero: ", paste(hess.eigen[which(eigen.diffs < eigen.tol)], collapse=", "))
    }

    return(hess)
}

computeResiduals <- function(X, beta, Z, u, y){
    mu_hat <- (X %*% beta) + (Z %*% u)
    return(y - mu_hat)
}

computeYbar <- function(X, beta, Z, u, R, v){
    #  comoyte a vector of working/estimated Y values for NR solutions
    Rinv <- solve(R)

    (X %*% beta) + (Z %*% u) + (Rinv %*% v)
}

computeAlpha <- function(zi, y){
    y_mean <- mean(y)
    denom <- sum((y - y_mean) **2)/(length(y) - 1)
    ((zi ** 2)/denom) - zi
}


solveMME <- function(X, Z, hessian, y_bar, R){
    ## solve the mixed-model equations at theta using Cholesky factorization
    # hessian %*% theta_est = [t(X)Ry_bar t(Z)Ry_bar]^T
    # Ax = b
    # solve Ly=b then L*x=y

    b.x <- t(X) %*% R %*% y_bar
    b.z <- t(Z) %*% R %*% y_bar
    b <- rbind(b.x, b.z)

    # solving using cholesky
    hess.chol <- chol(hessian)
    z <- backsolve(hess.chol, hessian %*% b, transpose=TRUE) # solve the backwards equation
    new_theta <- backsolve(hess.chol, z)
    new_theta
}


# laplaceApprox <- function(alpha, y, zi, G, Ginv, u, tau.hessian){
#     ## compute the Laplace approximation to the full extended likelihood
#     nb.liklihood <- nbLogLikelihood(alpha, y, zi)
#     norm.liklihood <- normLogLikelihood(G, Ginv, u)
#
#     det.hess <- det(tau.hessian)
#     as.numeric(nb.liklihood + norm.liklihood - (0.5 * det.hess))
# }
#
#
# nbLogLikelihood <- function(alpha, y, zi){
#     ## compute the negative binomial log likelihood over our variables and observations
#     n <- length(y)
#     n  * (alpha * log(alpha) - log(gamma(alpha)) + sum(log(gamma(y + alpha)) - (alpha * log(zi + alpha)) + (y * log(zi)) - log(zi + alpha)))
# }
#
#
# normLogLikelihood <- function(G, Ginv, u){
#     ## compute the normal log likelihood over our variables and observations
#     det.G <- det(G)
#     u.mat <- as.matrix(u)
#
#     -0.5 * log(det.G) - (0.5 * (t(u) %*% Ginv %*% u))
# }
#
#
#
#
