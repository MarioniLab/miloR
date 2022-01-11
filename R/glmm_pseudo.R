#' Perform differential abundance testing using a NB-generalised linear mixed model
#' 
#' This function will perform DA testing on all nhoods using a negative binomial generalised linear mixed model
#'
#' @param x A \code{\linkS4class{Milo}} object with a non-empty
#' \code{nhoodCounts} slot.
#' @param error.model A string vector dictating the type of error model to use for the LMM - either
#' 'normal' (default) or 'negbinom'. The former should only be used for approximately normally distributed
#' input variables, and the latter for overdispersed counts.
#'
#'
#' @importFrom MASS ginv
#' @export
runGLMM <- function(X, Z, y, init.theta=NULL, crossed=FALSE, random.levels=NULL,
                    glmm.control=list(det.tol=1e-10, cond.tol=1e-12, theta.tol=1e-6,
                                      likli.tol=1e-6, max.iter=100, lambda=1e-1, laplace.int="fe"),
                    dispersion = 0.5)){
    # model components
    # X - fixed effects model matrix
    # Z - random effects model matrix
    # A - genetic relationship matrix
    # y - observed phenotype
    
    # check the input X and Z variables for factor levels. Can we assume that if the Z's are discrete integers then
    # these are group levels, otherwise, treat them as continuous
    
    # I need to keep track of the levels/clusters in each random effect
    theta.conv <- glmm.control[["theta.tol"]] # convergence for the parameters
    loglihood.eps <- glmm.control[["likli.tol"]] # convergence for the loglihood
    
    # create full Z with expanded random effect levels
    full.Z <- initializeFullZ(Z=Z, cluster_levels=random.levels)
    
    if(is.null(init.theta)){
        # random value initiation from runif
        # curr_u <- matrix(rnorm(ncol(full.Z), mean=0, sd=1), ncol=1)
        curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        rownames(curr_u) <- colnames(full.Z)
        
        curr_beta <- ginv((t(X) %*% X)) %*% t(X) %*% log(y + 1) # OLS for the betas is usually a good starting point for NR
        #curr_beta <- matrix(c(4.7, 0.2, 0.2), ncol = 1)
        rownames(curr_beta) <- colnames(X)
    } else{
        curr_beta <- matrix(init.theta[["beta"]], ncol=1)
        rownames(curr_beta) <- colnames(X)
        
        curr_u <- matrix(init.theta[["rand"]] , ncol=1)
        rownames(curr_u) <- colnames(full.Z)
    }
    
    # compute sample variances of the us
    init.sigma <- matrix(unlist(lapply(mapUtoIndiv(full.Z, curr_u, random.levels=random.levels),
                                       FUN=function(Bj){
                                           (1/(length(Bj)-1)) * crossprod(Bj, Bj)
                                       })), ncol=1)
    
    rownames(init.sigma) <- colnames(Z)
    curr_sigma <- init.sigma
    
    init.u <- curr_u
    init.beta <- curr_beta
    init.theta <- do.call(rbind, list(init.beta, init.u))
    curr_theta <- init.theta
    
    # failure mode when the exp(Zu) estimates are infinite <- call this a convergence failure?
    inf.zu <- any(is.infinite(exp(full.Z %*% curr_u))) | any(is.na(exp(full.Z %*% curr_u)))
    
    if(inf.zu){
        stop("Infinite estimates of u")
    }
    
    mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))
    
    # use y_bar as the sample mean and s_hat as the sample variance
    y_bar <- mean(y)
    #s_hat <- var(y)
    new.r <- dispersion
    #new.r <- computeDispersion(mu.vec, s_hat) # methods of moments based estimate <- wrong!!
    
    max.hit <- glmm.control[["max.iter"]]
    
    theta_diff <- rep(Inf, nrow(curr_theta))
    abs_diff <- abs(theta_diff)
    sigma_diff <- Inf
    loglihood.diff <- Inf
    lambda <- glmm.control[["lambda"]] # lambda used in the Levenberg-Marquardt adjustment
    #init.vars <- runif(length(random.levels)) # sample variance components from ~U(0, 1)
    
    init.G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=init.sigma)
    curr_G <- init.G
    G_partials <- computeGPartials(curr_G, curr_sigma)
    
    # do a quick loglihood evaluation
    G_inv <- ginv(curr_G)
    init.loglihood  <- 0
    loglihood <- init.loglihood
    
    data_shat <- var(y)
    init.res.var <- (data_shat - colSums(curr_sigma))
    init.var.comps <- c(curr_sigma[, 1], init.res.var)/data_shat
    names(init.var.comps) <- c(rownames(curr_sigma), "residual")
    
    var.comps <- init.var.comps
    conv.list <- list()
    iters <- 1
    
    meet.conditions <- !((all(abs_diff < theta.conv)) & (loglihood.diff < loglihood.eps) & (sigma_diff < theta.conv) | iters >= max.hit)
    while(meet.conditions){
        D_inv <- computeDinv(mu.vec)
        V0 <- computeV0(mu=mu.vec, r=new.r)
        #V <- V0 + (D_inv %*% full.Z %*% curr_G %*% t(full.Z) %*% D_inv)
        V_inv <- computeVinv(V0=V0, D_inv=D_inv, Z=full.Z, G=curr_G)
        B <- computeB(y=y, r=new.r, mu=mu.vec)
        W <- computeW(mu=mu.vec, r=new.r, Z=full.Z, G=curr_G, D_inv=D_inv)
        Q <- computeQ(mu=mu.vec, r=new.r)
        
        # compute dG\dus
        Gu_partials <- computeGuPartials(curr_G=curr_G, u_hat=curr_u, cluster_levels=random.levels, sigmas=curr_sigma)
        score_beta <- betaScore(X=X, D_inv=D_inv, V_inv=V_inv, mu=mu.vec, y=y, r=new.r)
        score_u <- randScore(Z=full.Z, D_inv=D_inv, V_inv=V_inv, G_inv=G_inv, mu=mu.vec, y=y, r=new.r, u_hat=curr_u, Gu_partials=Gu_partials)
        
        full.score <- do.call(rbind, list(score_beta, score_u))
        full.hess <- jointHess(X=X, Z=full.Z, D_inv=D_inv, V_inv=V_inv, G_inv=G_inv, B=B, W=W, Q=Q, Gu_partials=Gu_partials)
        
        theta_nr.out <- singleNR(score_vec=full.score, hess_mat=full.hess, theta_hat=curr_theta)
        theta_update <- theta_nr.out$theta
        hess_theta <- theta_nr.out$hessian
        
        theta_diff <- theta_update - curr_theta # does this needs to be all negative? No, just _very_ small
        curr_theta <- theta_update
        curr_beta <- curr_theta[colnames(X), , drop=FALSE]
        curr_u <- curr_theta[colnames(full.Z), , drop=FALSE]
        
        # update sigmas _after_ thetas
        # Lets try BFGS
        indiv.u <- mapUtoIndiv(full.Z, curr_u, random.levels=random.levels)
        
        # we have estimated the parameters, and now we need to estimate the components of G (variance components)
        
        #----   INSERT NEW CODE HERE!!!   ---####
        V_star <- computeV_star(full.Z-full.Z, curr_G=curr_G, D_inv=D_inv, V0=V0)
        V_star_inv <- computeV_star_inv(V_star=V_star)
        y_star <- computey_star(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u)
        score_sigma <- sigmaScore(G_inv=G_inv, curr_u=curr_u, G_sigma_partials=G_sigma_partials, y_star=y_star)
        hessian_sigma <- sigmaHessian(G=curr_G, G_inv=G_inv, curr_u=curr_u, G_sigma_partials=G_sigma_partials)
        sigma_nr.out <- singleNR(score_vec=score_sigma, hess_mat=hessian_sigma, theta_hat=curr_sigma)
        sigma_update <- sigma_nr.out$theta

        sigma_diff <- sigma_update - curr_sigma
        curr_sigma <- sigma_update
        mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))
        
        # update G with new sigmas
        curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)
        
        if(det(curr_G) < 1e-10 | nrow(curr_sigma) == 1){
            # use a generalized inverse if numerically singular
            # this will occur if there is a single random effect too
            G_inv <- ginv(curr_G)
        } else{
            G_inv <- solve(curr_G)
        }
        
        abs_diff <- abs(theta_diff)
        
        # need to sum all of the variances
        res.var <- (data_shat - colSums(curr_sigma))
        curr_var.comps <- c(curr_sigma[, 1], res.var)/data_shat
        names(curr_var.comps) <- c(rownames(curr_sigma), "residual")
        
        # compute dispersions
        #new.r <- computeDispersion(mu.vec, s_hat) # methods of moments based estimate
        
        # do we even need to do this if we aren't using Laplace for the sigmas?
        # loglihood integrating over the random effects only
        
        if(glmm.control$laplace.int %in% c("fe")){
            fe.hess <- betaHess(X=X, D_inv=D_inv, V_inv=V_inv, B=B, W=W, Q=Q) # this is used to integrate over FEs for var comp estimation
            la.y <- indivLaplace(mu.vec, y, new.r, curr_u, curr_sigma, fe.hess, random.levels, full.Z)
            # curr.loglihood <- laplaceApprox(mu.vec, y, new.r, curr_G, G_inv, curr_u, fe.hess) # integrating over just the FEs
            curr.loglihood <- sum(la.y)
        } else if(glmm.control$laplace.int %in% c("re")){
            re.hess <- randHess(Z=full.Z, G_inv=G_inv, Gu_partials=Gu_partials, D_inv=D_inv, V_inv=V_inv, B=B, W=W, Q=Q) # this is used to integrate over FEs for var comp estimation
            la.y <- indivLaplace(mu.vec, y, new.r, curr_u, curr_sigma, re.hess, random.levels, full.Z)
            # curr.loglihood <- laplaceApprox(mu.vec, y, new.r, curr_G, G_inv, curr_u, fe.hess) # integrating over just the FEs
            curr.loglihood <- sum(la.y)
        } else if(glmm.control$laplace.int %in% c("full")){
            la.y <- indivLaplace(mu.vec, y, new.r, curr_u, curr_sigma, full.hess, random.levels, full.Z)
            # curr.loglihood <- laplaceApprox(mu.vec, y, new.r, curr_G, G_inv, curr_u, full.hess) # integrating over both the RE and FEs
            curr.loglihood <- sum(la.y)
        }
        
        loglihood.diff <- curr.loglihood - loglihood
        
        meet.conditions <- !((all(abs_diff < theta.conv)) & (abs(loglihood.diff) < loglihood.eps) &
                                 (all((sigma_diff) < theta.conv))| iters >= max.hit)
        
        loglihood <- curr.loglihood
        var.comps.diff <- curr_var.comps - var.comps
        
        # with these variance components we can evaluate the _full_ likelihood
        full.loglihood <- nbLogLikelihood(mu=mu.vec, r=new.r, y=y) + normLogLikelihood(curr_G, G_inv, curr_u)
        
        curr.u_bars <- do.call(rbind, lapply(random.levels, FUN=function(X) mean(curr_theta[X, ])))
        
        R <- V_inv - (full.Z %*% curr_G %*% t(full.Z)) #changed to V_inv?
        
        conv.list[[paste0(iters)]] <- list("Iter"=iters, "Theta"=curr_theta, "Mu"=mu.vec, "Residual"=y - mu.vec, "Loglihood"=loglihood,
                                           "Hessian"=hess_theta, "Dispersion"=new.r, "Score"=full.score, "Theta.Diff"=theta_diff, "G"=curr_G,
                                           "Rand.Mean"=curr.u_bars, "Sigmas"=curr_sigma, "LA.Y"=la.y,
                                           "R"=R,
                                           "Var.Comps"=curr_var.comps, "Var.Comp.Diff"=var.comps.diff, "Full.Loglihood"=full.loglihood)
        iters <- iters + 1
        
        var.comps <- curr_var.comps
        
        # failure mode when the exp(Zu) estimates are infinite <- call this a convergence failure?
        inf.zu <- any(is.infinite(exp(full.Z %*% curr_u))) | any(is.na(exp(full.Z %*% curr_u)))
        
        if(inf.zu){
            warning("Infinite estimates of u - estimates are diverging. Restarting estimation")
            break
        }
    }
    
    # compute SEs for final estimates
    theta_se <- computeSE(hess_theta)
    converged <- ((all(abs_diff < theta.conv)) & (abs(loglihood.diff) < loglihood.eps) & (all(abs(sigma_diff) < theta.conv)))
    
    final.list <- list("FE"=curr_theta[colnames(X), ], "RE"=curr_theta[colnames(full.Z), ], "Loglihood"=loglihood,
                       "Theta.Converged"=abs_diff < theta.conv, "LA.Y"=la.y,
                       "Loglihood.Converged"=abs(loglihood.diff) < loglihood.eps,
                       "Sigma.Converged"=sigma_diff < theta.conv,
                       "VarComp"=var.comps, "converged"=converged, "SE"=theta_se,
                       "Iters"=iters, "Dispersion"=new.r,
                       "Sigmas"=curr_sigma, "Iterations"=conv.list)
    return(final.list)
}

#done
computeV_star <- function(full.Z=full.Z, curr_G=curr_G, D_inv=D_inv, V0=V0){
    V_star = full.Z %*% curr_G %*% t(full.Z) + D_inv %*% V0 %*% D_inv
    return(V_star)
}

#done
computeV_star_inv <- function(V_star=V_star){
    V_star_inv <- if(det(V_star) < 1e-10){ 
        V_star_inv <- ginv(V_star)
    } else{
        V_star_inv <- solve(V_star)}
    return(V_star_inv)
}

#done
computey_star <- function(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u){
    y_star <- (X %*% curr_beta) + (full.Z %*% curr_u) + D_inv %*% (y - exp((X %*% curr_beta) + (full.Z %*% curr_u)))
    return(y_star)
}

#missing partial derivative
sigmaScore <- function(G_inv=G_inv, curr_u=curr_u, G_sigma_partials=G_sigma_partials){
    for i in seq_len(length(???)) {
        LHS <- -0.5*matrix.trace(V_star_inv %*% ???)
        RHS <- 0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% ??? %*% V_star_inv %*% (y_star - X %*% curr_beta)
    }
    return(LHS + RHS)
}

#missing partial derivatives
sigmaHessian <- function(G=curr_G, G_inv=G_inv, curr_u=curr_u, G_sigma_partials=G_sigma_partials) {

    hessian <- matrix(0L, ncol=length(G_sigma_partials), nrow=length(G_sigma_partials))
    for(i in seq_len(length(G_sigma_partials))) {
        for (j in seq_len(length(G_sigma_partials))){
            hessian[i, j] <- -0.5*matrix.trace(V_star_inv %*% ????? - V_star_inv %*% ??? %*% V_star_inv %*% ???) +
                0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% (????? - 2*??? %*% V_star_inv %*% ???) %*% V_star_inv %*% (y_star - X %*% curr_beta)
        }
    }
    return(hessian)
}

### write function for V*sigma
???

computeSE <- function(hessian, det.tol=1e-10){
    # compute the parameter estimate standard errors from the hessian
    hess.kappa <- 1/kappa(hessian) # small reciprocal condition number indicates need for generalised inverse
    
    if(hess.kappa <= det.tol){
        hess.se <- sqrt(diag(ginv(hessian)))
    } else{
        hess.se <- sqrt(diag(solve(hessian)))
    }
    
    names(hess.se) <- rownames(hessian)
    return(hess.se)
}

computeDispersion <- function(mu, s_hat){
    ## use my terrible methods of moments to compute a value for 'r'
    n <- length(mu)
    mu_bar <- mean(mu)
    mu_square <- sum((mu ** 2))/n
    
    r <- ((mu_square)/(s_hat + (mu_bar * (1 + mu_bar)))) - 1
    
    return(r)
}

computeGPartials <- function(curr_G, sigmas){
    # compute the partial derivatives dG/dsigma
    # return a list
    
    partial.list <- list()
    for(x in seq_len(nrow(sigmas))){
        partial.list[[x]] <- (curr_G == sigmas[x, ]) + 0 # log transforms to model scale
    }
    
    return(partial.list)
}


computeGuPartials <- function(curr_G, u_hat, cluster_levels, sigmas){
    # compute the partial derivatives dG/du
    # return a list
    n.re <- length(cluster_levels)
    partial.list <- list()
    dim.names <- unlist(cluster_levels)
    # log.sigmas <- matrix(apply(sigmas, 1, log), ncol=1)
    # rownames(log.sigmas) <- rownames(sigmas)
    # log.sigmas[is.infinite(log.sigmas), ] <- 0
    
    for(x in seq_len(n.re)){
        x.G <- matrix(0L, ncol=ncol(curr_G), nrow=nrow(curr_G),
                      dimnames=list(dim.names, dim.names))
        for(i in seq_along(cluster_levels[[x]])){
            i.x <- cluster_levels[[x]][i]
            x.G[i.x, i.x] <- 2*u_hat[cluster_levels[[x]][i], ]
        }
        
        partial.list[[x]] <- x.G
    }
    
    return(partial.list)
}

initialiseG <- function(Z, cluster_levels, sigmas){
    # construct the correct size of G given the random effects and variance components
    # names of cluster_levels and columns of Z must match
    # the independent sigmas go on the diagonal and the off-diagonal are the crossed/interactions
    # sigmas must be named
    sum.levels <- sum(unlist(lapply(cluster_levels, length)))
    G <- sparseMatrix(i=sum.levels, j=sum.levels, repr="C", x=0L)
    dimnames(G) <- list(unlist(cluster_levels), unlist(cluster_levels))
    i <- j <- 1
    # log_sigmas <- matrix(apply(sigmas, 1, log), ncol=1)
    log_sigmas <- sigmas
    log_sigmas[!is.finite(log_sigmas), ] <- 0
    rownames(log_sigmas) <- rownames(sigmas)
    
    for(x in seq_len(nrow(sigmas))){
        x.q <- length(cluster_levels[[rownames(sigmas)[x]]])
        diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- log_sigmas[x, ] # is this sufficient to transform the sigma to the model scale?
        i <- j <- i+x.q
    }
    return(as.matrix(G))
}

initializeFullZ <- function(Z, cluster_levels, stand.cols=FALSE){
    # construct the full Z with all random effect levels
    n.cols <- ncol(Z)
    col.classes <- apply(Z, 2, class)
    i.z.list <- list()
    for(i in seq_len(n.cols)){
        i.class <- col.classes[i]
        if(i.class %in% c("factor")){ # treat as factors
            i.levels <- levels(Z[, i, drop=FALSE])
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
        } else if(i.class %in% c("character")){
            i.levels <- unique(Z[, i, drop=FALSE])
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
        } else if(i.class %in% c("numeric")){ # split into unique levels if integer levels
            i.mod <- all(Z[, i, drop=FALSE] %% 1 == 0)
            if(isTRUE(i.mod)){
                i.levels <- unique(Z[, i])
                i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
            } else{
                i.z <- Z[, i, drop=FALSE] # if float then treat as continuous
            }
        } else if(i.class %in% c("integer")){
            i.levels <- unique(Z[, i])
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
        }
        colnames(i.z) <- cluster_levels[[colnames(Z)[i]]]
        
        # to standardise or not?
        if(isTRUE(stand.cols)){
            q <- ncol(i.z)
            i.ident <- diag(1L, nrow=nrow(i.z), ncol=nrow(i.z))
            i.star <- i.z - ((i.ident %*% i.z)/q)
            i.z <- i.star
        }
        
        i.z.list[[colnames(Z)[i]]] <- i.z
    }
    full.Z <- do.call(cbind, i.z.list)
    return(full.Z)
}

indivLaplace <- function(mu, y, r, curr_u, curr_sigma, hessian, random.levels, full.Z){
    ## compute the Laplace approximation to the marginal loglihood for the variance components
    ## this returns the per-individual likelihood as an (approximately) Gaussian variable
    
    # get a list of the G matrix for each random effect broadcast out to the number of observations
    G <- maptoG(curr_sigma=curr_sigma)
    Ginv <- ginv(G)
    
    indiv_u <- mapUtoIndiv(full.Z=full.Z, curr_u=curr_u, random.levels=random.levels)
    
    # I think I'll need to broadcast eveything out to the same dimensionality as G,
    # i.e. c * n by c * n, where c is the number of random effects and n is the number of observations
    
    nb.liklihood <- indivNegBinLogLikelihood(mu=mu, r=r, y=y)
    norm.liklihood <- indivNormLogLikelihood(G, Ginv, indiv_u)
    
    det.hess <- det(hessian)
    if(det.hess > 0){
        log.det.hess <- log(det.hess)
    } else{
        log.det.hess <- 0
    }
    nb.liklihood + norm.liklihood - ((1/2) * (log.det.hess/(2*pi)))
}


singleNR <- function(score_vec, hess_mat, theta_hat, lambda=1e-5, det.tol=1e-10, cond.tol=1e-15){
    # sequentially update the parameter using the Newton-Raphson algorithm
    # theta ~= theta_hat - hess^-1 * score
    # this needs to be in a direction of descent towards a minimum
    
    # need to check that hess is positive definite - if not then use a modification similar to the
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
    
    if(pos.eigens){
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
            complex.eigens <- any(is.complex(eigen(try_hess)$values))
            
            if(isTRUE(complex.eigens)){
                is.pd <- FALSE
            } else{
                is.pd <- all(eigen(try_hess)$values > 0)
            }
            
            if(isFALSE(is.pd)){
                lambda <- lambda + (2*lambda) # double lambda each time, this is very crude
            } else{
                new_hess <- try_hess
            }
        }
        
        hess.condition <- 1/kappa(new_hess)
        if(det(new_hess) < det.tol | hess.condition < cond.tol){
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


computeW <- function(mu, r, Z=full.Z, G=curr_G, D_inv=D_inv){
    # diagonal matrix containing elements of dV^-1/dmu
    w <- (((mu**2)*(1-(2*mu)))/(r)) - ((2*(mu**3))/(r**2)) - mu
    W0 <- diag(length(mu))
    diag(W0) <- w
    
    W <- W0 + (Z %*% G %*% t(Z) %*% D_inv) + (D_inv %*% Z %*% G %*% t(Z))
    return(W)
}


computeB <- function(y, r, mu){
    # diagonal matrix containing elements of d(y - db(theta)/dtheta)/dmu
    n <- length(y)
    b <- y - (n*mu) + (n*r)/(1 - (r*(mu**(-1))))
    B <- diag(n)
    diag(B) <- b
    return(B)
}


computeQ <- function(r, mu){
    # diagonal matrix containing elements of d(y - db(theta)/dtheta)/du
    n <- length(mu)
    q <- -n*(1 + (r**2)/((mu**2)*((1-r*(mu**-1))**2)))
    Q <- diag(n)
    diag(Q) <- q
    return(Q)
}


## setup functions
computeG <- function(u_hats, cluster_levels, curr_G, G_inv, sigmas, diag=FALSE){
    
    # cluster_levels should be a list with the component clusters for each random effect
    u_bars <- do.call(rbind, lapply(cluster_levels,
                                    FUN=function(UX) {
                                        mean(u_hats[UX, ])
                                    }))
    rownames(u_bars) <- names(cluster_levels)
    
    detG <- det(curr_G)
    G.partials <- computeGPartials(curr_G, sigmas=sigmas)
    G_score <- varScore(G_inv=G_inv, u_hat=u_bars, G_partials=G.partials, n.comps=nrow(sigmas))
    G_hess <- varHess(G_inv=G_inv, u_hat=u_bars, G_partials=G.partials, det.G=detG)

    # check G.hess is PD
    G.out <- singleNR(score_vec=G_score, hess_mat=G_hess, theta_hat=sigmas)
    G <- G.out$theta
    G.hess <- G.out$hessian
    
    # constrain the sigmas to max(0, G_update)
    # diag(G) <- as.numeric(sapply(diag(G), FUN=function(gx) max(0, gx)))
    
    return(list("G"=G, "G.Hess"=G.hess))
}

computeV0 <- function(mu, r){
    # compute diagonal matrix of variances
    v.vec <- ((mu**2/r)) - mu
    V0 <- matrix(0L, ncol=length(mu), nrow=length(mu))
    diag(V0) <- v.vec
    return(V0)
}

computeVinv <- function(V0, D_inv, Z, G){
    # Compute V^-1
    # need to check that V is not singular
    V <- V0 + (D_inv %*% Z %*% G %*% t(Z) %*% D_inv)
    
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


computeC <- function(G_inv, Gu_partials){
    # compute the cxc diagonal matrix of trace partial derivatives
    
    part.4.trace <- matrix(0L, ncol=ncol(G_inv), nrow=nrow(G_inv))
    for(i in seq_len(length(Gu_partials))){
        i.partial <- Gu_partials[[i]]
        for(j in seq_len(length(Gu_partials))){
            j.partial <- Gu_partials[[j]]
            part.4.trace[i, j] <- 0.5 * matrix.trace((G_inv %*% i.partial) %*% (G_inv %*% j.partial))
        }
    }
    return(part.4.trace)
}


betaScore <- function(X, D_inv, V_inv, mu, y, r){
    # score vector for fixed effects
    # kx1 vector
    n <- length(y)
    rhs <- (n*mu) - (n*r)/(1 - (r*(mu**-1)))
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


betaUHess <- function(X, Z, D_inv, V_inv, B, W, Q){
    # compute the d S(beta)/du
    part.1 <- t(X) %*% D_inv %*% V_inv %*% B %*% Z
    part.2 <- t(X) %*% D_inv %*% W %*% D_inv %*% B %*% Z
    part.3 <- t(X) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% Z
    
    hess <- part.1 + part.2 + part.3
    return(hess)
}


uBetaHess <- function(X, Z, D_inv, V_inv, B, W, Q, G_inv, C){
    # compute the d S(u)/dbeta
    # does this return a Moore-Penrose inverse?
    Z_Tinv <- ginv(Z) # note this also returns the transpose
    
    part.1 <- t(Z) %*% D_inv %*% V_inv %*% B %*% X
    part.2 <- t(Z) %*% D_inv %*% W %*% D_inv %*% B %*% X
    part.3 <- t(Z) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% X
    part.4 <- C %*% Z_Tinv %*% X
    part.5 <- G_inv %*% Z_Tinv %*% X
    
    hess <- part.1 + part.2 + part.3 - part.4 - part.5
    return(hess)
}


betaSigmaHess <- function(X, D_inv, V_inv, B, W, Q, M){
    # compute the d S(beta)\dSigma
    part.1 <- t(X) %*% D_inv %*% V_inv %*% B %*% Z %*% M
    part.2 <- t(X) %*% D_inv %*% W %*% D_inv %*% B %*% Z %*% M
    part.3 <- t(X) %*% D_inv %*% V_inv %*% Q %*% D %*% Z %*% M
    
    hess <- part.1 + part.2 + part.3
    return(hess)
}


uSigmaHess <- function(Z, D_inv, V_inv, G_inv, B, W, Q, M, U){
    #compute the dS(u)/dSigma
    part.1 <- t(Z) %*% D_inv %*% V_inv %*% B %*% Z %*% M
    part.2 <- t(Z) %*% D_inv %*% W %*% D_inv %*% B %*% Z %*% M
    part.3 <- t(Z) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% Z %*% M
    part.4 <- 0.5 * (U %*% G_inv %*% G_inv)
    
    hess <- part.1 + part.2 + part.3 - part.4
    return(hess)
}


sigmaBetaHess <- function(X, Z_inv, M_inv, G_inv, u_hat){
    # compute the dS(sigma)/dbeta
    (0.5 * (G_inv %*% G_inv) - (u_hat %*% t(u_hat) %*% G_inv %*% M_inv)) %*% t(Z_inv) %*% X
}


sigmaUHess <- function(Z_inv, M_inv, G_inv, u_hat){
    # compute the dS(sigma)/du
    0.5 * (G_inv %*% G_inv) - (u_hat %*% t(u_hat) %*% G_inv %*% M_inv)
}


jointHess <- function(X, Z, D_inv, V_inv, G_inv, B, W, Q, M, u_hat, Gu_partials){
    # construct the joint hessian for the beta's and u's
    C <- computeC(G_inv, Gu_partials)
    beta.beta <- betaHess(X, D_inv, V_inv, B, W, Q)
    rand.rand <- randHess(Z, D_inv, V_inv, G_inv, B, W, Q, Gu_partials=Gu_partials)
    
    beta.rand <- betaUHess(X, Z, D_inv, V_inv, B, W, Q)
    rand.beta <- uBetaHess(X, Z, D_inv, V_inv, B, W, Q, G_inv=G_inv, C=C)
    
    top.hess <- cbind(beta.beta, beta.rand)
    midd.hess <- cbind(rand.beta, rand.rand)
    
    full.hess <- do.call(rbind, list(top.hess, midd.hess))
    return(full.hess)
}


randScore <- function(Z, D_inv, V_inv, G_inv, mu, y, r, u_hat, Gu_partials){
    # score function for random effects
    # c X 1 vector
    n <- length(y)
    y_diff <- y - n*(mu - (r/(1 - (r/mu))))
    LHS <- t(Z) %*% D_inv %*% V_inv %*% y_diff
    RHS <- matrix(0L, ncol=1, nrow=nrow(u_hat))
    for(i in seq_len(length(Gu_partials))){
        RHS[i, ] <- 0.5 * matrix.trace(G_inv %*% Gu_partials[[i]])
    }
    
    RHS <- RHS - G_inv %*% u_hat
    # RHS <- 0.5 * (G_inv %*% u_hat) - (G_inv %*% u_hat)
    return(LHS - RHS)
}


randHess <- function(Z, D_inv, V_inv, G_inv, B, W, Q, Gu_partials){
    # compute Hessian for the random effects
    # c x c matrix
    part.1 <- t(Z) %*% D_inv %*% V_inv %*% B %*% Z
    part.2 <- t(Z) %*% D_inv %*% W %*% D_inv %*% B %*% Z
    part.3 <- t(Z) %*% D_inv %*% V_inv %*% Q %*% D_inv %*% Z
    
    part.4.trace <- matrix(0L, ncol=ncol(G_inv), nrow=nrow(G_inv))
    for(i in seq_len(length(Gu_partials))){
        i.partial <- Gu_partials[[i]]
        for(j in seq_len(length(Gu_partials))){
            j.partial <- Gu_partials[[j]]
            part.4.trace[i, j] <- 0.5 * matrix.trace((G_inv %*% i.partial) %*% (G_inv %*% j.partial))
        }
    }
    
    part.4 <- part.4.trace - G_inv
    
    hess <- part.1 + part.2 + part.3 - part.4
    
    return(hess)
}

maptoG <- function(curr_sigma){
    # map the matrix of variances to the full nxn covariance matrix
    G <- matrix(0L, ncol=nrow(curr_sigma), nrow=nrow(curr_sigma))
    rand.levels <- rownames(curr_sigma)
    
    for(j in seq_along(rand.levels)){
        G[j, j] <- curr_sigma[j, ]
    }
    
    return(G)
}

varScore <- function(sigma, random.levels, indiv.u, mu, y, r, curr_u, curr_hess, full.Z, G_inv, G_partials){
    # note some of these are unused because the fn and gradient take the same set of arguments
    # compute the score function for the sigmas
    curr_sigma <- matrix(sigma, ncol=1)
    rownames(curr_sigma) <- names(random.levels)
    
    svec <- matrix(0L, nrow=nrow(curr_sigma), ncol=1)
    
    for(i in seq_len(length(random.levels))){
        i.rand <- names(random.levels)[i]
        svec[i,] <- (-0.5 * matrix.trace(G_inv %*% G_partials[[i]])) + (0.5 * t(curr_u) %*% G_inv %*% G_partials[[i]] %*% G_inv %*% curr_u)
        # svec[i, ] <- (1/curr_sigma[i, ]) + ((4*sum(indiv.u[[i.rand]]**2))/(curr_sigma[i,]**3))
    }
    return(svec)
}


varHess <- function(curr_G, det.G, G_inv, u_hat, G_partials){
    # compute the Hessian for the sigmas
    # G_partials is a list of single-entry matrices for each dG/dsigma
    n.dims <- length(G_partials)
    var.hess <- matrix(0L, ncol=n.dims, nrow=n.dims)
    
    inv.det <- 1/(2*det.G)
    for(i in seq_len(n.dims)){
        i.partial <- G_partials[[i]]
        for(j in seq_len(n.dims)){
            j.partial <- G_partials[[j]]
            # this is almost there, but which dG\d sigma should it be that multiplies the 1/2|G| ?
            part.1 <- (0.5 * matrix.trace(G_inv %*% j.partial)) * ((matrix.trace(G_inv %*% i.partial) * matrix.trace(G_inv %*% j.partial)) - matrix.trace((G_inv %*% i.partial) %*% (G_inv %*% j.partial)))
            part.2 <- - t(u_hat) %*% G_inv %*% i.partial %*% G_inv %*% j.partial %*% G_inv %*% u_hat
            ij.hess <- part.1 + part.2
            var.hess[i, j] <- ij.hess
        }
    }
    
    return(var.hess)
}

indivNegBinLogLikelihood <- function(mu, r, y){
    ## compute the negative binomial log likelihood over our variables for each observation
    ## need to use lgamma because gamma() can't handle larger integers
    
    ## _NB_ DOUBLE CHECK THIS!!
    (y * log(1 - (r/mu))) - (r * log(1 - (r/mu))) + r * log(r) - r*log(mu) + (lgamma(y+1)-log(gamma(r)))
}


nbLogLikelihood <- function(mu, r, y){
    ## compute the negative binomial log likelihood over our variables and observations
    n <- length(y)
    y_bar <- mean(y)
    
    ## need to use lgamma because gamma() can't handle larger integers
    sum((n * y_bar * log(1 - (r/mu))) - (n * r * log(1 - (r/mu))) + r * log(r) - r*log(mu) + (lgamma(y+1)-log(gamma(r))))
}


mapUtoIndiv <- function(full.Z, curr_u, random.levels){
    # map the vector of random effects to the full nx1 vector
    rand.levels <- names(random.levels)
    indiv.u.list <- list()
    
    for(j in seq_along(rand.levels)){
        j.G <- matrix(0L, ncol=nrow(full.Z), nrow=nrow(full.Z))
        j.re <- rand.levels[j]
        j.levels <- random.levels[[j.re]]
        j.b <- full.Z[, j.levels] %*% curr_u[j.levels, ]
        indiv.u.list[[j.re]] <- j.b
    }
    
    return(indiv.u.list)
}

indivNormLogLikelihood <- function(G, Ginv, u){
    ## compute the normal log likelihood over our variables for each observation
    ## this needs to map G to the nxn matrix
    ## need to do this for each random effect
    n <- unique(unlist(lapply(u, length)))
    c <- nrow(G)
    
    big.u <- t(matrix(do.call(cbind, u), ncol=c, nrow=n))
    
    det.G <- det(G)
    if(det.G > 0){
        log.detG <- log(det.G)
    } else{
        log.detG <- 0
    }
    
    c <- length(u)
    norm.loglihoods <- c()
    
    for(j in seq_len(n)){
        norm.loglihoods <- c(norm.loglihoods,
                             sum(-((c/2) * log(2*pi)) - (0.5 * log.detG) - (0.5 * (t(big.u[, j, drop=FALSE]) %*% Ginv %*% big.u[, j, drop=FALSE]))))
    }
    return(norm.loglihoods)
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
    
    sum(-((c/2) * log(2*pi)) - (0.5 * log.detG) - (0.5 * (t(u) %*% Ginv %*% u)))
}

### utility functions

matrix.trace <- function(x){
    # check is square matrix first
    x.dims <- dim(x)
    if(x.dims[1] != x.dims[2]){
        stop("matrix is not square")
    } else{
        return(sum(diag(x)))
    }
}


#### utils
glmmControl.defaults <- function(...){
    # return the default glmm control values
    return(list(det.tol=1e-10, cond.tol=1e-12, theta.tol=1e-6,
                likli.tol=1e-6, max.iter=100, lambda=1e-1, laplace.int="fe"))
}