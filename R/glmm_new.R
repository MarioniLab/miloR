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
                                      likli.tol=1e-6, max.iter=100, lambda=1e-1, laplace.int="fe"), dispersion = 0.5){
                                        
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
        curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        curr_u <- matrix(c(-0.1315006, 0.3109791, -0.2375283), ncol = 1)
        rownames(curr_u) <- colnames(full.Z)
                            
        curr_beta <- ginv((t(X) %*% X)) %*% t(X) %*% log(y + 1) # OLS for the betas is usually a good starting point for NR
        curr_beta <- matrix(c(2, 0.04), ncol = 1)
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
    #set sigma to known value
    init.sigma <- matrix(c(0.04), ncol = 1)
        
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
    data_shat <- var(y)
    new.r <- dispersion
                        
    max.hit <- glmm.control[["max.iter"]]
                        
    theta_diff <- rep(Inf, nrow(curr_theta))
    abs_diff <- abs(theta_diff)
    sigma_diff <- Inf
    loglihood.diff <- Inf
    lambda <- glmm.control[["lambda"]] # lambda used in the Levenberg-Marquardt adjustment

    init.G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=init.sigma)
    curr_G <- init.G

    # do a quick loglihood evaluation
    G_inv <- ginv(curr_G)
    init.loglihood  <- 0
    loglihood <- init.loglihood
    
    #double check these definitions!!                    
                        init.res.var <- (data_shat - colSums(curr_sigma))
                        init.var.comps <- c(curr_sigma[, 1], init.res.var)/data_shat
                        names(init.var.comps) <- c(rownames(curr_sigma), "residual")
                        
                        var.comps <- init.var.comps
                        conv.list <- list()
                        iters <- 1
                        
   meet.conditions <- !((all(abs_diff < theta.conv)) & (loglihood.diff < loglihood.eps) & (sigma_diff < theta.conv) | iters >= max.hit)
   
   while(meet.conditions){
       
       #---- First estimate variance components with Newton Raphson procedure ---####
       D_inv <- computeDinv(mu.vec)
       D <- computeD(D_inv=D_inv)
       V <- computeV0(mu=mu.vec, r=new.r)
       W <- computeW(D=D, V=V)
       
       V_star <- computeV_star(full.Z=full.Z, curr_G=curr_G, W=W)
       V_star_inv <- computeV_star_inv(V_star=V_star)
       y_star <- computey_star(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u)
       
       score_sigma <- sigmaScore(G_inv=G_inv, curr_u=curr_u, G_sigma_partials=G_sigma_partials, y_star=y_star)
       hessian_sigma <- sigmaHessian(G=curr_G, G_inv=G_inv, curr_u=curr_u, G_sigma_partials=G_sigma_partials)
       sigma_nr.out <- singleNR(score_vec=score_sigma, hess_mat=hessian_sigma, theta_hat=curr_sigma)
       sigma_update <- sigma_nr.out$theta
       
       sigma_diff <- abs(sigma_update - curr_sigma)
       
       # update sigma, G, and G_inv
       curr_sigma <- sigma_update
       curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)
       
       if(det(curr_G) < 1e-10 | nrow(curr_sigma) == 1){
           # use a generalized inverse if numerically singular
           # this will occur if there is a single random effect too
           G_inv <- ginv(curr_G)
       } else{
           G_inv <- solve(curr_G)
       }
       
       #---- Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u---####
       theta_update <- solve_equations(X=X, W=W, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star)
                            
       theta_diff <- theta_update - curr_theta # does this needs to be all negative? No, just _very_ small
       abs_diff <- abs(theta_diff)
       
       # update B, u and mu_vec to determine new values of score and hessian matrices
       curr_theta <- theta_update
       curr_beta <- curr_theta[colnames(X), , drop=FALSE]
       curr_u <- curr_theta[colnames(full.Z), , drop=FALSE]
       mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))
       
       meet.conditions <- !((all(abs_diff < theta.conv)) & (all((sigma_diff) < theta.conv))| iters >= max.hit)
        
   }
   
}
                    
                            
       #  # need to sum all of the variances
       #  res.var <- (data_shat - colSums(curr_sigma))
       #  curr_var.comps <- c(curr_sigma[, 1], res.var)/data_shat
       #  names(curr_var.comps) <- c(rownames(curr_sigma), "residual")
       # 
       #  meet.conditions <- !((all(abs_diff < theta.conv)) & (abs(loglihood.diff) < loglihood.eps) &
       #                           (all((sigma_diff) < theta.conv))| iters >= max.hit)
       #                      
       #  var.comps.diff <- curr_var.comps - var.comps
       #  
       #  conv.list[[paste0(iters)]] <- list("Iter"=iters, "Theta"=curr_theta, "Mu"=mu.vec, "Residual"=y - mu.vec, "Loglihood"=loglihood,
       #                                                         "Hessian"=hess_theta, "Dispersion"=new.r, "Score"=full.score, "Theta.Diff"=theta_diff, "G"=curr_G,
       #                                                         "Rand.Mean"=curr.u_bars, "Sigmas"=curr_sigma, "LA.Y"=la.y,
       #                                                         "R"=R,
       #                                                         "Var.Comps"=curr_var.comps, "Var.Comp.Diff"=var.comps.diff, "Full.Loglihood"=full.loglihood)
       # iters <- iters + 1
       #                      
       # var.comps <- curr_var.comps
       #                      
       # # failure mode when the exp(Zu) estimates are infinite <- call this a convergence failure?
       # inf.zu <- any(is.infinite(exp(full.Z %*% curr_u))) | any(is.na(exp(full.Z %*% curr_u)))
       #                      
       #                      if(inf.zu){
       #                          warning("Infinite estimates of u - estimates are diverging. Restarting estimation")
       #                          break
       #                      }
       #                  }
       #                  
       #                  # compute SEs for final estimates
       #                  theta_se <- computeSE(hess_theta)
       #                  converged <- ((all(abs_diff < theta.conv)) & (abs(loglihood.diff) < loglihood.eps) & (all(abs(sigma_diff) < theta.conv)))
       #                  
       #                  final.list <- list("FE"=curr_theta[colnames(X), ], "RE"=curr_theta[colnames(full.Z), ], "Loglihood"=loglihood,
       #                                     "Theta.Converged"=abs_diff < theta.conv, "LA.Y"=la.y,
       #                                     "Loglihood.Converged"=abs(loglihood.diff) < loglihood.eps,
       #                                     "Sigma.Converged"=sigma_diff < theta.conv,
       #                                     "VarComp"=var.comps, "converged"=converged, "SE"=theta_se,
       #                                     "Iters"=iters, "Dispersion"=new.r,
       #                                     "Sigmas"=curr_sigma, "Iterations"=conv.list)
       #                  return(final.list)
       #              }

computeV_star <- function(full.Z=full.Z, curr_G=curr_G, W=W){
    W_inv <- ginv(W)
    V_star = full.Z %*% curr_G %*% t(full.Z) + W_inv
    return(V_star)
}

computeV_star_inv <- function(V_star=V_star){
    V_star_inv <- if(det(V_star) < 1e-10){ 
        V_star_inv <- ginv(V_star)
    } else{
        V_star_inv <- solve(V_star)}
    return(V_star_inv)
}

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
            i.levels <- as.factor(paste(sort(as.integer(i.levels))))
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
        } else if(i.class %in% c("character")){
            i.levels <- unique(Z[, i, drop=FALSE])
            i.levels <- as.factor(paste(sort(as.integer(i.levels))))
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
        } else if(i.class %in% c("numeric")){ # split into unique levels if integer levels
            i.mod <- all(Z[, i, drop=FALSE] %% 1 == 0)
            if(isTRUE(i.mod)){
                i.levels <- unique(Z[, i])
                i.levels <- as.factor(paste(sort(as.integer(i.levels))))
                i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
            } else{
                i.z <- Z[, i, drop=FALSE] # if float then treat as continuous
            }
        } else if(i.class %in% c("integer")){
            i.levels <- (unique(Z[, i]))
            i.levels <- as.factor(paste(sort(as.integer(i.levels))))
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


solve_equations <- function(X=X, W=W, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star){
    UpperLeft <- t(X) %*% W %*% X
    UpperRight <- t(X) %*% W %*% full.Z
    LowerLeft <- t(full.Z) %*% W %*% X
    LowerRight <- t(full.Z) %*% W %*% full.Z + G_inv
    
    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W %*% y_star), (t(full.Z) %*% W %*% y_star))
    
    theta_update <- ginv(LHS) %*% RHS
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
    # Dinv is diag(mu_i)
    Dinv <- matrix(0L, ncol=length(mu), nrow=length(mu))
    diag(Dinv) <- mu
    return(Dinv)
}


computeW <- function(D=D, V=V){
    
    W0 = (D %*% V %*% D)
    
    # now we have to take the inverse of W0
    d.det <- det(W0)
    d.kappa <- tryCatch({
        1/kappa(W0)
    },
    error=function(err){
        warning("D is computationally singular - cannot estimate condition number")
        1e-20 # arbitrarily small number
    },
    warning=function(warn) {
        1e-20
    },
    finally={
    })
    
    is.singular <- (d.det == 0) | is.infinite(d.det) | d.kappa < 1e-15
    
    if(isFALSE(is.singular)){
        W <- solve(W0)
    } else if(is.infinite(det(W0))){
        # warning("V has an infinite determinate, using a generalized inverse")
        W <- ginv(W0)
    } else {
        # warning("V is computationally singular, using a generalized inverse")
        W <- ginv(W0)
    }
    return(W)
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

computeD <- function(D_inv){
    # Compute D from D^-1
    # need to check that D is not singular
    
    d.det <- det(D_inv)
    d.kappa <- tryCatch({
        1/kappa(D)
    },
    error=function(err){
        warning("D is computationally singular - cannot estimate condition number")
        1e-20 # arbitrarily small number
    },
    warning=function(warn) {
        1e-20
    },
    finally={
    })
    
    is.singular <- (d.det == 0) | is.infinite(d.det) | d.kappa < 1e-15
    
    if(isFALSE(is.singular)){
        D <- solve(D_inv)
    } else if(is.infinite(det(D_inv))){
        # warning("V has an infinite determinate, using a generalized inverse")
        D <- ginv(D_inv)
    } else {
        # warning("V is computationally singular, using a generalized inverse")
        D <- ginv(D_inv)
    }
    
    return(D)
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