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
                    
    theta.conv <- glmm.control[["theta.tol"]] # convergence for the parameters
    max.hit <- glmm.control[["max.iter"]]

    # create full Z with expanded random effect levels
    full.Z <- initializeFullZ(Z=Z, cluster_levels=random.levels)
                        
    if(is.null(init.theta)){
        # random value initiation from runif
        curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        #curr_u <- matrix(c(0.09918169, -0.16222715, 0.06304546), ncol = 1)
        rownames(curr_u) <- colnames(full.Z)
        
        # OLS for the betas is usually a good starting point for NR                    
        curr_beta <- ginv((t(X) %*% X)) %*% t(X) %*% log(y + 1) 
        #curr_beta <- matrix(c(2, 0.25), ncol = 1)
        rownames(curr_beta) <- colnames(X)
        
    } else{
        curr_beta <- matrix(init.theta[["beta"]], ncol=1)
        rownames(curr_beta) <- colnames(X)
                            
        curr_u <- matrix(init.theta[["rand"]] , ncol=1)
        rownames(curr_u) <- colnames(full.Z)
    }
                        
    # compute sample variances of the us
    curr_sigma <- matrix(unlist(lapply(mapUtoIndiv(full.Z, curr_u, random.levels=random.levels),
                                       FUN=function(Bj){
                                           (1/(length(Bj)-1)) * crossprod(Bj, Bj)
                                           })), ncol=1)
    #curr_sigma <- matrix(c(260), ncol = 1)
    rownames(curr_sigma) <- colnames(Z)
    
    #create a single variable for the thetas
    curr_theta <- do.call(rbind, list(curr_beta, curr_u))

    # failure mode when the exp(Zu) estimates are infinite <- call this a convergence failure?
    inf.zu <- any(is.infinite(exp(full.Z %*% curr_u))) | any(is.na(exp(full.Z %*% curr_u)))
    if(inf.zu){
        stop("Infinite estimates of u")
    }
    
    #compute mu.vec using inverse link function                   
    mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))
                        
    # use y_bar as the sample mean and s_hat as the sample variance
    data_shat <- var(y)
    new.r <- dispersion
                        
    theta_diff <- rep(Inf, nrow(curr_theta))
    sigma_diff <- Inf

    #compute variance-covariance matrix G
    curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)
    G_inv <- computeG_inv(curr_G=curr_G)

    #variance as percentage of total variance                    
    init.res.var <- (data_shat - colSums(curr_sigma))
    init.var.comps <- c(curr_sigma[, 1], init.res.var)/data_shat
    names(init.var.comps) <- c(rownames(curr_sigma), "residual")
    var.comps <- init.var.comps
    
    conv.list <- list()
    iters <- 1
    
    meet.conditions <- !((all(theta_diff < theta.conv)) & (sigma_diff < theta.conv) | iters >= max.hit)
   
    while(meet.conditions){
        
        print(iters)
        conv.list[[paste0(iters)]] <- list("Iter"=iters, "Theta"=curr_theta, "Sigma"=curr_sigma,
                                           "Theta.Diff"=theta_diff, "Sigma.Diff" = sigma_diff,
                                           "Theta.Converged"=theta_diff < theta.conv,
                                           "Sigma.Converged"=sigma_diff < theta.conv)
        
        #compute all matrices - information about them found within the respective functions
        D <- computeD(mu=mu.vec)
        D_inv <- computeD_inv(D=D)
        y_star <- computey_star(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u)
        V <- computeV0(mu=mu.vec, r=new.r)
        W_inv <- computeW_inv(D=D_inv, V=V)
        W <- computeW(W_inv=W_inv)
        V_star <- computeV_star(full.Z=full.Z, curr_G=curr_G, W_inv=W_inv)
        V_star_inv <- computeV_star_inv(V_star=V_star)
        #V_partial <- computeV_partial(full.Z=full.Z, mu.vec=mu.vec, D=D, r=new.r, V=V, D_inv=D_inv, curr_u=curr_u)
        V_partial <- computeV_partial_alt(full.Z=full.Z, curr_G=curr_G)
        
        #---- First estimate variance components with Newton Raphson procedure ---#
        print(curr_sigma)
        score_sigma <- sigmaScore(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta)
        hessian_sigma <- sigmaHessian(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta)
        sigma_update <- singleNR(score_vec=score_sigma, hess_mat=hessian_sigma, theta_hat=curr_sigma)
        sigma_diff <- abs(sigma_update - curr_sigma)
        
        # update sigma, G, and G_inv
        curr_sigma <- sigma_update
        curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)
        G_inv <- computeG_inv(curr_G=curr_G)
        
        #---- Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u---####
        theta_update <- solve_equations(X=X, W=W, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star) 
        theta_diff <- abs(theta_update - curr_theta)
        
        # update B, u and mu_vec to determine new values of score and hessian matrices
        curr_theta <- theta_update
        rownames(curr_theta) <- c(colnames(X), colnames(full.Z))
        curr_beta <- curr_theta[colnames(X), , drop=FALSE]
        curr_u <- curr_theta[colnames(full.Z), , drop=FALSE]
        mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))
   
        #compute percentage variances
        res.var <- (data_shat - colSums(curr_sigma))
        curr_var.comps <- c(curr_sigma[, 1], res.var)/data_shat
        names(curr_var.comps) <- c(rownames(curr_sigma), "residual")
        var.comps <- curr_var.comps
           
        iters <- iters + 1
        meet.conditions <- !((all(theta_diff < theta.conv)) & (all((sigma_diff) < theta.conv))| iters >= max.hit)
    }
    
    converged <- ((all(theta_diff < theta.conv)) & (all(abs(sigma_diff) < theta.conv)))
    final.list <- list("FE"=curr_theta[colnames(X), ], "RE"=curr_theta[colnames(full.Z), ],
                       "Sigma"=curr_sigma,
                       "Theta.Converged"=theta_diff < theta.conv,
                       "Sigma.Converged"=sigma_diff < theta.conv,
                       "VarComp"=var.comps, 
                       "converged"=converged,
                       "Iters"=iters,
                       "Dispersion"=new.r, 
                       "Iterations"=conv.list)
    return(final.list)
}
                    
computeV_star <- function(full.Z=full.Z, curr_G=curr_G, W_inv=W_inv){
    V_star = full.Z %*% curr_G %*% t(full.Z) + W_inv
    return(V_star)
}

computeV_star_inv <- function(V_star=V_star){
    if(det(V_star) < 1e-10){ 
        V_star_inv <- ginv(V_star)
    } else{
        V_star_inv <- solve(V_star)}
    return(V_star_inv)
}

computey_star <- function(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u){
    y_star <- ((X %*% curr_beta) + (full.Z %*% curr_u)) + D_inv %*% (y - exp((X %*% curr_beta) + (full.Z %*% curr_u)))
    return(y_star)
}

sigmaScore <- function(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta){
    LHS <- -0.5*matrix.trace(V_star_inv %*% V_partial)
    RHS <- 0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% V_partial %*% V_star_inv %*% (y_star - X %*% curr_beta)
    return(LHS + RHS)
}

sigmaHessian <- function(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta) {
    LHS <- -0.5*matrix.trace(-V_star_inv %*% V_partial %*% V_star_inv %*% V_partial)
    RHS <- 0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% (-2*V_partial %*% V_star_inv %*% V_partial) %*% V_star_inv %*% (y_star - X %*% curr_beta)
    return(LHS + RHS)
}

computeV_partial <- function(full.Z=full.Z, mu.vec=mu.vec, D=D, r=new.r, V=V, D_inv=D_inv, curr_u=curr_u, random.levels=random.levels){
    
    d.vec <- -1/(mu.vec^2)
    d0 <- matrix(0L, ncol=length(mu.vec), nrow=length(mu.vec))
    diag(d0) <- d.vec
    
    v.vec <- (2*mu.vec)/r - 1
    v0 <- matrix(0L, ncol=length(mu.vec), nrow=length(mu.vec))
    diag(v0) <- v.vec
    
    db_dsigma0 <- matrix((apply((mapUtoIndiv2(full.Z, curr_u, random.levels=random.levels)), MARGIN = 1,
                         FUN=function(x){
                             (1/(length(x)-1) * 2 * sum(x))^-1
                         })), ncol=3)
    s0 <- matrix(0L, ncol=length(mu.vec), nrow=length(mu.vec))
    diag(s0) <- db_dsigma0
    
    D_deriv <- d0 %*% D %*% full.Z %*% s0
    V_deriv <- v0 %*% D %*% full.Z %*% s0
        
    V_partial <- full.Z %*% t(full.Z) + (D_deriv %*% V %*% D_inv) + (D_inv %*% V_deriv %*% D_inv) + (D_inv %*% V %*% D_deriv)
    return(V_partial)
}

computeV_partial_alt <- function(full.Z=full.Z, curr_G=curr_G){
    V_partial <- full.Z %*% t(full.Z)
    return(V_partial)
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

computeG_inv <- function(curr_G=curr_G){
    if(det(curr_G) < 1e-10){
        # use a generalized inverse if numerically singular
        # this will occur if there is a single random effect too
        G_inv <- ginv(curr_G)
    } else{
        G_inv <- solve(curr_G)
    }
    return(G_inv)
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
    
    if(det(hess_mat) < 1e-10){
        theta_new <- theta_hat - ginv(hess_mat) %*% score_vec
    } else{
        theta_new <- theta_hat - solve(hess_mat) %*% score_vec
    }
    
    return(theta_new)
}

computeW <- function(W_inv=W_inv){
    # now we have to take the inverse of W_inv
    if(det(W_inv) < 1e-10){ 
        W <- ginv(W_inv)
    } else{
        W <- solve(W_inv)}
    return(W)
}

computeW_inv <- function(D_inv=D_inv, V=V){
    W_inv = (D_inv %*% V %*% D_inv)
    return(W_inv)
}

computeV0 <- function(mu, r){
    # compute diagonal matrix of variances
    v.vec <- ((mu**2/r)) - mu
    V0 <- matrix(0L, ncol=length(mu), nrow=length(mu))
    diag(V0) <- v.vec
    return(V0)
}

computeD <- function(mu=mu.vec){
    # D is diag(mu_i)
    D <- matrix(0L, ncol=length(mu), nrow=length(mu))
    diag(D) <- mu
    return(D)
}

computeD_inv <- function(D){
    # Compute D^-1 from D
    # need to check that D is not singular
    if(det(D) < 1e-10){ 
        D_inv <- ginv(D)
    } else{
        D_inv <- solve(D)}
    return(D_inv)
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

mapUtoIndiv2 <- function(full.Z, curr_u, random.levels){
    # map the vector of random effects to the full nx1 vector
    rand.levels <- unlist(random.levels)
    j.G <- matrix(0L, ncol=nrow(full.Z), nrow=length(rand.levels))
    
    for(j in seq_along(rand.levels)){
        j.re <- rand.levels[j]
        j.b <- full.Z[, j.re, drop = FALSE] %*% curr_u[j.re, ]
        j.G[j ,] <- j.b
    }
    
    return(j.G)
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

glmmControl.defaults <- function(...){
    # return the default glmm control values
    return(list(det.tol=1e-10, cond.tol=1e-12, theta.tol=1e-6,
                likli.tol=1e-6, max.iter=100, lambda=1e-1, laplace.int="fe"))
}