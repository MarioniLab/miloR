### Components list for the GLMM

# dispersion estimation - approximate profile likelihood a la edgeR?
# estimable/predictable functions - LR, Wald-type or score test?
# random effects definitions/covariance definition - need to be able to provide a strict covariance matrix
# define likelihood functions
# implement and test both LMM and NB-LMM


# Newton-Raphson to estimate theta (beta and u), and tau (variance components)
# need theta_hat and tau_hat, will need to iterate between the two to acheive convergence?


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
                                      likli.tol=1e-6, max.iter=100, lambda=1e-1, laplace.int="fe")){
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
        curr_u <- matrix(rnorm(ncol(full.Z), mean=1, sd=1), ncol=1)
        rownames(curr_u) <- colnames(full.Z)

        curr_beta <- ginv((t(X) %*% X)) %*% t(X) %*% log(y + 1) # OLS for the betas is usually a good starting point for NR
        rownames(curr_beta) <- colnames(X)

    } else{
        curr_beta <- matrix(init.theta[["beta"]], ncol=1)
        rownames(curr_beta) <- colnames(X)

        curr_u <- matrix(init.theta[["rand"]] , ncol=1)
        rownames(curr_u) <- colnames(full.Z)
    }

    # initialise using solutions of variance components from ANOVA
    init.sigma <- matrix(estimateInitialSigmas(y, Z), ncol=1) # the last one is the residual variance component
    init.sigma <- init.sigma[-nrow(init.sigma), , drop=FALSE]
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
    s_hat <- var(y)
    r.val <- (y_bar **2)/(s_hat + y_bar) # methods of moments based estimate
    max.hit <- glmm.control[["max.iter"]]

    theta_diff <- rep(Inf, nrow(curr_theta))
    abs_diff <- abs(theta_diff)
    sigma_diff <- Inf
    loglihood.diff <- Inf
    lambda <- glmm.control[["lambda"]] # lambda used in the Levenberg-Marquardt adjustment
    init.vars <- runif(length(random.levels)) # sample variance components from ~U(0, 1)

    init.G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=init.sigma)
    curr_G <- init.G
    G_partials <- computeGPartials(curr_G, curr_sigma)

    # do a quick loglihood evaluation
    G_inv <- ginv(curr_G)
    init.loglihood  <- 0
    loglihood <- init.loglihood

    init.res.var <- (s_hat - sum(diag(init.G)))
    init.var.comps <- c(diag(init.G), init.res.var)/s_hat
    names(init.var.comps) <- c(rownames(init.G), "residual")

    var.comps <- init.var.comps
    conv.list <- list()
    iters <- 1

    meet.conditions <- !((all(abs_diff < theta.conv)) & (loglihood.diff < loglihood.eps) & (sigma_diff < theta.conv) | iters >= max.hit)

    while(meet.conditions){
        D_inv <- computeDinv(mu.vec)
        V_inv <- computeVinv(mu=mu.vec, y=y, r=r.val)
        B <- computeB(y=y, r=r.val, mu=mu.vec)
        W <- computeW(mu=mu.vec, r=r.val)
        Q <- computeQ(mu=mu.vec, r=r.val)

        # compute dG\dus
        Gu_partials <- computeGuPartials(curr_G=curr_G, u_hat=curr_u, cluster_levels=random.levels, sigmas=curr_sigma)
        score_beta <- betaScore(X=X, D_inv=D_inv, V_inv=V_inv, mu=mu.vec, y=y, r=r.val)
        score_u <- randScore(Z=full.Z, D_inv=D_inv, V_inv=V_inv, G_inv=G_inv, mu=mu.vec, y=y, r=r.val, u_hat=curr_u, Gu_partials=Gu_partials)

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
        ## Using Mike's appallingly shonky ANOVA-like approach (MASALA)
        beta_ss <- computeFixedSumSquares(curr_beta, X, y_bar)
        u_ss <- computeRandomSumSquares(curr_u, random.levels, full.Z, y_bar)
        res_ss <- computeResidualSimSquares(y, y_bar, beta_ss, u_ss)
        beta_df <- computeFixedDf(X)
        rownames(beta_df) <- colnames(X[, -1, drop=FALSE])
        u_df <- computeRandomDf(random.levels)
        rownames(u_df) <- names(random.levels)

        res_df <- computeResDf(length(y), beta_df, u_df)
        rownames(res_df) <- "residual"

        beta_ms <- computeFixedMeanSquares(beta_ss, beta_df)
        rownames(beta_ms) <- colnames(X[, -1, drop=FALSE])

        u_ms <- computeRandomMeanSquares(u_ss, u_df)
        rownames(u_ms) <- colnames(Z)

        res_ms <- res_ss/res_df
        rownames(res_ms) <- "residual"

        ss <- rbind(beta_ss, u_ss, res_ss)
        ms <- rbind(beta_ms, u_ms, res_ms)
        df <- rbind(beta_df, u_df, res_df)

        ems <- buildEMS(X, curr_beta, random.levels, y, df)
        # solve the system of equations
        ems.chol <- chol(ems)
        ems.solve <- backsolve(ems.chol, ms)
        rownames(ems.solve) <- rownames(ms)

        # although we don't use the fixed effects we should note that this
        # actually gives us the squared fixed effects
        sigma_update <- ems.solve[names(random.levels), , drop=FALSE]
        res_sigma <- ems.solve["residual", ]

        sigma_diff <- sigma_update - curr_sigma
        curr_sigma <- sigma_update
        mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))

        # update G with new sigmas
        curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)

        if(det(curr_G) < 1e-10){
            # use a generalized inverse if numerically singular
            G_inv <- ginv(curr_G)
        } else{
            G_inv <- solve(curr_G)
        }

        abs_diff <- abs(theta_diff)

        # need to sum all of the variances
        res.var <- (s_hat - colSums(curr_sigma))
        curr_var.comps <- c(curr_sigma[, 1], res.var)/s_hat
        names(curr_var.comps) <- c(rownames(curr_sigma), "residual")

        # loglihood integrating over the random effects only
        if(glmm.control$laplace.int %in% c("fe")){
            fe.hess <- betaHess(X=X, D_inv=D_inv, V_inv=V_inv, B=B, W=W, Q=Q) # this is used to integrate over FEs for var comp estimation
            curr.loglihood <- laplaceApprox(mu.vec, y, r.val, curr_G, G_inv, curr_u, fe.hess) # integrating over just the FEs
        } else if(glmm.control$laplace.int %in% c("full")){
            curr.loglihood <- laplaceApprox(mu.vec, y, r.val, curr_G, G_inv, curr_u, full.hess) # integrating over both the RE and FEs
        }

        loglihood.diff <- curr.loglihood - loglihood

        # the loglihood can bimble along at a small value but not converge with enough tolerance <- a local optimum?
        # message("Theta convergence: ", all(abs_diff < theta.conv))
        # message("Loglihood convergence:", abs(loglihood.diff) < loglihood.eps)
        # message("Delta Loglihood: ", abs(loglihood.diff))
        # message("VC convergence: ", all(sigma_diff < theta.conv))
        # print(iters >= max.hit)

        meet.conditions <- !((all(abs_diff < theta.conv)) & (abs(loglihood.diff) < loglihood.eps) &
                                 (all(sigma_diff < theta.conv))| iters >= max.hit)

        loglihood <- curr.loglihood
        var.comps.diff <- curr_var.comps - var.comps

        # with these variance components we can evaluate the _full_ likelihood
        full.loglihood <- nbLogLikelihood(mu=mu.vec, r=r.val, y=y) + normLogLikelihood(curr_G, G_inv, curr_u)

        curr.u_bars <- do.call(rbind, lapply(random.levels, FUN=function(X) mean(curr_theta[X, ])))

        bigV <- matrix(0L, ncol=length(mu.vec), nrow=length(mu.vec))
        v.vec <- (mu.vec**2 * (1/r.val)) - mu.vec
        diag(bigV) <- v.vec
        R <- bigV - (full.Z %*% curr_G %*% t(full.Z))

        conv.list[[paste0(iters)]] <- list("Iter"=iters, "Theta"=curr_theta, "Mu"=mu.vec, "Residual"=y - mu.vec, "Loglihood"=loglihood,
                                           "Hessian"=hess_theta, "Dispersion"=r.val, "Score"=full.score, "Theta.Diff"=theta_diff, "G"=curr_G,
                                           "Rand.Mean"=curr.u_bars, "Sigmas"=curr_sigma, #"Sigma.Hess"=hess_sigma,
                                           "V"=bigV, "R"=R,
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

    final.list <- list("FE"=curr_theta[colnames(X), ], "RE"=curr_theta[colnames(full.Z), ], "Loglihood"=loglihood,
                       "VarComp"=curr_var.comps,
                       "Sigmas"=curr_sigma, "Iterations"=conv.list)
    return(final.list)
}


buildEMS <- function(X, curr_beta, random.levels, y, df){
    # Construct a matrix N that will be used to map the
    # expected mean squares to the model parameters
    # order is FE, RE, residual
    tot.vars <- ncol(X) - 1 + length(random.levels) + 1
    M <- matrix(0L, nrow=tot.vars, ncol=tot.vars)
    M[, ncol(M)] <- 1 # the residual term is present for _all_ parameters

    for(j in seq_len(ncol(X)-1)){
        j.var <- colnames(X)[j+1]
        j.scalar <- apply(df[-which(j.var == rownames(df)), , drop=FALSE], 2, prod)

        M[j, j] <- j.scalar
    }

    # this currently doesn't handle interaction terms
    for(l in seq_len(length(random.levels))){
        l.var <- names(random.levels)[l]
        l.scalar <- apply(df[-which(l.var == rownames(df)), , drop=FALSE], 2, prod)
        l.idx <- ncol(X) - 1 + l
        M[l.idx, l.idx] <- l.scalar
    }

    return(M)
}



computeResidualSimSquares <- function(y, y_bar, beta_ss, u_ss){
    # compute the residual SS from the total - sum(beta_ss + u_ss)
    total_ss <- sum((y - y_bar)**2)
    residual_ss <- matrix(total_ss - (colSums(beta_ss) + colSums(u_ss)), ncol=1, nrow=1)
    return(residual_ss)
}


computeFixedSumSquares <- function(curr_beta, X, y_bar){
    # compute the fixed effect mean squares
    # X is a matrix, curr_beta is a column vector, y_bar is a scalar
    # note the betas are on a _log_ scale
    ss <- matrix(0L, ncol=1, nrow=ncol(X) - 1)
    for(j in seq_len(ncol(X)-1)){
        ss.j <- ((colSums(exp(curr_beta[c(1, j+1), , drop=FALSE])) - y_bar)**2)
        ss[j, ] <- ss.j
    }
    return(ss)
}


computeRandomSumSquares <- function(curr_u, random.levels, full.Z, y_bar){
    # compute the random effect mean squares
    # Z is a matrix, curr_u is a matrix with all of the estimated RE coefficients -
    # these are the group means for each level.
    # note the u's are on a _log_ scale
    ss <- matrix(0L, ncol=1, nrow=length(random.levels))
    for(k in seq_len(length(random.levels))){
        k.name <- names(random.levels)[k]
        k.levels <- random.levels[[k.name]]
        k.ss <- colSums((exp(curr_u[k.levels, , drop=FALSE]) - y_bar)**2)
        ss[k, ] <- k.ss
    }

    return(ss)
}


computeFixedDf <- function(X){
    df <- matrix(0L, ncol=1, nrow=ncol(X)-1)

    for(j in seq_len(ncol(X)-1)){
        if(is.factor(X[, j+1])){
            df.j <- length(levels(X[, j+1])) - 1
        } else if(is.character(X[, j+1])){
            df.j <- length(unique(X[, j+1])) - 1
        } else{
            df.j <- 1
        }

        df[j, ] <- df.j
    }
    return(df)
}


computeRandomDf <- function(random.levels){
    df <- matrix(0L, ncol=1, nrow=length(random.levels))
    re.levels <- names(random.levels)
    for(k in seq_len(length(re.levels))){
        k.re <- re.levels[k]
        if(length(random.levels[[k.re]]) > 1){
            df[k, ] <- length(random.levels[[k.re]]) - 1
        } else{
            df[k, ] <- 1
        }
    }
    return(df)
}


computeResDf <- function(n, beta_df, u_df){
    # compute the residual dfs
    return(matrix(n - (colSums(beta_df) + colSums(u_df)) - 1, ncol=1, nrow=1))
}


computeFixedMeanSquares <- function(beta_ss, beta_df){
    # compute the fixed effect mean squares, i.e. SS/df
    if(nrow(beta_ss) != nrow(beta_df)){
        stop("Rows are not concodrant")
    } else{
        return(beta_ss/beta_df)
    }
}


computeRandomMeanSquares <- function(u_ss, u_df){
    # compute the random effect mean squares
    if(nrow(u_ss) != nrow(u_df)){
        stop("Rows are not concodrant")
    } else{
        return(u_ss/u_df)
    }
}



computeGPartials <- function(curr_G, sigmas){
    # compute the partial derivatives dG/dsigma
    # return a list

    partial.list <- list()
    for(x in seq_len(nrow(sigmas))){
        partial.list[[x]] <- (curr_G == sigmas[x, ]) + 0
    }

    return(partial.list)
}

computeGuPartials <- function(curr_G, u_hat, cluster_levels, sigmas){
    # compute the partial derivatives dG/du
    # return a list
    n.re <- length(cluster_levels)
    partial.list <- list()
    for(x in seq_len(n.re)){
        x.G <- matrix(0L, ncol=ncol(curr_G), nrow=nrow(curr_G))
        x.G[curr_G == sigmas[x, ]] <- 2*u_hat[cluster_levels[[x]], ]
        partial.list[[x]] <- x.G
    }

    return(partial.list)
}


estimateInitialSigmas <- function(y, Z){
    ## compute the EMS given the input formula and solve for the ANOVA variance components
    ## the plan is to use VCA here - it'll save trying to reimplement it myself
    require(VCA)

    design <- as.data.frame(Z)
    design$y <- y
    form <- as.formula(paste("y ~ ", paste(colnames(Z), collapse=" + ")))
    var.comps <- anovaVCA(form, design)$VCoriginal
    names(var.comps) <- c(colnames(Z), "residual")
    return(var.comps)
}


computeQuadRoots <- function(X, curr_theta, r, S_hat){
    # solve the quadratic formula to get the potential roots to solve for sigmas
    a = exp(2 * X %*% curr_theta[colnames(X), ])
    b = -exp(2 * X %*% curr_theta[colnames(X), ])
    c = (1/r) * (exp(2 * X %*% curr_theta[colnames(X),]) - (r * exp(X %*% curr_theta[colnames(X),]))) - S_hat

    x_1 = (-b + sqrt(b**2 - 4*a*c))/(2 * a)
    x_2 = (-b - sqrt(b**2 - 4*a*c))/(2 * a)

    return(list(x_1, x_2))
}


sigmaSolve <- function(roots, Z){
    # Z is a vector now
    sigma <- log(roots)/as.vector(t(Z) %*% Z)
    return(sigma)
}



initialiseG <- function(Z, cluster_levels, sigmas){
    # construct the correct size of G given the random effects and variance components
    # names of cluster_levels and columns of Z must match
    # the independent sigmas go on the diagonal and the off-diagonal are the crossed/interactions
    # sigmas must be named
    sum.levels <- sum(unlist(lapply(cluster_levels, length)))
    G <- sparseMatrix(i=sum.levels, j=sum.levels, repr="C", x=0L)
    i <- j <- 1

    for(x in seq_len(nrow(sigmas))){
        x.q <- length(cluster_levels[[rownames(sigmas)[x]]])
        diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ]
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
            i.mod <- all(Z[, i, drop=FALSE] %% 1) == 0
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
    y_diff <- y - (n*mu) - (n*r)/(1 - (r*(mu**-1)))
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


varScore <- function(G_inv, u_hat, G_partials, n.comps){
    # compute the score function for the sigmas
    svec <- matrix(0L, nrow=n.comps, ncol=1)
    for(i in seq_len(n.comps)){
        svec[i,] <- (-0.5 * matrix.trace(G_inv %*% G_partials[[i]])) + (0.5 * t(u_hat) %*% G_inv %*% G_partials[[i]] %*% G_inv %*% u_hat)
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


laplaceApprox <- function(mu, y, r, G, G_inv, curr_u, hessian){
    ## compute the Laplace approximation to the marginal loglihood for the variance components
    nb.liklihood <- nbLogLikelihood(mu=mu, r=r, y=y)
    norm.liklihood <- normLogLikelihood(G, G_inv, curr_u)

    det.hess <- det(hessian)
    if(det.hess > 0){
        log.det.hess <- log(det.hess)
    } else{
        log.det.hess <- 0
    }
    nb.liklihood + norm.liklihood - ((1/2) * (log.det.hess/(2*pi)))
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

computeModes <- function(x) {
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)]
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
