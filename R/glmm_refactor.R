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
runGLMM <- function(X, Z, y, init.theta=NULL, crossed=FALSE, random.levels=NULL, REML=FALSE,
                    glmm.control=list(theta.tol=1e-6, max.iter=100),
                    dispersion = 0.5){

    # model components
    # X - fixed effects model matrix
    # Z - random effects model matrix
    # A - genetic relationship matrix
    # y - observed phenotype

    theta.conv <- glmm.control[["theta.tol"]] # convergence for the parameters
    max.hit <- glmm.control[["max.iter"]]

    # create full Z with expanded random effect levels
    full.Z <- initializeFullZ(Z=Z, cluster_levels=random.levels)

    # random value initiation from runif
    curr_u <- Matrix(runif(ncol(full.Z), 0, 1), ncol=1)
    rownames(curr_u) <- colnames(full.Z)

    # OLS for the betas is usually a good starting point for NR
    curr_beta <- Matrix::solve((t(X) %*% X)) %*% t(X) %*% log(y + 1)
    rownames(curr_beta) <- colnames(X)

    # compute sample variances of the us
    curr_sigma <- Matrix(lapply(lapply(mapUtoIndiv(full.Z, curr_u, random.levels=random.levels),
                                       FUN=function(Bj){
                                           (1/(length(Bj)-1)) * crossprod(Bj, Bj)
                                       }), function(y){attr(y, 'x')}), ncol=1, sparse=TRUE)
    rownames(curr_sigma) <- colnames(Z)

    #create a single variable for the thetas
    curr_theta <- do.call(rbind, list(curr_beta, curr_u))

    #compute mu.vec using inverse link function
    mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))

    # use y_bar as the sample mean and s_hat as the sample variance
    new.r <- dispersion

    theta_diff <- rep(Inf, nrow(curr_theta))
    sigma_diff <- Inf

    #compute variance-covariance matrix G
    curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)
    G_inv <- computeInv(curr_G)

    conv.list <- list()
    iters <- 1
    meet.conditions <- !((all(theta_diff < theta.conv)) & (sigma_diff < theta.conv) | iters >= max.hit)

    while(meet.conditions){

        conv.list[[paste0(iters)]] <- list("Iter"=iters, "Theta"=curr_theta, "Sigma"=curr_sigma,
                                           "Theta.Diff"=theta_diff, "Sigma.Diff" = sigma_diff,
                                           "Theta.Converged"=theta_diff < theta.conv,
                                           "Sigma.Converged"=sigma_diff < theta.conv)

        #compute all matrices - information about them found within their respective functions
        D <- computeD(mu=mu.vec)
        D_inv <- computeInv(D)
        y_star <- computey_star(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u, y=y)
        V <- computeV(mu=mu.vec, r=new.r)
        W <- computeW(D_inv=D_inv, V=V)
        W_inv <- computeInv(W)
        V_star <- computeV_star(full.Z=full.Z, curr_G=curr_G, W=W)
        V_star_inv <- computeVstar_inverse(full.Z=full.Z, curr_G=curr_G, W_inv=W_inv)
        V_partial <- computeV_partial(full.Z=full.Z, random.levels=random.levels)

        # precompute all of the necessary matrices
        matrix_list <- preComputeMatrices(V_star_inv=V_star_inv, V_partial=V_partial, X=X,
                                          curr_beta=curr_beta, full.Z=full.Z, curr_u=curr_u, y_star=y_star)

        #---- First estimate variance components with Newton Raphson procedure ---#
        if (isFALSE(REML)) {
            score_sigma <- sigmaScore(matrix_list=matrix_list, V_star_inv=V_star_inv, random.levels=random.levels)
            information_sigma <- sigmaInformation(V_star_inv=V_star_inv, V_partial=V_partial, random.levels=random.levels)
        } else if (isTRUE(REML)) {
            P <- computeP_REML(V_star_inv=V_star_inv, X=X) # this is the other bottleneck - why?
            matrix_list[["PVSTARi"]] <- lapply(V_partial, function(i) Matrix::crossprod(P, i))
            score_sigma <- sigmaScoreREML(matrix_list=matrix_list, V_star_inv=V_star_inv, y_star=y_star, P=P, random.levels=random.levels)
            information_sigma <- sigmaInformationREML(matrix_list=matrix_list, random.levels=random.levels)
        }
        sigma_update <- FisherScore(score_vec=score_sigma, hess_mat=information_sigma, theta_hat=curr_sigma)
        sigma_diff <- abs(sigma_update - curr_sigma)

        # update sigma, G, and G_inv
        curr_sigma <- sigma_update
        curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)
        G_inv <- computeInv(curr_G)

        #---- Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u---####
        theta_update <- solve_equations(X=X, W_inv=W_inv, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star)
        theta_diff <- abs(theta_update - curr_theta)

        # update B, u and mu_vec to determine new values of score and hessian matrices
        curr_theta <- theta_update
        rownames(curr_theta) <- c(colnames(X), colnames(full.Z))
        curr_beta <- curr_theta[colnames(X), , drop=FALSE]
        curr_u <- curr_theta[colnames(full.Z), , drop=FALSE]
        mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))

        iters <- iters + 1
        meet.conditions <- !((all(theta_diff < theta.conv)) & (all((sigma_diff) < theta.conv))| iters >= max.hit)
    }

    converged <- ((all(theta_diff < theta.conv)) & (all(abs(sigma_diff) < theta.conv)))
    final.list <- list("FE"=curr_theta[colnames(X), ], "RE"=curr_theta[colnames(full.Z), ],
                       "Sigma"=matrix(curr_sigma),
                       "Theta.Converged"=theta_diff < theta.conv,
                       "Sigma.Converged"=sigma_diff < theta.conv,
                       "converged"=converged,
                       "Iters"=iters,
                       "Dispersion"=new.r,
                       "VSTAR"=V_star,
                       "G"=curr_G,
                       "Hessian"=information_sigma,
                       "Iterations"=conv.list)
    return(final.list)
}

computeW <- function(D_inv, V){
    W = D_inv %*% V %*% D_inv
    return(W)
}

computeV <- function(mu, r){
    # compute diagonal matrix of variances
    v.vec <- ((mu**2/r)) + mu
    V <- Matrix(0L, ncol=length(mu), nrow=length(mu), sparse = TRUE)
    diag(V) <- v.vec
    return(V)
}

computeD <- function(mu){
    # D is diag(mu_i)
    D <- Matrix(0L, ncol=length(mu), nrow=length(mu), sparse = TRUE)
    diag(D) <- mu
    return(D)
}

computeV_star <- function(full.Z, curr_G, W){
    V_star <- full.Z %*% Matrix::crossprod(curr_G, t(full.Z)) + W
    return(V_star)
}

computey_star <- function(X, curr_beta, full.Z, D_inv, curr_u, y){
    y_star <- ((X %*% curr_beta) + (full.Z %*% curr_u)) + D_inv %*% (y - exp((X %*% curr_beta) + (full.Z %*% curr_u)))
    return(y_star)
}

computeV_partial <- function(full.Z, random.levels){
    V_partial_vec <- list()
    j <- 1
    for (i in seq_along(random.levels)) {
        Z.temp <- full.Z[ , random.levels[[i]]]
        V_partial_vec[[j]] <- Z.temp %*% t(Z.temp)
        j <- j + 1
    }
    return(V_partial_vec)
}


computeVstar_inverse <- function(full.Z, curr_G, W_inv){
    # compute the inverse of V_star using Henderson-adjusted Woodbury formula, equation (18)
    # (A + UBU^T)^-1 = A^-1 - A^-1UB[I + U^TA^-1UB]^-1U^TA^-1
    # Only requires A^-1, where B = ZGZ^T, A=W, U=Z
    I <- Matrix(0L, nrow=ncol(full.Z), ncol=ncol(full.Z))
    diag(I) <- 1

    l.1 <- W_inv %*% full.Z
    left.p <-  l.1 %*% curr_G
    mid.1 <- t(full.Z) %*% W_inv
    mid.2 <- mid.1 %*% full.Z
    mid.inv <- Matrix::solve(I + mid.2 %*% curr_G)

    return(W_inv - (left.p %*% mid.inv %*% t(full.Z) %*% W_inv))
}


preComputeMatrices <- function(V_star_inv, V_partial, X, curr_beta, full.Z, curr_u, y_star){
    # precompute certain matrices from matrix multiplications that are needed > once
    mat.list <- list()
    mat.list[["XBETA"]] <- X %*% curr_beta
    mat.list[["ZU"]] <- full.Z %*% curr_u
    mat.list[["VSTARDi"]] <- lapply(V_partial, FUN=function(i) V_star_inv %*% i) # this is also a list of matrices
    mat.list[["YSTARMINXB"]] <- y_star - mat.list[["XBETA"]]
    mat.list[["XTVSTAR"]] <- t(X) %*% V_star_inv
    mat.list[["VSTARX"]] <- V_star_inv %*% X

    return(mat.list)
}


sigmaScore <- function(matrix_list, V_star_inv, random.levels){
    score_vec <- NA
    for (i in seq_along(random.levels)) {
        LHS <- -0.5*matrix.trace((matrix_list[["VSTARDi"]])[[i]])
        rhs.1 <- t(matrix_list[["YSTARMINXB"]]) %*% (matrix_list[["VSTARDi"]])[[i]]
        rhs.2 <- rhs.1 %*% V_star_inv
        RHS <- 0.5* rhs.2 %*% (matrix_list[["YSTARMINXB"]])
        score_vec[i] <- LHS + RHS
    }
    return(score_vec)
}

sigmaInformation <- function(V_star_inv, V_partial, random.levels) {

    sigma_info <- Matrix(0L, ncol=length(V_partial), nrow=length(V_partial))

    for(i in seq_along(V_partial)){
        for(j in seq_along(V_partial)){
            inner.1 <- V_star_inv %*% V_partial[[i]]
            inner.2 <- inner.1 %*% V_star_inv
            sigma_info[i, j] <- 0.5*matrix.trace(inner.2 %*% V_partial[[j]])
        }
    }
    return(sigma_info)
}


## some of these matrix multiplications are used multiple times throughout the code - perhaps we should store these to reduce the number?

sigmaScoreREML <- function(matrix_list, V_star_inv, y_star, P, random.levels){
    score_vec <- Matrix(0L, ncol=1, nrow=length(random.levels), sparse=FALSE)

    for (i in seq_along(random.levels)) {
        LHS <- -0.5 * matrix.trace(matrix_list[["PVSTARi"]][[i]])
        rhs.1 <- t(y_star) %*% matrix_list[["PVSTARi"]][[i]]
        rhs.2 <- rhs.1 %*% P
        RHS <- 0.5 * (rhs.2 %*% y_star)
        score_vec[i, ] <- LHS + RHS
    }

    return(score_vec)
}


sigmaInformationREML <- function(matrix_list, random.levels) {
    # this should be a matrix
    sigma_info <- Matrix(0L, ncol=length(random.levels), nrow=length(random.levels))

    for(i in seq_along(random.levels)){
        for(j in seq_along(random.levels)){
            sigma_info[i, j] <- 0.5*matrix.trace(Matrix::crossprod(matrix_list[["PVSTARi"]][[i]], matrix_list[["PVSTARi"]][[j]]))
        }
    }

    return(sigma_info)
}

computeP_REML <- function(V_star_inv, X) {
    # breaking these down to individual steps speeds up the operations considerably
    tx.m <- t(X) %*% V_star_inv
    x.inv <- computeInv(tx.m %*% X)
    tx.inv <- t(X) %*% V_star_inv
    vx <- V_star_inv %*% X
    Pminus <- vx %*% x.inv
    P <- V_star_inv - Pminus %*% tx.inv

    return(P)
}

FisherScore <- function(score_vec, hess_mat, theta_hat, lambda=1e-5, det.tol=1e-10, cond.tol=1e-15){
    # sequentially update the parameter using the Newton-Raphson algorithm
    # theta ~= theta_hat + hess^-1 * score
    # this needs to be in a direction of descent towards a minimum


    theta_new <- tryCatch({
        theta_hat + Matrix::solve(hess_mat) %*% score_vec
    }, error=function(cond){
        message("Hessian is singular. Original error message:")
        error(cond)
        return(NULL)
    }, finally={

    })

    return(theta_new)
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
    for(x in seq_len(nrow(sigmas))){
        x.q <- length(cluster_levels[[rownames(sigmas)[x]]])
        diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ] # is this sufficient to transform the sigma to the model scale?
        i <- j <- i+x.q
    }
    return((G))
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
            i.levels <- as.factor(unique(Z[, i, drop=FALSE])) # ordering is arbitrary
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
            print(head(i.z))
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
        print(head(i.z))
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
    full.Z <- Matrix(full.Z, sparse = TRUE)
    return(full.Z)
}

solve_equations <- function(X, W_inv, full.Z, G_inv, curr_beta, curr_u, y_star){

    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% full.Z
    LowerLeft <- t(full.Z) %*% W_inv %*% X
    LowerRight <- t(full.Z) %*% W_inv %*% full.Z + G_inv

    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W_inv %*% y_star), (t(full.Z) %*% W_inv %*% y_star))

    theta_update <- Matrix::solve(LHS) %*% RHS
    return(theta_update)
}

mapUtoIndiv <- function(full.Z, curr_u, random.levels){
    # map the vector of random effects to the full nx1 vector
    rand.levels <- names(random.levels)
    indiv.u.list <- list()

    for(j in seq_along(rand.levels)){
        j.G <- Matrix(0L, ncol=nrow(full.Z), nrow=nrow(full.Z))
        j.re <- rand.levels[j]
        j.levels <- random.levels[[j.re]]
        j.b <- full.Z[, j.levels] %*% curr_u[j.levels, ]
        indiv.u.list[[j.re]] <- j.b
    }

    return(indiv.u.list)
}

### utility functions
computeInv <- function(x){
    # Compute x^-1 from x
    # need to check that x is not singular - use tryCatch - if matrix is singular then report error message

    x_inv <- tryCatch(expr={
        Matrix::solve(x)
    },
    error=function(cond){
        message("Matrix cannot be inverted - most likely singular")
        message(cond)
        return(NULL)
    },
    finally={

    })

    return(x_inv)
}

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
    return(list(theta.tol=1e-6, max.iter=100))
}
