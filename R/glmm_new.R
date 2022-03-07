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
    curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
    rownames(curr_u) <- colnames(full.Z)

    # OLS for the betas is usually a good starting point for NR
    curr_beta <- solve((t(X) %*% X)) %*% t(X) %*% log(y + 1)
    rownames(curr_beta) <- colnames(X)

    # compute sample variances of the us
    curr_sigma <- Matrix(lapply(lapply(mapUtoIndiv(full.Z, curr_u, random.levels=random.levels),
                                       FUN=function(Bj){
                                           (1/(length(Bj)-1)) * crossprod(Bj, Bj)
                                           }), function(y){attr(y, 'x')}), ncol=1, sparse=TRUE)
    # curr_sigma <- matrix(unlist(lapply(mapUtoIndiv(full.Z, curr_u, random.levels=random.levels),
    #                                    FUN=function(Bj){
    #                                        (1/(length(Bj)-1)) * crossprod(Bj, Bj)
    #                                    })), ncol=1)

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
    G_inv <- solve(curr_G)

    conv.list <- list()
    iters <- 1
    meet.conditions <- !((all(theta_diff < theta.conv)) & (sigma_diff < theta.conv) | iters >= max.hit)

    while(meet.conditions){

        print(iters)
        conv.list[[paste0(iters)]] <- list("Iter"=iters, "Theta"=curr_theta, "Sigma"=curr_sigma,
                                           "Theta.Diff"=theta_diff, "Sigma.Diff" = sigma_diff,
                                           "Theta.Converged"=theta_diff < theta.conv,
                                           "Sigma.Converged"=sigma_diff < theta.conv)

        #compute all matrices - information about them found within their respective functions
        D <- computeD(mu=mu.vec)
        D_inv <- solve(D)
        y_star <- computey_star(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u, y=y)
        V <- computeV(mu=mu.vec, r=new.r)
        W <- computeW(D_inv=D_inv, V=V)
        W_inv <- solve(W)
        V_star <- computeV_star(full.Z=full.Z, curr_G=curr_G, W=W)
        V_star_inv <- solve(V_star)
        V_partial <- computeV_partial(full.Z=full.Z, random.levels=random.levels, curr_sigma=curr_sigma)

        #---- First estimate variance components with Newton Raphson procedure ---#
        if (isFALSE(REML)) {
            score_sigma <- sigmaScore(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, random.levels=random.levels)
            information_sigma <- sigmaInformation(V_star_inv=V_star_inv, V_partial=V_partial, random.levels=random.levels)
        } else if (isTRUE(REML)) {
            P <- computeP_REML(V_star_inv=V_star_inv, X=X)
            score_sigma <- sigmaScoreREML(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, P=P, random.levels=random.levels)
            information_sigma <- sigmaInformationREML(V_star_inv=V_star_inv, V_partial=V_partial, P=P, random.levels=random.levels)
        }
        sigma_update <- FisherScore(score_vec=score_sigma, hess_mat=information_sigma, theta_hat=curr_sigma, random.levels=random.levels)
        sigma_diff <- abs(sigma_update - curr_sigma)

        # update sigma, G, and G_inv
        curr_sigma <- sigma_update
        curr_G <- initialiseG(full.Z, cluster_levels=random.levels, sigmas=curr_sigma)
        G_inv <- solve(curr_G)

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

    SE <- calculateSE(X=X, full.Z=full.Z, W_inv=W_inv, G_inv=G_inv)
    Zscore <- calculateZscore(curr_beta=curr_beta, SE=SE)
    df <- Satterthwaite_df(X=X, SE=SE, REML=REML, W_inv=W_inv, full.Z=full.Z, curr_sigma=curr_sigma, curr_beta=curr_beta, random.levels=random.levels, V_partial=V_partial, V_star_inv=V_star_inv, G_inv=G_inv)
    Pvalue <- computePvalue(Zscore=Zscore, df=df)

    converged <- ((all(theta_diff < theta.conv)) & (all(abs(sigma_diff) < theta.conv)))
    final.list <- list("FE"=curr_theta[colnames(X), ], "RE"=curr_theta[colnames(full.Z), ],
                       "Sigma"=matrix(curr_sigma),
                       "Theta.Converged"=theta_diff < theta.conv,
                       "Sigma.Converged"=sigma_diff < theta.conv,
                       "converged"=converged,
                       "Iters"=iters,
                       "Dispersion"=new.r,
                       "SE"=SE,
                       "Zscore"=Zscore,
                       "Pvalue"=Pvalue,
                       "Iterations"=conv.list)
    return(final.list)
}

computeW <- function(D_inv=D_inv, V=V){
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

computeD <- function(mu=mu.vec){
    # D is diag(mu_i)
    D <- Matrix(0L, ncol=length(mu), nrow=length(mu), sparse = TRUE)
    diag(D) <- mu
    return(D)
}

computeV_star <- function(full.Z=full.Z, curr_G=curr_G, W=W){
    V_star = full.Z %*% curr_G %*% t(full.Z) + W
    return(V_star)
}

computey_star <- function(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u, y=y){
    y_star <- ((X %*% curr_beta) + (full.Z %*% curr_u)) + D_inv %*% (y - exp((X %*% curr_beta) + (full.Z %*% curr_u)))
    return(y_star)
}

computeV_partial <- function(full.Z=full.Z, random.levels=random.levels, curr_sigma=curr_sigma){
    V_partial_vec <- list()
    j <- 1
    for (i in random.levels) {
        Z.temp <- full.Z[ , i]
        V_partial_vec[[j]] <- Z.temp %*% t(Z.temp)
        j <- j + 1
    }
    names(V_partial_vec) <- rownames(curr_sigma)
    return(V_partial_vec)
}

sigmaScore <- function(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, random.levels=random.levels){
    score_vec <- list()
    count <- 1
    for (k in names(random.levels)){
      temp_Vpartial <- V_partial[grep(k, names(V_partial))]
      for (i in 1:length(temp_Vpartial)) {
        LHS <- -0.5*matrix.trace(V_star_inv %*% temp_Vpartial[[i]])
        RHS <- 0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% temp_Vpartial[[i]] %*% V_star_inv %*% (y_star - X %*% curr_beta)
        score_vec[[count]]<- LHS + RHS
        count <- count + 1
      }
    }
    return(score_vec)
}

sigmaInformation <- function(V_star_inv=V_star_inv, V_partial=V_partial, random.levels=random.levels) {
    info_vec <- list()
    temp_Vpartial <- list()
    count <- 1
    for (k in names(random.levels)){
      temp_Vpartial <- V_partial[grep(k, names(V_partial))]
      temp_mat <- matrix(NA, nrow=length(temp_Vpartial), ncol=length(temp_Vpartial))
      for (i in 1:length(temp_Vpartial)) {
        for (j in 1:length(temp_Vpartial)) {
          temp_mat[i,j] <- 0.5*matrix.trace(V_star_inv %*% temp_Vpartial[[i]] %*% V_star_inv %*% temp_Vpartial[[j]])
          }
      }
      info_vec[[count]] <- temp_mat
      count <- count + 1
    }
    return(info_vec)
}

sigmaScoreREML <- function(V_star_inv=V_star_inv, V_partial=V_partial, y_star=y_star, X=X, curr_beta=curr_beta, P=P, random.levels=random.levels){
    score_vec <- NA
    for (i in 1:length(random.levels)) {
        LHS <- -0.5*matrix.trace(P %*% V_partial[[i]])
        RHS <- 0.5*t(y_star - X %*% curr_beta) %*% V_star_inv %*% V_partial[[i]] %*% V_star_inv %*% (y_star - X %*% curr_beta)
        score_vec[i] <- LHS + RHS
    }
    return(score_vec)
}

sigmaInformationREML <- function(V_star_inv=V_star_inv, V_partial=V_partial, P=P, random.levels=random.levels) {
    info_vec <- NA
    info_vec <- sapply(1:length(random.levels), function(i){
        0.5*matrix.trace(P %*% V_partial[[i]] %*% P %*% V_partial[[i]])})
    return(info_vec)
}

computeP_REML <- function(V_star_inv=V_star_inv, X=X) {
   P <- V_star_inv - V_star_inv %*% X %*% solve(t(X) %*% V_star_inv %*% X) %*% t(X) %*% V_star_inv
   return(P)
}

FisherScore <- function(score_vec, hess_mat, theta_hat, random.levels, lambda=1e-5, det.tol=1e-10, cond.tol=1e-15){
    # sequentially update the parameter using the Newton-Raphson algorithm
    # theta ~= theta_hat + hess^-1 * score
    # this needs to be in a direction of descent towards a minimum
    theta_new <- NA
    for (i in 1:length(random.levels)) {
      theta_new[i] <- theta_hat[i] + solve(hess_mat[[i]]) * score_vec[[i]]
    }

    theta_new <- Matrix(theta_new, sparse = TRUE)
    rownames(theta_new) <- rownames(theta_hat)

    # if(det(hess_reformat) < 1e-10){
    #     theta_new <- theta_hat + ginv(hess_reformat) %*% score_reformat
    # } else{
    #     theta_new <- theta_reformat + solve(hess_reformat) %*% score_reformat
    #     theta_new <- diag(theta_new)
    # }

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
    full.Z <- Matrix(full.Z, sparse = TRUE)
    return(full.Z)
}

solve_equations <- function(X=X, W_inv=W_inv, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star){

    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% full.Z
    LowerLeft <- t(full.Z) %*% W_inv %*% X
    LowerRight <- t(full.Z) %*% W_inv %*% full.Z + G_inv

    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W_inv %*% y_star), (t(full.Z) %*% W_inv %*% y_star))

    theta_update <- solve(LHS) %*% RHS
    return(theta_update)
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

calculateSE <- function(X=X, full.Z=full.Z, W_inv=W_inv, G_inv=G_inv) {
    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% full.Z
    LowerLeft <- t(full.Z) %*% W_inv %*% X
    LowerRight <- t(full.Z) %*% W_inv %*% full.Z + G_inv
    se <- sqrt(diag(solve(UpperLeft - UpperRight %*% solve(LowerRight) %*% LowerLeft)))
    return(se)
}

calculateZscore <- function(curr_beta=curr_beta, SE=SE) {
    Zscore <- as.matrix(curr_beta)/SE
    return(Zscore)
}

computePvalue <- function(Zscore=Zscore, df=df) {
    pval <- 2*pt(abs(Zscore), df, lower.tail=FALSE)
    print(pval)
    return(pval)
}

Satterthwaite_df <- function(X=X, SE=SE, REML=REML, W_inv=W_inv, full.Z=full.Z, curr_sigma=curr_sigma, curr_beta=curr_beta, random.levels=random.levels, V_partial=V_partial, V_star_inv=V_star_inv, G_inv=G_inv) {

  ###---- first calculate g = derivative of C with respect to sigma ----
    function_jac <- function(x, X.fun=as.matrix(X), W_inv.fun=as.matrix(W_inv), full.Z.fun=as.matrix(full.Z)) {
      UpperLeft <- t(X.fun) %*% W_inv.fun %*% X.fun
      UpperRight <- t(X.fun) %*% W_inv.fun %*% full.Z.fun
      LowerLeft <- t(full.Z.fun) %*% W_inv.fun %*% X.fun
      LowerRight <- t(full.Z.fun) %*% W_inv.fun %*% full.Z.fun
      n <- length(random.levels)
      diag(LowerRight) <- diag(LowerRight) + rep(1/x, times=lengths(random.levels)) #when extending to random slopes, this needs to be changed to a matrix and added to LowerRight directly
      C <- solve(UpperLeft - UpperRight %*% solve(LowerRight) %*% LowerLeft)
    }

    jac <- numDeriv::jacobian(func=function_jac, x=as.vector(curr_sigma))
    jac_list <- lapply(1:ncol(jac), function(i)
      array(jac[, i], dim=rep(length(curr_beta), 2))) #when extending to random slopes, this would have to be reformatted into list, where each element belongs to one random effect

    #next, calculate V_a, the asymptotic covariance matrix of the estimated covariance parameters
    #given by formula below
    P <- computeP_REML(V_star_inv=V_star_inv, X=X)
    V_a <- list()
    count <- 1
    for (k in names(random.levels)){
      temp_Vpartial <- V_partial[grep(k, names(V_partial))]
      for (i in 1:length(temp_Vpartial)) {
        for (j in 1:length(temp_Vpartial)) {
          V_a[[count]] <- 2*(1/(matrix.trace(P %*% temp_Vpartial[[i]] %*% P %*% temp_Vpartial[[j]])))
          count <- count + 1
        }
      }
    }
    V_a <- bdiag(V_a)

    # P <- computeP_REML(V_star_inv=V_star_inv, X=X)
    # V_a <- matrix(NA, nrow=length(random.levels), ncol=length(random.levels))
    # for (i in 1:length(random.levels)) {
    #   for (j in 1:length(random.levels)) {
    #     V_a[i,j] <- 2*1/(matrix.trace(P %*% V_partial[[i]] %*% P %*% V_partial[[j]]))
    #   }
    # }
    # print(V_a)

    df <- rep(NA, length(curr_beta))
    for (i in 1:length(curr_beta)) {
      jac_var_beta <- unlist(lapply(lapply(jac_list, diag), `[[`, i))
      denom <- t(jac_var_beta) %*% (V_a) %*% jac_var_beta #g' Va g
      df[i] <- 2*((SE[i]^2)^2)/denom
    }
    return(as.matrix(df))
}
