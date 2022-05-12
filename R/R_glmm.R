#' Perform differential abundance testing using a NB-generalised linear mixed model
#'
#' This function will perform DA testing per-nhood using a negative binomial generalised linear mixed model
#'
#' @param X A matrix containing the fixed effects of the model.
#' @param Z A matrix containing the random effects of the model.
#' @param y A matrix containing the observed phenotype over each neighborhood.
#' @param offsets A vector containing the (log) offsets to apply normalisation for different numbers of cells across samples.
#' @param init.theta A column vector (m X 1 matrix) of initial estimates of fixed and random effect coefficients
#' @param Kin A n x n covariance matrix to explicitly model variation between observations
#' @param REML A logical value denoting whether REML (Restricted Maximum Likelihood) should be run. Default is TRUE.
#' @param random.levels A list describing the random effects of the model, and for each, the different unique levels.
#' @param glmm.control A list containing parameter values specifying the theta tolerance of the model and the maximum number of iterations to be run.
#' @param dispersion A scalar value for the dispersion of the negative binomial.
#' @param geno.only A logical value that flags the model to use either just the \code{matrix} `Kin` or the supplied random effects.
#'
#' @details
#' This function runs a negative binomial generalised linear mixed effects model. If mixed effects are detected in testNhoods,
#' this function is run to solve the model.
#'
#'
#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix solve crossprod
#' @export
fitGLMM <- function(X, Z, y, offsets, init.theta=NULL, Kin=NULL,
                    random.levels=NULL, REML=FALSE,
                    glmm.control=list(theta.tol=1e-6, max.iter=100),
                    dispersion = 0.5, geno.only=FALSE){

    # model components
    # X - fixed effects model matrix
    # Z - random effects model matrix
    # y - observed phenotype

    # check all dimensions conform
    if(nrow(X) != nrow(Z) | nrow(X) != length(y)){
        stop("Dimensions of y, X and Z are discordant. y: ", length(y), "x1, X:",
             nrow(X), "x", ncol(X), ", Z:", nrow(Z), "x", ncol(Z))
    }

    theta.conv <- glmm.control[["theta.tol"]] # convergence for the parameters
    max.hit <- glmm.control[["max.iter"]]

    # OLS for the betas is usually a good starting point for NR
    curr_beta <- solve((t(X) %*% X)) %*% t(X) %*% log(y + 1)
    rownames(curr_beta) <- colnames(X)
    
    # create full Z with expanded random effect levels
    full.Z <- initializeFullZ(Z, cluster_levels = random.levels)
    colnames(full.Z) <- unlist(random.levels)

    # sample random value for RE us
    curr_u <- Matrix(runif(ncol(full.Z), 0, 1), ncol=1, sparse = TRUE)

    # sample variances of the us
    curr_sigma <- Matrix(runif(ncol(Z), 0, 1), ncol=1, sparse = TRUE)
    rownames(curr_sigma) <- colnames(Z)

    u_indices <- sapply(seq_along(random.levels),
                        FUN=function(RX) which(colnames(full.Z) %in% random.levels[[RX]]),
                        simplify=FALSE)

    #create a single variable for the thetas
    curr_theta <- do.call(rbind, list(curr_beta, curr_u))

    #compute mu.vec using inverse link function
    mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))

    theta_diff <- rep(Inf, nrow(curr_theta))
    sigma_diff <- Inf

    #compute variance-covariance matrix G
    curr_G <- initialiseG(cluster_levels=random.levels, sigmas=curr_sigma)
    G_inv <- computeInv(curr_G)

    conv.list <- list()
    iters <- 1
    meet.conditions <- !((all(theta_diff < theta.conv)) & (sigma_diff < theta.conv) | iters >= max.hit)

    while(meet.conditions){
        #compute all matrices - information about them found within their respective functions
        D <- computeD(mu=mu.vec)
        D_inv <- solve(D)
        y_star <- computey_star(X=X, curr_beta = curr_beta, full.Z = full.Z, D_inv = D_inv, curr_u = curr_u, y=y)
        V <- computeV(mu=mu.vec, r=dispersion)
        W <- computeW(D_inv=D_inv, V=V)
        W_inv <- solve(W)
        V_star <- computeV_star(full.Z=full.Z, curr_G=curr_G, W=W)
        tryCatch({V_star_inv <- solve(V_star, tol = 1e-150)}, error=function(cond){
                   message("V_star won't invert.")
                   return(final.list <- list("FE"=c(NA, NA),
                                    "RE"=NA,
                                    "Sigma"=NA,
                                    "Theta.Converged"=NA,
                                    "Sigma.Converged"=NA,
                                    "converged"=NA,
                                    "Iters"=NA,
                                    "Dispersion"=NA,
                                    "SE"=c(NA, NA),
                                    "Zscore"=c(NA, NA),
                                    "df" = c(NA, NA),
                                    "G"=NA,
                                    "VSTAR"=NA,
                                    "Hessian"=NA,
                                    "pvalue"=c(NA, NA)))
                 })
        V_partial <- computeV_partial(full.Z=full.Z, random.levels=random.levels, u_indices=u_indices)

        matrix.list <- preComputeMatrices(V_star_inv, V_partial, X, curr_beta, full.Z, curr_u, y_star)
        #---- First estimate variance components with Newton Raphson procedure ---#
        if (isFALSE(REML)) {
            P <- computeP_REML(V_star_inv=V_star_inv, X=X) #needed for computeVarCovar
            if (isTRUE(P)){
              return(final.list <- list("FE"=c(NA, NA),
                                        "RE"=NA,
                                        "Sigma"=NA,
                                        "Theta.Converged"=NA,
                                        "Sigma.Converged"=NA,
                                        "converged"=NA,
                                        "Iters"=NA,
                                        "Dispersion"=NA,
                                        "SE"=c(NA, NA),
                                        "Zscore"=c(NA, NA),
                                        "df" = c(NA, NA),
                                        "G"=NA,
                                        "VSTAR"=NA,
                                        "Hessian"=NA,
                                        "pvalue"=c(NA, NA)))
            }
            PV <- computePV(V_partial=V_partial, P=P) #needed for computeVarCovar
            score_sigma <- sigmaScore(matrix_list=matrix.list, V_star_inv=V_star_inv, V_partial=V_partial, random.levels=random.levels)
            information_sigma <- sigmaInformation(V_star_inv=V_star_inv, V_partial=V_partial, random.levels=random.levels)
        } else if (isTRUE(REML)) {
            P <- computeP_REML(V_star_inv=V_star_inv, X=X)
            if (isTRUE(P)){
              return(final.list <- list("FE"=c(NA, NA),
                                        "RE"=NA,
                                        "Sigma"=NA,
                                        "Theta.Converged"=NA,
                                        "Sigma.Converged"=NA,
                                        "converged"=NA,
                                        "Iters"=NA,
                                        "Dispersion"=NA,
                                        "SE"=c(NA, NA),
                                        "Zscore"=c(NA, NA),
                                        "df" = c(NA, NA),
                                        "G"=NA,
                                        "VSTAR"=NA,
                                        "Hessian"=NA,
                                        "pvalue"=c(NA, NA)))
            }
            PV <- computePV(V_partial=V_partial, P=P)
            score_sigma <- sigmaScoreREML(PV=PV, P=P, y_star=y_star, random.levels=random.levels)
            information_sigma <- sigmaInformationREML(PV=PV, random.levels=random.levels)
        }

        # update sigma, G, and G_inv
        curr_sigma <- sigma_update
        curr_G <- initialiseG(cluster_levels=random.levels, sigmas=curr_sigma)
        G_inv <- solve(curr_G)

        #---- Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u---####
        theta_update <- solve_equations(X=X, W_inv=W_inv, full.Z=full.Z, G_inv=G_inv, curr_beta=curr_beta, curr_u=curr_u, y_star=y_star)
        if (isTRUE(theta_update)) {
          warning("Hessian is computationally singular - cannot solve GLMM")
          return(final.list <- list("FE"=c(NA, NA),
                                    "RE"=NA,
                                    "Sigma"=NA,
                                    "Theta.Converged"=NA,
                                    "Sigma.Converged"=NA,
                                    "converged"=NA,
                                    "Iters"=NA,
                                    "Dispersion"=NA,
                                    "SE"=c(NA, NA),
                                    "Zscore"=c(NA, NA),
                                    "df" = c(NA, NA),
                                    "G"=NA,
                                    "VSTAR"=NA,
                                    "Hessian"=NA,
                                    "pvalue"=c(NA, NA)))
        }

        theta_diff <- abs(theta_update - curr_theta)

        # update B, u and mu_vec to determine new values of score and hessian matrices
        curr_theta <- theta_update
        rownames(curr_theta) <- c(colnames(X), colnames(full.Z))
        curr_beta <- curr_theta[colnames(X), , drop=FALSE]
        curr_u <- curr_theta[colnames(full.Z), , drop=FALSE]
        mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))

        if (any(is.infinite(mu.vec))) {
          warning("Estimates increasing to infinity - cannot solve GLMM.")
          return(final.list <- list("FE"=c(NA, NA),
                                    "RE"=NA,
                                    "Sigma"=NA,
                                    "Theta.Converged"=NA,
                                    "Sigma.Converged"=NA,
                                    "converged"=NA,
                                    "Iters"=NA,
                                    "Dispersion"=NA,
                                    "SE"=c(NA, NA),
                                    "Zscore"=c(NA, NA),
                                    "df" = c(NA, NA),
                                    "G"=NA,
                                    "VSTAR"=NA,
                                    "Hessian"=NA,
                                    "pvalue"=c(NA, NA)))
        }

        iters <- iters + 1
        meet.conditions <- !((all(theta_diff < theta.conv)) & (all((sigma_diff) < theta.conv))| iters >= max.hit)
    }

    SE <- calculateSE(X=X, full.Z=full.Z, W_inv=W_inv, G_inv=G_inv)
    Zscore <- calculateZscore(curr_beta=curr_beta, SE=SE)
    Va <- computeVarCovar(random.levels, PV)
    mint <- nrow(curr_beta)
    cint <- nrow(curr_u)
    coeff.matrix <- makeCoefMatrix(X=X, full.Z=full.Z, W_inv=W_inv, G_inv=G_inv)
    df <- Satterthwaite_df(coeff.mat=coeff.matrix, mint=mint, cint=cint, SE=SE, V_a=Va,
                           V_partial=V_partial, G_inv=G_inv, curr_sigma=curr_sigma, curr_beta=curr_beta, random.levels=random.levels)
    Pvalue <- computePvalue(Zscore=Zscore, df=df)

    converged <- ((all(theta_diff < theta.conv)) & (all(abs(sigma_diff) < theta.conv)))
    print(converged)
    if (isFALSE(converged)){
      message(paste("Model has not converged after", iters, "iterations. Consider increasing max.iter or drop a random effect."))
    }
    
    final.list <- list("FE"=as.vector(curr_beta),
                       "RE"=as.vector(curr_u),
                       "Sigma"=as.vector(curr_sigma),
                       "Theta.Converged"=theta_diff < theta.conv,
                       "Sigma.Converged"=sigma_diff < theta.conv,
                       "converged"=converged,
                       "Iters"=iters,
                       "Dispersion"=dispersion,
                       "SE"=SE,
                       "Zscore"=Zscore,
                       "df" = df,
                       "G"=curr_G,
                       "VSTAR"=V_star,
                       "Hessian"=information_sigma,
                       "pvalue"=Pvalue)

    return(final.list)
}


#' @importMethodsFrom Matrix %*%
#' @export
computeW <- function(D_inv, V){
    W = D_inv %*% V %*% D_inv
    return(W)
}


#' @importFrom Matrix Matrix
#' @export
computeV <- function(mu, r){
    # compute diagonal matrix of variances
    v.vec <- ((mu**2/r)) + mu
    V <- Matrix(0L, ncol=length(mu), nrow=length(mu), sparse = TRUE)
    diag(V) <- v.vec
    return(V)
}


#' @importFrom Matrix Matrix
#' @export
computeD <- function(mu){
    # D is diag(mu_i)
    D <- Matrix(0L, ncol=length(mu), nrow=length(mu), sparse = TRUE)
    diag(D) <- mu
    return(D)
}


#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix crossprod
#' @export
computeV_star <- function(full.Z, curr_G, W){
    # V_star_R <- (full.Z %*% curr_G %*% t(full.Z)) + W
    V_star_C <- computeVStar(as.matrix(full.Z), as.matrix(curr_G), as.matrix(W))
    return(V_star_C)
}
        
#' @importMethodsFrom Matrix %*%
#' @export
computeP_REML <- function(V_star_inv, X) {
    # breaking these down to individual steps speeds up the operations considerably
    tx.m <- t(X) %*% V_star_inv
    x.inv <- computeInv(tx.m %*% X)
    if (isTRUE(x.inv)) {
      return(P=TRUE)
    }
    vx <- V_star_inv %*% X
    Pminus <- vx %*% x.inv
    tx.inv <- t(X) %*% V_star_inv
    P <- V_star_inv - Pminus %*% tx.inv
    return(P)
}



#' @importFrom Matrix sparseMatrix diag
#' @export
initialiseG <- function(cluster_levels, sigmas, Kin=NULL){
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
        if(!is.null(Kin)){
            diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ] # is this sufficient to transform the sigma to the model scale?
        } else{
            if(rownames(sigmas[x, , drop=FALSE]) %in% c("Genetic")){
                diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ] * Kin
            }else{
                diag(G[c(i:(i+x.q-1)), c(i:(i+x.q-1)), drop=FALSE]) <- sigmas[x, ] # is this sufficient to transform the sigma to the model scale?
            }
        }

        i <- j <- i+x.q
    }
    return(as.matrix(G))
}


#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix diag
#' @export
initializeFullZ <- function(Z, cluster_levels, stand.cols=FALSE){
    # construct the full Z with all random effect levels

    # check that all of the levels are present in random.levels AND the
    # entries of Z
    all.present <- unlist(sapply(seq_along(cluster_levels), FUN=function(PX){
        all(cluster_levels[[PX]] %in% unique(Z[, PX]))
        }, simplify=FALSE))

    if(!all(all.present)){
        stop("Columns of Z are discordant with input random effect levels")
    }

    n.cols <- ncol(Z)
    z.names <- colnames(Z)
    if(is.null(z.names)){
        stop("Columns of Z must have valid names")
    }

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
            i.levels <- unique(Z[, i])
            i.levels <- as.factor(paste(sort(as.integer(i.levels))))
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
        }

        colnames(i.z) <- paste0(colnames(Z)[i], cluster_levels[[colnames(Z)[i]]])

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

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
computePV <- function(V_partial, P){
  PV <- list()
  for (i in 1:length(V_partial)) {
    PV[[i]] <- P %*% V_partial[[i]]
  }
  return(PV)
}

## @importFrom Matrix Matrix
## @export
#initializeFullZ <- function(Z) {
#  full.Z <- matrix(0L, nrow=nrow(Z), ncol = 0)
#  for (i in 1:ncol(Z)) {
#    temp.Z <- Matrix(table(seq_along(1:nrow(Z)), Z[,i]), sparse = TRUE)
#    full.Z <- cbind(full.Z, temp.Z)
#  }
#  return(full.Z)
#}

#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix solve
#' @export
solve_equations <- function(X, W_inv, full.Z, G_inv, curr_beta, curr_u, y_star){
  
    UpperLeft <- t(X) %*% W_inv %*% X
    UpperRight <- t(X) %*% W_inv %*% full.Z
    LowerLeft <- t(full.Z) %*% W_inv %*% X
    LowerRight <- t(full.Z) %*% W_inv %*% full.Z + G_inv

    LHS <- rbind(cbind(UpperLeft, UpperRight), cbind(LowerLeft, LowerRight))
    RHS <- rbind((t(X) %*% W_inv %*% y_star), (t(full.Z) %*% W_inv %*% y_star))

    theta_update <- tryCatch({solve(LHS) %*% RHS}
             , error = function(e){
               exit <- TRUE
               }
    )
    return(theta_update)
}


#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix
#' @export
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


### utility functions
#' @importFrom Matrix solve
#' @export
computeInv <- function(x){
    # Compute x^-1 from x
    # need to check that x is not singular - use tryCatch - if matrix is singular then report error message
    x_inv <- tryCatch(expr={
        solve(x)
    },
    error=function(cond){
        message("Matrix cannot be inverted - most likely singular")
        message(cond)
        exit <- TRUE
    })
    return(x_inv)
}

#' Wrapper function for fitPLGlmm
#'
#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix solve crossprod
#' @export
fitGLMM <- function(X, Z, y, init.theta=NULL, crossed=FALSE, random.levels=NULL, REML=FALSE,
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
    curr_sigma <- matrix(unlist(lapply(mapUtoIndiv(full.Z, curr_u, random.levels=random.levels),
                                       FUN=function(Bj){
                                           (1/(length(Bj)-1)) * crossprod(Bj, Bj)
                                           })), ncol=1)
    rownames(curr_sigma) <- colnames(Z)

    #create a single variable for the thetas
    curr_theta <- do.call(rbind, list(curr_beta, curr_u))

    #compute mu.vec using inverse link function
    mu.vec <- exp((X %*% curr_beta) + (full.Z %*% curr_u))

    # use y_bar as the sample mean and s_hat as the sample variance
    new.r <- dispersion

    # theta_diff <- rep(Inf, nrow(curr_theta))
    # sigma_diff <- rep(80000, nrow(curr_sigma))

    #compute variance-covariance matrix G
    curr_G <- initialiseG(cluster_levels=random.levels, sigmas=curr_sigma)

    conv.list <- list()
    iters <- 1

    u_indices <- sapply(seq_along(random.levels),
                        FUN=function(RX) which(colnames(full.Z) %in% random.levels[[RX]]),
                        simplify=FALSE)

    # flatten column matrices to vectors
    mu.vec <- mu.vec[, 1]
    curr_beta <- curr_beta[, 1]
    curr_u <- curr_u[, 1]
    curr_theta <- curr_theta[, 1]
    curr_sigma <- curr_sigma[, 1]

    final.list <- fitPLGlmm(Z=full.Z, X=X, muvec=mu.vec, curr_beta=curr_beta,
                            curr_theta=curr_theta, curr_u=curr_u, curr_sigma=curr_sigma,
                            curr_G=curr_G, y=y, u_indices=u_indices, theta_conv=theta.conv, rlevels=random.levels,
                            curr_disp=new.r, REML=TRUE, maxit=15)

    # compute Z scores, DF and P-values
    mint <- length(curr_beta)
    cint <- length(curr_u)
    dfs <- Satterthwaite_df(final.list[["COEFF"]], mint, cint, final.list[["SE"]], final.list[["Sigma"]], final.list[["FE"]],
                            final.list[["Vpartial"]], final.list[["VCOV"]], final.list[["Ginv"]], random.levels)
    pvals <- computePvalue(final.list[["t"]], dfs)

    final.list[["DF"]] <- dfs
    final.list[["PVALS"]] <- pvals

    return(final.list)
}

#' @importFrom Matrix diag
#' @export
matrix.trace <- function(x){
    # check is square matrix first
    x.dims <- dim(x)
    if(x.dims[1] != x.dims[2]){
        stop("matrix is not square")
    } else{
        return(sum(diag(x)))
    }
}


#' @importFrom MASS ginv
#' @export
inv <- function(x){
  if (det(x) > 8.313969e-20) {
    solve(x)
  } else {
    stop("det is 0")
  }
}


#' @export
glmmControl.defaults <- function(...){
    # return the default glmm control values
    return(list(theta.tol=1e-6, max.iter=100))
}


#' @importFrom stats pt
#' @export
computePvalue <- function(Zscore, df) {
    pval <- 2*pt(abs(Zscore), df, lower.tail=FALSE)
    return(pval)
}


#' @importMethodsFrom Matrix %*% t
#' @importFrom Matrix solve diag
###---- first calculate g = derivative of C with respect to sigma ----
function_jac <- function(x, coeff.mat, mint, cint, G_inv, random.levels) {
    UpperLeft <- coeff.mat[c(1:mint), c(1:mint)]
    UpperRight <- coeff.mat[c(1:mint), c((mint+1):(mint+cint))]
    LowerLeft <- coeff.mat[c((mint+1):(mint+cint)), c(1:mint)]
    LowerRight <- coeff.mat[c((mint+1):(mint+cint)), c((mint+1):(mint+cint))] - G_inv

    n <- length(random.levels)
    diag(LowerRight) <- diag(LowerRight) + rep(1/x, times=lengths(random.levels)) #when extending to random slopes, this needs to be changed to a matrix and added to LowerRight directly
    C <- solve(UpperLeft - UpperRight %*% solve(LowerRight) %*% LowerLeft)
}


#' @importMethodsFrom Matrix %*% t
#' @importFrom Matrix solve diag
#' @importFrom numDeriv jacobian
#' @export
Satterthwaite_df <- function(coeff.mat, mint, cint, SE, curr_sigma, curr_beta, V_partial, V_a, G_inv, random.levels) {

  if(any(class(curr_sigma) %in% c("Matrix", "matrix", "dgeMatrix", "dgCMatrix"))){
    curr_sigma <- as.vector(curr_sigma)
  } else{
      stop("Class of curr_sigma not recognised")
  }

  jac <- jacobian(func=function_jac, x=curr_sigma, coeff.mat=coeff.mat, mint=mint, cint=cint, G_inv=G_inv, random.levels=random.levels)
  jac_list <- lapply(1:ncol(jac), function(i)
    array(jac[, i], dim=rep(length(curr_beta), 2))) #when extending to random slopes, this would have to be reformatted into list, where each element belongs to one random effect

  # V_a is provided externally
  df <- rep(NA, length(curr_beta))
  for (i in 1:length(curr_beta)) {
    jac_var_beta <- matrix(unlist(lapply(lapply(jac_list, diag), `[[`, i)), ncol=1) # could this be done with AD?
    denom <- t(jac_var_beta) %*% (V_a) %*% jac_var_beta #g' Va g
    df[i] <- 2*((SE[i]^2)^2)/denom
  }
  return(as.matrix(df))
}
