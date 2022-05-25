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

    if(isFALSE(geno.only) & !is.null(Kin)){
        # Kin must be nXn
        if(nrow(Kin) != ncol(Kin)){
            stop("Input covariance matrix is not square: ", nrow(Kin), "x", ncol(Kin))
        }

        if(nrow(Kin) != nrow(Z)){
            stop("Input covariance matrix and Z design matrix are discordant: ",
                 nrow(Z), "x", ncol(Z), ", ", nrow(Kin), "x", ncol(Kin))
        }

        # create full Z with expanded random effect levels
        full.Z <- initializeFullZ(Z=Z, cluster_levels=random.levels)

        # random value initiation from runif
        curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        rownames(curr_u) <- colnames(full.Z)

        # compute sample variances of the us
        curr_sigma <- Matrix(runif(ncol(Z), 0, 1), ncol=1, sparse = TRUE)
        rownames(curr_sigma) <- colnames(Z)

        ## add the genetic components
        ## augment Z with I
        geno.I <- diag(nrow(full.Z))
        colnames(geno.I) <- paste0("CovarMat", seq_len(ncol(geno.I)))
        full.Z <- do.call(cbind, list(full.Z, geno.I))

        # add a genetic variance component
        sigma_g <- Matrix(runif(1, 0, 1), ncol=1, nrow=1, sparse=TRUE)
        rownames(sigma_g) <- "CovarMat"
        curr_sigma <- do.call(rbind, list(curr_sigma, sigma_g))

        # add genetic BLUPs
        g_u <- matrix(runif(nrow(full.Z), 0, 1), ncol=1)
        rownames(g_u) <- colnames(geno.I)
        curr_u <- do.call(rbind, list(curr_u, g_u))

        random.levels <- c(random.levels, list("CovarMat"=colnames(geno.I)))

        #compute variance-covariance matrix G
        curr_G <- initialiseG(cluster_levels=random.levels, sigmas=curr_sigma, Kin=Kin)
    } else if(isTRUE(geno.only) & !is.null(Kin)){
        # Kin must be nXn
        if(nrow(Kin) != ncol(Kin)){
            stop("Input covariance matrix is not square: ", nrow(Kin), "x", ncol(Kin))
        }

        # this means Z is an identity matrix - need to check this
        if(ncol(Z) != nrow(Z)){
            stop("Z matrix is not square: ", nrow(Z), "x", ncol(Z))
        }

        full.Z <- Z
        colnames(full.Z) <- paste0(names(random.levels), seq_len(ncol(full.Z)))

        # random value initiation from runif
        curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        rownames(curr_u) <- colnames(full.Z)

        # compute sample variances of the us
        curr_sigma <- Matrix(runif(1, 0, 1), ncol=1, sparse = TRUE)
        rownames(curr_sigma) <- names(random.levels)

        #compute variance-covariance matrix G
        curr_G <- initialiseG(cluster_levels=random.levels, sigmas=curr_sigma, Kin=Kin)
    } else if(is.null(Kin)){
        # create full Z with expanded random effect levels
        full.Z <- initializeFullZ(Z=Z, cluster_levels=random.levels)

        # random value initiation from runif
        curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        rownames(curr_u) <- colnames(full.Z)

        # compute sample variances of the us
        curr_sigma <- Matrix(runif(ncol(Z), 0, 1), ncol=1, sparse = TRUE)
        rownames(curr_sigma) <- colnames(Z)

        #compute variance-covariance matrix G
        curr_G <- initialiseG(cluster_levels=random.levels, sigmas=curr_sigma)
    }

    # create a single variable for the thetas
    curr_theta <- do.call(rbind, list(curr_beta, curr_u))
    if(nrow(curr_theta) != sum(c(ncol(X), ncol(full.Z)))){
        stop("Number of parameters does not match columns of X and Z")
    }

    #compute mu.vec using inverse link function
    mu.vec <- exp(offsets + (X %*% curr_beta) + (full.Z %*% curr_u))
    if(any(is.infinite(mu.vec))){
        stop("Infinite values in initial estimates - reconsider model")
    }

    if(isTRUE(any(is.na(mu.vec[, 1])))){
        if(isTRUE(any(is.na(offsets)))){
            stop("NA values in offsets - remove these samples before re-running model")
        } else{
            stop("NAs values in initial estimates - remove these samples before re-running model")
        }
    }

    # be careful here as the colnames of full.Z might match multiple RE levels <- big source of bugs!!!
    u_indices <- sapply(seq_along(names(random.levels)),
                        FUN=function(RX) {
                            which(colnames(full.Z) %in% random.levels[[RX]])
                        }, simplify=FALSE)

    if(sum(unlist(lapply(u_indices, length))) != ncol(full.Z)){
        stop("Non-unique column names in Z - please ensure these are unique")
    }

    # flatten column matrices to vectors
    mu.vec <- mu.vec[, 1]
    curr_beta <- curr_beta[, 1]

    curr_sigma <- curr_sigma[, 1]
    curr_u <- curr_u[, 1]
    curr_theta <- curr_theta[, 1]

    if(is.null(Kin)){
        final.list <- fitPLGlmm(Z=full.Z, X=X, muvec=mu.vec, offsets=offsets, curr_beta=curr_beta,
                                curr_theta=curr_theta, curr_u=curr_u, curr_sigma=curr_sigma,
                                curr_G=as.matrix(curr_G), y=y, u_indices=u_indices, theta_conv=theta.conv, rlevels=random.levels,
                                curr_disp=dispersion, REML=TRUE, maxit=max.hit)
    } else{
        final.list <- fitGeneticPLGlmm(Z=full.Z, X=X, K=Kin, offsets=offsets,
                                       muvec=mu.vec, curr_beta=curr_beta,
                                       curr_theta=curr_theta, curr_u=curr_u, curr_sigma=curr_sigma,
                                       curr_G=curr_G, y=y, u_indices=u_indices, theta_conv=theta.conv, rlevels=random.levels,
                                       curr_disp=dispersion, REML=TRUE, maxit=max.hit)
    }
<<<<<<< HEAD

    if(isFALSE(final.list$converged)){
        warning("Model has not converged after ", final.list$Iters,
                " iterations. Consider increasing max.iter or drop a random effect")
    }

=======
    
    # if(isFALSE(final.list$converged)){
    #     warning("Model has not converged after ", final.list$Iters,
    #             " iterations. Consider increasing max.iter or drop a random effect")
    # }
    
>>>>>>> 9f5431d (resolving conflicts)
    # compute Z scores, DF and P-values
    mint <- length(curr_beta)
    cint <- length(curr_u)
    dfs <- Satterthwaite_df(final.list[["COEFF"]], mint, cint, final.list[["SE"]], final.list[["Sigma"]], final.list[["FE"]],
                            final.list[["Vpartial"]], final.list[["VCOV"]], final.list[["Ginv"]], random.levels)
    pvals <- computePvalue(final.list[["t"]], dfs)

    if(any(is.infinite(pvals))){
        stop("Setting infinite p-values to NA")
        pvals[is.infinite(pvals)] <- NA
    }

    final.list[["DF"]] <- dfs
    final.list[["PVALS"]] <- pvals

    # final checks
    na.params <- is.na(c(final.list[["Sigma"]], final.list[["FE"]], final.list[["RE"]]))
    if(sum(na.params) > 0){
        stop("NA parameter estimates - reconsider model")
    }

    inf.params <- is.infinite(c(final.list[["Sigma"]], final.list[["FE"]], final.list[["RE"]]))
    if(sum(inf.params) > 0){
        stop("Infinite parameter estimates - reconsider model or increase sample size")
    }

    return(final.list)
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
        if(is.numeric(unique(Z[, PX]))){
            all(cluster_levels[[PX]] %in% paste0(names(cluster_levels)[PX], unique(Z[, PX])))
        } else{
            all(cluster_levels[[PX]] %in% unique(Z[, PX]))
        }
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
