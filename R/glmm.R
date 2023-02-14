#' Perform differential abundance testing using a NB-generalised linear mixed model
#'
#' This function will perform DA testing per-nhood using a negative binomial generalised linear mixed model
#' @param X A matrix containing the fixed effects of the model.
#' @param Z A matrix containing the random effects of the model.
#' @param y A matrix containing the observed phenotype over each neighborhood.
#' @param offsets A vector containing the (log) offsets to apply normalisation for different numbers of cells across samples.
#' @param init.theta A column vector (m X 1 matrix) of initial estimates of fixed and random effect coefficients
#' @param Kin A n x n covariance matrix to explicitly model variation between observations
#' @param REML A logical value denoting whether REML (Restricted Maximum Likelihood) should be run. Default is TRUE.
#' @param random.levels A list describing the random effects of the model, and for each, the different unique levels.
#' @param glmm.control A list containing parameter values specifying the theta tolerance of the model, the maximum number of iterations to be run,
#' initial parameter values for the fixed (init.beta) and random effects (init.u), and glmm solver (see details).
#' @param dispersion A scalar value for the dispersion of the negative binomial.
#' @param geno.only A logical value that flags the model to use either just the \code{matrix} `Kin` or the supplied random effects.
#' @param solver a character value that determines which optmisation algorithm is used for the variance components. Must be either
#' HE (Haseman-Elston regression) or Fisher (Fisher scoring).
#' @param var.dist A character value that determines the form of the variance - either Negative Binomial (NB) or Poisson (P). The latter is used
#' to model overdispersion using random effects alone, and therefore ignores the \code{dispersion} parameter. To avoid warnings while using the
#' Poisson variance form set \code{dispersion = NA}.
#'
#' @details
#' This function runs eithr a negative binomial of Poisson generalised linear mixed effects model. If mixed effects are detected in testNhoods,
#' this function is run to solve the model. The solver defaults to the \emph{Fisher} optimiser, and in the case of negative variance estimates
#' it will switch to the non-negative least squares (NNLS) Haseman-Elston solver. This behaviour can be pre-set by passing
#' \code{glmm.control$solver="HE"} for Haseman-Elston regression, which is the recommended solver when a covariance matrix is provided,
#' or \code{glmm.control$solver="HE-NNLS"} which is the constrained HE optimisation algorithm. \code{var.dist} controls the form of the variance function
#' \(either NB or P\). For the NB-GLMM case the dispersion parameter handles some of the counts overdispersion. However, the Poisson model
#' instead includes a "residual" variance parameter which is also reported as \emph{Resid.Sigma}.
#'
#' @return  A list containing the GLMM output, including inference results. The list elements are as follows:
#' \describe{
#' \item{\code{FE}:}{\code{numeric} vector of fixed effect parameter estimates.}
#' \item{\code{RE}:}{\code{list} of the same length as the number of random effect variables. Each slot contains the best
#' linear unbiased predictors (BLUPs) for the levels of the corresponding RE variable.}
#' \item{\code{Sigma:}}{\code{numeric} vector of variance component estimates, 1 per random effect variable.}
#' \item{\code{converged:}}{\code{logical} scalar of whether the model has reached the convergence tolerance or not.}
#' \item{\code{Iters:}}{\code{numeric} scalar with the number of iterations that the model ran for. Is strictly <= \code{max.iter}.}
#' \item{\code{Dispersion:}}{\code{numeric} scalar of the dispersion estimate computed off-line}
#' \item{\code{Hessian:}}{\code{matrix} of 2nd derivative elements from the fixed and random effect parameter inference.}
#' \item{\code{SE:}}{\code{matrix} of standard error estimates, derived from the hessian, i.e. the square roots of the diagonal elements.}
#' \item{\code{t:}}{\code{numeric} vector containing the compute t-score for each fixed effect variable.}
#' \item{\code{COEFF:}}{\code{matrix} containing the coefficient matrix from the mixed model equations.}
#' \item{\code{P:}}{\code{matrix} containing the elements of the REML projection matrix.}
#' \item{\code{Vpartial:}}{\code{list} containing the partial derivatives of the (pseudo)variance matrix with respect to each variance
#' component.}
#' \item{\code{Ginv:}}{\code{matrix} of the inverse variance components broadcast to the full Z matrix.}
#' \item{\code{Vsinv:}}{\code{matrix} of the inverse pseudovariance.}
#' \item{\code{Winv:}}{\code{matrix} of the inverse elements of W = D^-1 V D^-1}
#' \item{\code{VCOV:}}{\code{matrix} of the variance-covariance for all model fixed and random effect variable parameter estimates.
#' This is required to compute the degrees of freedom for the fixed effect parameter inference.}
#' \item{\code{DF:}}{\code{numeric} vector of the number of inferred degrees of freedom. For details see \link{Satterthwaite_df}.}
#' \item{\code{PVALS:}}{\code{numeric} vector of the compute p-values from a t-distribution with the inferred number of degrees of
#' freedom.}
#' \item{\code{ERROR:}}{\code{list} containing Rcpp error messages - used for internal checking.}
#' }
#' @author Mike Morgan
#'
#' @examples
#' data(sim_nbglmm)
#' random.levels <- list("RE1"=paste("RE1", levels(as.factor(sim_nbglmm$RE1)), sep="_"),
#'                       "RE2"=paste("RE2", levels(as.factor(sim_nbglmm$RE2)), sep="_"))
#' X <- as.matrix(data.frame("Intercept"=rep(1, nrow(sim_nbglmm)), "FE2"=as.numeric(sim_nbglmm$FE2)))
#' Z <- as.matrix(data.frame("RE1"=paste("RE1", as.numeric(sim_nbglmm$RE1), sep="_"),
#'                           "RE2"=paste("RE2", as.numeric(sim_nbglmm$RE2), sep="_")))
#' y <- sim_nbglmm$Mean.Count
#' dispersion <- 0.5
#'
#' glmm.control <- glmmControl.defaults()
#' glmm.control$theta.tol <- 1e-6
#' glmm.control$max.iter <- 15
#' model.list <- fitGLMM(X=X, Z=Z, y=y, offsets=rep(0, nrow(X)), random.levels=random.levels,
#'                       REML = TRUE, glmm.control=glmm.control, dispersion=dispersion, solver="Fisher")
#' model.list
#'
#' @name fitGLMM
#'
#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix solve crossprod
#' @importFrom stats runif
#' @export
fitGLMM <- function(X, Z, y, offsets, init.theta=NULL, Kin=NULL,
                    random.levels=NULL, REML=FALSE,
                    glmm.control=list(theta.tol=1e-6, max.iter=100,
                                      init.sigma=NULL, init.beta=NULL,
                                      init.u=NULL, solver=NULL),
                    dispersion = 0.5, geno.only=FALSE,
                    solver=NULL, var.dist="NB"){

    if(!glmm.control$solver %in% c("HE", "Fisher", "HE-NNLS")){
        stop(glmm.control$solver, " not recognised - must be HE, HE-NNLS or Fisher")
    }

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
    if(is.null(glmm.control[["init.beta"]])){
        curr_beta <- solve((t(X) %*% X)) %*% t(X) %*% log(y + 1)
    } else{
        curr_beta = matrix(glmm.control[["init.beta"]], ncol=1)
    }
    rownames(curr_beta) <- colnames(X)

    if(!var.dist %in% c("NB", "P")){
        stop("var.dist value ", var.dist, " not recognised - must be either NB or P")
    }

    if(var.dist %in% c("P") & !is.na(dispersion)){
        warning("Using dispersion with Poisson variance - dispersion will be ignored")
    }

    if(var.dist %in% c("P") & is.null(Kin)){
        # this only works here when no kinship - need to add later if kinship
        message("Adding residual variance component variable")
        eye <- matrix(seq_len(nrow(Z)), ncol=1, nrow=nrow(Z))
        colnames(eye) <- c("Resid")
        Z <- cbind(Z, eye)
        random.levels[["Resid"]] <- paste0("Resid", seq_len(nrow(Z)))
    }

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

        # add the I matrix for P var form and kinship matrix
        if(var.dist %in% c("P")){
            eye <- matrix(0L, nrow=nrow(Z), ncol=nrow(Z))
            colnames(eye) <- paste0("Resid", seq_len(ncol(eye)))
            full.Z <- cbind(full.Z, eye)
            random.levels[["Resid"]] <- colnames(eye)
        }

        # random value initiation from runif
        if(is.null(glmm.control[["init.u"]])){
            curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        } else{
            curr_u <- matrix(glmm.control[["init.u"]], ncol=1)
        }
        rownames(curr_u) <- colnames(full.Z)

        # compute sample variances of the us
        if(is.null(glmm.control[["init.sigma"]])){
            curr_sigma <- Matrix(runif(ncol(Z), 0, 1), ncol=1, sparse = TRUE)
        } else{
            curr_sigma <- Matrix(glmm.control[["init.sigma"]], ncol=1, sparse=TRUE)
        }
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

        # if we only have a GRM then Z _is_ full.Z?
        # full.Z <- initializeFullZ(Z, cluster_levels=random.levels)
        full.Z <- Z
        # should this be the matrix R?
        colnames(full.Z) <- paste0(names(random.levels), seq_len(ncol(full.Z)))

        # add the I matrix for P var form and kinship matrix
        if(var.dist %in% c("P")){
            eye <- matrix(0L, nrow=nrow(Z), ncol=nrow(Z))
            colnames(eye) <- paste0("Resid", seq_len(ncol(eye)))
            full.Z <- cbind(full.Z, eye)
            random.levels[["Resid"]] <- colnames(eye)
        }

        # random value initiation from runif
        if(is.null(glmm.control[["init.u"]])){
            curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        } else{
            curr_u <- matrix(glmm.control[["init.u"]], ncol=1)
        }
        rownames(curr_u) <- colnames(full.Z)

        # compute sample variances of the us
        if(is.null(glmm.control[["init.sigma"]])){
            curr_sigma <- Matrix(runif(length(random.levels), 0, 1), ncol=1, sparse = TRUE)
        } else{
            curr_sigma <- Matrix(glmm.control[["init.sigma"]], ncol=1, sparse=TRUE)
        }

        rownames(curr_sigma) <- names(random.levels)

        #compute variance-covariance matrix G
        curr_G <- initialiseG(cluster_levels=random.levels, sigmas=curr_sigma, Kin=Kin)
    } else if(is.null(Kin)){
        # create full Z with expanded random effect levels
        full.Z <- initializeFullZ(Z=Z, cluster_levels=random.levels)

        # random value initiation from runif
        if(is.null(glmm.control[["init.u"]])){
            curr_u <- matrix(runif(ncol(full.Z), 0, 1), ncol=1)
        } else{
            curr_u <- matrix(glmm.control[["init.u"]], ncol=1)
        }
        rownames(curr_u) <- colnames(full.Z)

        # compute sample variances of the us
        if(is.null(glmm.control[["init.sigma"]])){
            curr_sigma <- Matrix(runif(ncol(Z), 0, 1), ncol=1, sparse = TRUE)
        } else{
            curr_sigma <- Matrix(glmm.control[["init.sigma"]], ncol=1, sparse=TRUE)
        }
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
        final.list <- tryCatch(fitPLGlmm(Z=full.Z, X=X, muvec=mu.vec, offsets=offsets, curr_beta=curr_beta,
                                         curr_theta=curr_theta, curr_u=curr_u, curr_sigma=curr_sigma,
                                         curr_G=as.matrix(curr_G), y=y, u_indices=u_indices, theta_conv=theta.conv, rlevels=random.levels,
                                         curr_disp=dispersion, REML=TRUE, maxit=max.hit, solver=glmm.control$solver, vardist=var.dist),
                               error=function(err){
                                   message(err)
                                   return(list("FE"=NA, "RE"=NA, "Sigma"=NA,
                                               "converged"=NA, "Iters"=NA, "Dispersion"=NA,
                                               "Hessian"=NA, "SE"=NA, "t"=NA, "PSVAR"=NA,
                                               "COEFF"=NA, "P"=NA, "Vpartial"=NA, "Ginv"=NA,
                                               "Vsinv"=NA, "Winv"=NA, "VCOV"=NA, "DF"=NA, "PVALS"=NA,
                                               "ERROR"=err))
                                   })
    } else{
        final.list <- tryCatch(fitGeneticPLGlmm(Z=full.Z, X=X, K=as.matrix(Kin), offsets=offsets,
                                                muvec=mu.vec, curr_beta=curr_beta,
                                                curr_theta=curr_theta, curr_u=curr_u, curr_sigma=curr_sigma,
                                                curr_G=curr_G, y=y, u_indices=u_indices, theta_conv=theta.conv, rlevels=random.levels,
                                                curr_disp=dispersion, REML=TRUE, maxit=max.hit, solver=glmm.control$solver, vardist=var.dist),
                               error=function(err){
                                   message(err)
                                   return(list("FE"=NA, "RE"=NA, "Sigma"=NA,
                                               "converged"=NA, "Iters"=NA, "Dispersion"=NA,
                                               "Hessian"=NA, "SE"=NA, "t"=NA, "PSVAR"=NA,
                                               "COEFF"=NA, "P"=NA, "Vpartial"=NA, "Ginv"=NA,
                                               "Vsinv"=NA, "Winv"=NA, "VCOV"=NA, "DF"=NA, "PVALS"=NA,
                                               "ERROR"=err))
                                   })
    }

    check_na_output <- any(is.na(unlist(lapply(names(final.list)[!grepl(names(final.list), pattern="ERROR")],
                                               function(CH) final.list[[CH]]))))
    if(check_na_output){
        traceback(3)
        error_cat <- paste(unlist(final.list[["ERROR"]]), sep=", ")
        stop("No GLMM results returned - see traceback for errors. ", error_cat)
    }

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


#' Construct the initial G matrix
#'
#' This function maps the variance estimates onto the full \code{c x q} levels for each random effect. This
#' ensures that the matrices commute in the NB-GLMM solver. This function is included for reference, and
#' should not be used directly
#' @param cluster_levels A \code{list} containing the random effect levels for each variable
#' @param sigmas A \code{matrix} of c X 1, i.e. a column vector, containing the variance component estimates
#' @param Kin A \code{matrix} containing a user-supplied covariance matrix
#'
#' @details Broadcast the variance component estimates to the full \code{c\*q x c\*q} matrix.
#'
#' @return \code{matrix} of the full broadcast variance component estimates.
#' @author Mike Morgan & Alice Kluzer
#'
#' @examples
#' data(sim_nbglmm)
#' random.levels <- list("RE1"=paste("RE1", levels(as.factor(sim_nbglmm$RE1)), sep="_"),
#'                       "RE2"=paste("RE2", levels(as.factor(sim_nbglmm$RE2)), sep="_"))
#' rand.sigma <- matrix(runif(2), ncol=1)
#' rownames(rand.sigma) <- names(random.levels)
#' big.G <- initialiseG(random.levels, rand.sigma)
#' dim(big.G)
#'
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


#' Construct the full Z matrix
#'
#' Using a simplified version of the \code{n x c} Z matrix, with one column per variable, construct the fully broadcast
#' \code{n x (c*q)} binary matrix that maps each individual onto the random effect variable levels. It is not intended
#' for this function to be called by the user directly, but it can be useful to debug mappings between random effect
#' levels and input variables.
#' @param Z A \code{n x c} matrix containing the numeric or character levels
#' @param cluster_levels A \code{list} that maps the column names of Z onto the individual levels
#' @param stand.cols A logical scalar that determines if Z* should be computed which is the row-centered and
#' scaled version of the full Z matrix
#'
#' @details
#' To make sure that matrices commute it is necessary to construct the full \code{n x c*q} matrix. This is a binary
#' matrix where each level of each random effect occupies a column, and the samples/observations are mapped onto
#' the correct levels based on the input Z.
#'
#' @return \code{matrix} Fully broadcast Z matrix with one column per random effect level for all random effect variables
#' in the model.
#' @author Mike Morgan & Alice Kluzer
#'
#' @examples
#' data(sim_nbglmm)
#' random.levels <- list("RE1"=paste("RE1", levels(as.factor(sim_nbglmm$RE1)), sep="_"),
#'                       "RE2"=paste("RE2", levels(as.factor(sim_nbglmm$RE2)), sep="_"))
#' Z <- as.matrix(data.frame("RE1"=paste("RE1", as.numeric(sim_nbglmm$RE1), sep="_"),
#'                           "RE2"=paste("RE2", as.numeric(sim_nbglmm$RE2), sep="_")))
#' fullZ <- initializeFullZ(Z, random.levels)
#' dim(Z)
#' dim(fullZ)
#'
#' @importMethodsFrom Matrix %*%
#' @importFrom Matrix Matrix diag
#' @export
initializeFullZ <- function(Z, cluster_levels, stand.cols=FALSE){

    # construct the full Z with all random effect levels
    n.cols <- ncol(Z)
    z.names <- colnames(Z)
    if(is.null(z.names)){
        stop("Columns of Z must have valid names")
    }

    # check that all of the levels are present in random.levels AND the
    # entries of Z
    all.present <- unlist(sapply(seq_along(cluster_levels), FUN=function(PX){
        if(is.numeric(unique(Z[, PX])) & ncol(Z) != nrow(Z)){
            all(cluster_levels[[PX]] %in% paste0(names(cluster_levels)[PX], unique(Z[, PX])))
        } else if(ncol(Z) == nrow(Z)){
            all(cluster_levels[[PX]] %in% colnames(Z))
        } else {
            all(cluster_levels[[PX]] %in% unique(Z[, PX]))
            }
    }, simplify=FALSE))


    if(!all(all.present)){
        stop("Columns of Z are discordant with input random effect levels")
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


#' Compute the trace of a matrix
#'
#' Exactly what it says on the tin - compute the sum of the matrix diagonal
#' @param x A \code{matrix}
#'
#' @details It computes the matrix trace of a square matrix.
#'
#' @return \code{numeric} scalar of the matrix trace.
#' @author Mike Morgan
#'
#' @examples
#' matrix.trace(matrix(runif(9), ncol=3, nrow=3))
#'
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

#' glmm control default values
#'
#' This will give the default values for the GLMM solver
#' @param ... see \code{fitGLMM} for details
#
#' @details The default values for the parameter estimation convergence is 1e-6, and the
#' maximum number of iterations is 100. In practise if the solver converges it generally does
#' so fairly quickly on moderately well conditioned problems. The default solver is Fisher
#' scoring, but this will switch (with a warning produced) to the NNLS Haseman-Elston solver
#' if negative variance estimates are found.
#'
#' @return \code{list} containing the default values GLMM solver. This can be saved in the
#' user environment and then passed to \link{testNhoods} directly to modify the convergence
#' criteria of the solver that is used.
#' \describe{
#' \item{\code{theta.tol:}}{\code{numeric} scalar that sets the convergence threshold for the
#' parameter inference - this is applied globally to fixed and random effect parameters, and
#' to the variance estimates.}
#' \item{\code{max.iter:}}{\code{numeric} scalar that sets the maximum number of iterations that
#' the NB-GLMM will run for.}
#' \item{\code{solver:}}{\code{character} scalar that sets the solver to use. Valid values are
#' \emph{Fisher}, \emph{HE} or \emph{HE-NNLS}. See \link{fitGLMM} for details.}
#' }
#' @author Mike Morgan
#' @examples
#' mmcontrol <- glmmControl.defaults()
#' mmcontrol
#' mmcontrol$solver <- "HE-NNLS"
#' mmcontrol
#'
#' @export
glmmControl.defaults <- function(...){
    # return the default glmm control values
    return(list(theta.tol=1e-6, max.iter=100, solver='Fisher'))
}


#' Compute the p-value for the fixed effect parameters
#'
#' Based on the asymptotic t-distribution, comptue the 2-tailed p-value that estimate != 0. This
#' function is not intended to be used directly, but is included for reference or if an alternative
#' estimate of the degrees of freedom is available.
#' @param Zscore A numeric vector containing the Z scores for each fixed effect parameter
#' @param df A numeric vector containing the estimated degrees of freedom for each fixed effect
#' parameter
#'
#' @details Based on sampling from a 2-tailed t-distribution with \code{df} degrees of freedom,
#' compute the probability that the calculated \code{Zscore} is greater than or equal to what would be
#' expected from random chance.
#' @return Numeric vector of p-values, 1 per fixed effect parameter
#' @author Mike Morgan & Alice Kluzer
#' @examples
#' NULL
#'
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


#' Compute degrees of freedom using Satterthwaite method
#'
#' This function is not intended to be called by the user, and is included for reference
#' @param coeff.mat A \code{matrix} class object containing the coefficient matrix from the mixed model equations
#' @param mint A numeric scalar of the number of fixed effect variables in the model
#' @param cint A numeric scalar of the number of random effect variables in the model
#' @param SE A \code{1 x mint} \code{matrix}, i.e. column vector, containing the standard errors of the fixed effect
#' parameter estimates
#' @param curr_sigma A \code{1 x cint matrix}, i.e. column vector, of the variance component parameter estimates
#' @param curr_beta A \code{1 x mint matrix}, i.e. column vector, of the fixed effect parameter estimates
#' @param V_partial A \code{list} of the partial derivatives for each fixed and random effect variable in the model
#' @param V_a A \code{c+m x c+m} variance-covariance matrix of the fixed and random effect variable parameter estimates
#' @param G_inv A \code{nxc X nxc} inverse matrix containing the variance component estimates
#' @param random.levels A \code{list} containing the mapping between the random effect variables and each respective set
#' of levels for said variable.
#'
#' @details The Satterthwaite degrees of freedom are computed, which estimates the numbers of degrees of freedom in the
#' NB-GLMM based on ratio of the squared standard errors and the product of the Jacobians of the variance-covariance matrix
#' from the fixed effect variable parameter estimation with full variance-covariance matrix. For more details see
#' Satterthwaite FE, Biometrics Bulletin (1946) Vol 2 No 6, pp110-114.
#'
#' @return \code{matrix} containing the inferred number of degrees of freedom for the specific model.
#' @author Mike Morgan & Alice Kluzer
#' @examples
#' NULL
#'
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
