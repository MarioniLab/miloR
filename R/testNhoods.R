#' Perform differential neighbourhood abundance testing
#'
#' This will perform differential neighbourhood abundance testing after cell
#' counting.
#' @param x A \code{\linkS4class{Milo}} object with a non-empty
#' \code{nhoodCounts} slot.
#' @param design A \code{formula} or \code{model.matrix} object describing the
#' experimental design for differential abundance testing. The last component
#' of the formula or last column of the model matrix are by default the test
#' variable. This behaviour can be overridden by setting the \code{model.contrasts}
#' argument
#' @param design.df A \code{data.frame} containing meta-data to which \code{design}
#' refers to
#' @param kinship (optional) An n X n \code{matrix} containing pair-wise relationships between
#' observations, such as expected relationships or computed from SNPs/SNVs/other genetic variants.
#' Row names and column names should correspond to the column names of \code{nhoods(x)} and rownames
#' of \code{design.df}.
#' @param min.mean A scalar used to threshold neighbourhoods on the minimum
#' average cell counts across samples.
#' @param model.contrasts A string vector that defines the contrasts used to perform
#' DA testing.
#' @param fdr.weighting The spatial FDR weighting scheme to use. Choice from max,
#' neighbour-distance, graph-overlap or k-distance (default). If \code{none} is passed no
#' spatial FDR correction is performed and returns a vector of NAs.
#' @param robust If robust=TRUE then this is passed to edgeR and limma which use a robust
#' estimation for the global quasilikelihood dispersion distribution. See \code{edgeR} and
#' Phipson et al, 2013 for details.
#' @param norm.method A character scalar, either \code{"logMS"}, \code{"TMM"} or \code{"RLE"}.
#' The \code{"logMS"} method normalises the counts across samples using the log columns sums of
#' the count matrix as a model offset. \code{"TMM"} uses the trimmed mean of M-values normalisation
#' as described in Robinson & Oshlack, 2010, whilst \code{"RLE"} uses the relative log expression
#' method by Anders & Huber, 2010, to compute normalisation factors relative to a reference computed from
#' the geometric mean across samples.  The latter methods provides a degree of robustness against false positives
#' when there are very large compositional differences between samples.
#' @param reduced.dim A character scalar referring to the reduced dimensional slot used to compute distances for
#' the spatial FDR. This should be the same as used for graph building.
#' @param REML A logical scalar that controls the variance component behaviour to use either restricted maximum
#' likelihood (REML) or maximum likelihood (ML). The former is recommened to account for the bias in the ML
#' variance estimates.
#' @param glmm.solver A character scalar that determines which GLMM solver is applied. Must be one of: Fisher, HE
#' or HE-NNLS. HE or HE-NNLS are recommended when supplying a user-defined covariance matrix.
#' @param max.iters A scalar that determines the maximum number of iterations to run the GLMM solver if it does
#' not reach the convergence tolerance threshold.
#' @param max.tol A scalar that deterimines the GLMM solver convergence tolerance. It is recommended to keep
#' this number small to provide some confidence that the parameter estimates are at least in a feasible region
#' and close to a \emph{local} optimum
#' @param subset.nhoods A character, numeric or logical vector that will subset the analysis to the specific nhoods. If
#' a character vector these should correspond to row names of \code{nhoodCounts}. If a logical vector then
#' these should have the same \code{length} as \code{nrow} of \code{nhoodCounts}. If numeric, then these are assumed
#' to correspond to indices of \code{nhoodCounts} - if the maximal index is greater than \code{nrow(nhoodCounts(x))}
#' an error will be produced.
#' @param fail.on.error A logical scalar the determines the behaviour of the error reporting. Used for debugging only.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the arguments for parallelisation. By default
#' this will evaluate using \code{SerialParam()}. See \code{details}on how to use parallelisation in \code{testNhoods}.
#'
#' @details
#' This function wraps up several steps of differential abundance testing using
#' the \code{edgeR} functions. These could be performed separately for users
#' who want to exercise more contol over their DA testing. By default this
#' function sets the \code{lib.sizes} to the colSums(x), and uses the
#' Quasi-Likelihood F-test in \code{glmQLFTest} for DA testing. FDR correction
#' is performed separately as the default multiple-testing correction is
#' inappropriate for neighbourhoods with overlapping cells.
#' The GLMM testing cannot be performed using \code{edgeR}, however, a separate
#' function \code{fitGLMM} can be used to fit a mixed effect model to each
#' nhood (see \code{fitGLMM} docs for details).
#'
#' Parallelisation is currently only enabled for the NB-LMM and uses the BiocParallel paradigm. In
#' general the GLM implementation in \code{glmQLFit} is sufficiently fast that it does not require
#' parallelisation. Parallelisation requires the user to pass a \linkS4class{BiocParallelParam} object
#' with the parallelisation arguments contained therein. This relies on the user to specify how
#' parallelisation - for details see the \code{BiocParallel} package.
#'
#' @return A \code{data.frame} of model results, which contain:
#' \describe{
#' \item{\code{logFC}:}{Numeric, the log fold change between conditions, or for
#' an ordered/continous variable the per-unit
#' change in (normalized) cell counts per unit-change in experimental variable.}
#' \item{\code{logCPM}:}{Numeric, the log counts per million (CPM), which equates
#' to the average log normalized cell counts
#' across all samples.}
#' \item{\code{F}:}{Numeric, the F-test statistic from the quali-likelihood F-test
#' implemented in \code{edgeR}.}
#' \item{\code{PValue}:}{Numeric, the unadjusted p-value from the quasi-likelihood F-test.}
#' \item{\code{FDR}:}{Numeric, the Benjamini & Hochberg false discovery weight
#' computed from \code{p.adjust}.}
#' \item{\code{Nhood}:}{Numeric, a unique identifier corresponding to the specific
#' graph neighbourhood.}
#' \item{\code{SpatialFDR}:}{Numeric, the weighted FDR, computed to adjust for spatial
#' graph overlaps between neighbourhoods. For details see \link{graphSpatialFDR}.}
#' }
#'
#' @author Mike Morgan
#'
#' @examples
#' library(SingleCellExperiment)
#' ux.1 <- matrix(rpois(12000, 5), ncol=400)
#' ux.2 <- matrix(rpois(12000, 4), ncol=400)
#' ux <- rbind(ux.1, ux.2)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#'
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'                             reducedDims=SimpleList(PCA=pca$x))
#'
#' milo <- Milo(sce)
#' milo <- buildGraph(milo, k=20, d=10, transposed=TRUE)
#' milo <- makeNhoods(milo, k=20, d=10, prop=0.3)
#' milo <- calcNhoodDistance(milo, d=10)
#'
#' cond <- rep("A", ncol(milo))
#' cond.a <- sample(1:ncol(milo), size=floor(ncol(milo)*0.25))
#' cond.b <- setdiff(1:ncol(milo), cond.a)
#' cond[cond.b] <- "B"
#' meta.df <- data.frame(Condition=cond, Replicate=c(rep("R1", 132), rep("R2", 132), rep("R3", 136)))
#' meta.df$SampID <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
#' milo <- countCells(milo, meta.data=meta.df, samples="SampID")
#'
#' test.meta <- data.frame("Condition"=c(rep("A", 3), rep("B", 3)), "Replicate"=rep(c("R1", "R2", "R3"), 2))
#' test.meta$Sample <- paste(test.meta$Condition, test.meta$Replicate, sep="_")
#' rownames(test.meta) <- test.meta$Sample
#' da.res <- testNhoods(milo, design=~Condition, design.df=test.meta[colnames(nhoodCounts(milo)), ], norm.method="TMM")
#' da.res
#'
#' @name testNhoods
NULL


#' @export
#' @importFrom Matrix colSums rowMeans
#' @importFrom MatrixGenerics colSums2
#' @importFrom utils tail
#' @importFrom stats dist median model.matrix
#' @importFrom limma makeContrasts
#' @importFrom BiocParallel bplapply SerialParam bptry bpok bpoptions bpstopOnError
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags calcNormFactors
testNhoods <- function(x, design, design.df, kinship=NULL,
                       fdr.weighting=c("k-distance", "neighbour-distance", "max", "graph-overlap", "none"),
                       min.mean=0, model.contrasts=NULL, robust=TRUE, reduced.dim="PCA", REML=TRUE,
                       norm.method=c("TMM", "RLE", "logMS"), max.iters = 50, max.tol = 1e-5, glmm.solver=NULL,
                       subset.nhoods=NULL, fail.on.error=FALSE, BPPARAM=SerialParam()){
    is.lmm <- FALSE
    geno.only <- FALSE
    if(is(design, "formula")){
        # parse to find random and fixed effects - need to double check the formula is valid
        parse <- unlist(strsplit(gsub(design, pattern="~", replacement=""), split= "+", fixed=TRUE))
        if(length(parse) < 1){
            stop("Forumla ", design, " not a proper formula - variables should be separated by '+'")
        }

        find_re <- any(grepl(parse, pattern="1\\W?\\|"))
        if(!find_re & any(grepl(parse, pattern="[0-9]?\\W?\\|\\W?"))){
            stop(design, " is an invalid formula for random effects. Use (1 | VARIABLE) format.")
        }

        if(find_re | !is.null(kinship)){
            message("Random effects found")

            is.lmm <- TRUE
            if(find_re | is.null(kinship)){
                # make model matrices for fixed and random effects
                z.model <- .parse_formula(design, design.df, vtype="re")
                rownames(z.model) <- rownames(design.df)
            } else if(find_re & !is.null(kinship)){
                if(!all(rownames(kinship) == rownames(design.df))){
                    stop("Genotype rownames do not match design.df rownames")
                }

                z.model <- .parse_formula(design, design.df, vtype="re")
                rownames(z.model) <- rownames(design.df)
            } else if(!find_re & !is.null(kinship)){
                z.model <- diag(nrow(kinship))
                colnames(z.model) <- paste0("Genetic", seq_len(nrow(kinship)))
                rownames(z.model) <- rownames(design.df)
                geno.only <- TRUE
            }

            # this will always implicitly include an intercept term - perhaps
            # this shouldn't be the case?
            x.model <- .parse_formula(design, design.df, vtype="fe")
            rownames(x.model) <- rownames(design.df)
            max.iters <- max.iters

            if(all(rownames(x.model) != rownames(z.model))){
                stop("Discordant sample names for mixed model design matrices")
            }

        } else{
            x.model <- model.matrix(design, data=design.df)
            rownames(x.model) <- rownames(design.df)
        }
    } else if(is(design, "matrix")){
        x.model <- design
        if(nrow(x.model) != nrow(design.df)){
            stop("Design matrix and model matrix are not the same dimensionality")
        }

        if(any(rownames(x.model) != rownames(design.df))){
            warning("Design matrix and model matrix dimnames are not the same")
            # check if rownames are a subset of the design.df
            check.names <- any(rownames(x.model) %in% rownames(design.df))
            if(isTRUE(check.names)){
                rownames(x.model) <- rownames(design.df)
            } else{
                stop("Design matrix and model matrix rownames are not a subset")
            }
        }
    } else{
        stop("design must be either a formula or model matrix")
    }

    if(!is(x, "Milo")){
        stop("Unrecognised input type - must be of class Milo")
    } else if(.check_empty(x, "nhoodCounts")){
        stop("Neighbourhood counts missing - please run countCells first")
    }

    if(!any(norm.method %in% c("TMM", "logMS", "RLE"))){
        stop("Normalisation method ", norm.method, " not recognised. Must be either TMM, RLE or logMS")
    }

    if(!reduced.dim %in% reducedDimNames(x)){
        stop(reduced.dim, " is not found in reducedDimNames. Avaiable options are ", paste(reducedDimNames(x), collapse=","))
    }

    subset.counts <- FALSE
    if(ncol(nhoodCounts(x)) != nrow(x.model)){
        # need to allow for design.df with a subset of samples only
        if(all(rownames(x.model) %in% colnames(nhoodCounts(x)))){
            message("Design matrix is a strict subset of the nhood counts")
            subset.counts <- TRUE
        } else{
            stop("Design matrix (", nrow(x.model), ") and nhood counts (", ncol(nhoodCounts(x)), ") are not the same dimension")
        }
    }

    if(is.lmm){
        if(ncol(nhoodCounts(x)) != nrow(z.model)){
            # need to allow for design.df with a subset of samples only
            if(all(rownames(z.model) %in% colnames(nhoodCounts(x)))){
                message("Random effects design matrix is a strict subset of the nhood counts")
                subset.counts <- TRUE
            } else{
                stop("Random effects design matrix (", nrow(z.model), ") and nhood counts (",
                     ncol(nhoodCounts(x)), ") are not the same dimension")
            }
        }
    }

    # assume nhoodCounts and model are in the same order
    # cast as DGEList doesn't accept sparse matrices
    # what is the cost of cast a matrix that is already dense vs. testing it's class
    if(min.mean > 0){
        if(isTRUE(subset.counts)){
            if(!is.null(subset.nhoods)){
                mean.keep <- rowMeans(nhoodCounts(x)[, rownames(x.model)]) >= min.mean
                if(is.character(subset.nhoods)){
                    if(all(subset.nhoods) %in% rownames(nhoodCounts(x))){
                        keep.nh <- rownames(nhoodCounts(x) %in% subset.nhoods)
                    } else{
                        stop("Nhood subsetting is illegal - use same names as in rownames of nhoodCounts")
                    }
                } else if(is.logical(subset.nhoods)){
                    if(length(subset.nhoods) != nrow(nhoodCounts(x))){
                        stop("Logical subset vector must be same length as nrow nhoodCounts")
                    } else{
                        keep.nh <- subset.nhoods
                    }
                } else if(is.numeric(subset.nhoods)){
                    if(max(subset.nhoods) > nrow(nhoodCounts(x))){
                        stop("Maximum index is out of bounds: ", max(subset.nhoods))
                    } else{
                        keep.nh <- seq_len(nrow(nhoodCounts(x))) %in% subset.nhoods
                    }
                } else{
                    stop("Subsetting vector type not recognised: ", type(subset.nhoods))
                }

                keep.nh <- mean.keep & keep.nh

            } else{
                keep.nh <- rowMeans(nhoodCounts(x)[, rownames(x.model)]) >= min.mean
            }
        } else if(!is.null(subset.nhoods)){
            mean.keep <- rowMeans(nhoodCounts(x)[, rownames(x.model)]) >= min.mean
            if(is.character(subset.nhoods)){
                if(all(subset.nhoods) %in% rownames(nhoodCounts(x))){
                    keep.nh <- rownames(nhoodCounts(x) %in% subset.nhoods)
                } else{
                    stop("Nhood subsetting is illegal - use same names as in rownames of nhoodCounts")
                }
            } else if(is.logical(subset.nhoods)){
                if(length(subset.nhoods) != nrow(nhoodCounts(x))){
                    stop("Logical subset vector must be same length as nrow nhoodCounts")
                } else{
                    keep.nh <- subset.nhoods
                }
            } else if(is.numeric(subset.nhoods)){
                if(max(subset.nhoods) > nrow(nhoodCounts(x))){
                    stop("Maximum index is out of bounds: ", max(subset.nhoods))
                } else{
                    keep.nh <- seq_len(nrow(nhoodCounts(x))) %in% subset.nhoods
                }
            } else{
                stop("Subsetting vector type not recognised: ", type(subset.nhoods))
            }

            keep.nh <- mean.keep & keep.nh
        } else{
            keep.nh <- rowMeans(nhoodCounts(x)) >= min.mean
        }
    } else if(!is.null(subset.nhoods)){
        if(is.character(subset.nhoods)){
            if(all(subset.nhoods) %in% rownames(nhoodCounts(x))){
                keep.nh <- rownames(nhoodCounts(x) %in% subset.nhoods)
            } else{
                stop("Nhood subsetting is illegal - use same names as in rownames of nhoodCounts")
            }
        } else if(is.logical(subset.nhoods)){
            if(length(subset.nhoods) != nrow(nhoodCounts(x))){
                stop("Logical subset vector must be same length as nrow nhoodCounts")
            } else{
                keep.nh <- subset.nhoods
            }
        } else if(is.numeric(subset.nhoods)){
            if(max(subset.nhoods) > nrow(nhoodCounts(x))){
                stop("Maximum index is out of bounds: ", max(subset.nhoods))
            } else{
                keep.nh <- seq_len(nrow(nhoodCounts(x))) %in% subset.nhoods
            }
        } else{
            stop("Subsetting vector type not recognised: ", type(subset.nhoods))
        }

    } else{
        if(isTRUE(subset.counts)){
            keep.nh <- rep(TRUE, nrow(nhoodCounts(x)[, rownames(x.model)]))
        }else{
            keep.nh <- rep(TRUE, nrow(nhoodCounts(x)))
        }
    }

    if(isTRUE(subset.counts)){
        keep.samps <- intersect(rownames(x.model), colnames(nhoodCounts(x)[keep.nh, ]))
    } else{
        keep.samps <- colnames(nhoodCounts(x)[keep.nh, ])
    }

    if(any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) != rownames(x.model)) & !any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) %in% rownames(x.model))){
        stop("Sample names in design matrix and nhood counts are not matched.
             Set appropriate rownames in design matrix.")
    } else if(any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) != rownames(x.model)) & any(colnames(nhoodCounts(x)[keep.nh, keep.samps]) %in% rownames(x.model))){
        warning("Sample names in design matrix and nhood counts are not matched. Reordering")
        x.model <- x.model[colnames(nhoodCounts(x)[keep.nh, keep.samps]), ]
        if(is.lmm){
            z.model <- z.model[colnames(nhoodCounts(x)[keep.nh, keep.samps]), , drop = FALSE]
        }
    }

    if(length(norm.method) > 1){
        message("Using TMM normalisation")
        dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                       lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
        dge <- calcNormFactors(dge, method="TMM")
    } else if(norm.method %in% c("TMM")){
        message("Using TMM normalisation")
        dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                       lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
        dge <- calcNormFactors(dge, method="TMM")
    } else if(norm.method %in% c("RLE")){
        message("Using RLE normalisation")
        dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                       lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
        dge <- calcNormFactors(dge, method="RLE")
    }else if(norm.method %in% c("logMS")){
        message("Using logMS normalisation")
        dge <- DGEList(counts=nhoodCounts(x)[keep.nh, keep.samps],
                       lib.size=colSums(nhoodCounts(x)[keep.nh, keep.samps]))
    }

    dge <- estimateDisp(dge, x.model)

    if (is.lmm) {
        message("Running GLMM model - this may take a few minutes")

        if(isFALSE(geno.only)){
            rand.levels <- lapply(seq_along(colnames(z.model)), FUN=function(X) {
                zx <- unique(z.model[, X])
                if(is.numeric(zx)){
                    paste0(colnames(z.model)[X], zx)
                } else{
                    zx
                }
            })
            names(rand.levels) <- colnames(z.model)
        } else{
            rand.levels <- list("Genetic"=colnames(z.model))
        }

        # extract tagwise dispersion for glmm
        # re-scale these to allow for non-zero variances
        dispersion <- dge$tagwise.dispersion

        # I think these need to be logged
        offsets <- log(dge$samples$norm.factors)
        glmm.cont <- list(theta.tol=max.tol, max.iter=max.iters, solver=glmm.solver)

        #wrapper function is the same for all analyses
        glmmWrapper <- function(Y, disper, Xmodel, Zmodel, off.sets, randlevels, reml, glmm.contr, genonly=FALSE, kin.ship=NULL, BPPARAM=BPPARAM, error.fail=FALSE){
            #bp.list <- NULL
            # this needs to be able to run with BiocParallel
            bp.list <- bptry({bplapply(seq_len(nrow(Y)), BPOPTIONS=bpoptions(stop.on.error = FALSE),
                                         FUN=function(i, Xmodel, Zmodel, Y, off.sets,
                                                      randlevels, disper, genonly,
                                                      kin.ship, glmm.contr, reml){
                                             fitGLMM(X=Xmodel, Z=Zmodel, y=Y[i, ], offsets=off.sets,
                                                     random.levels=randlevels, REML = reml,
                                                     dispersion=disper[i], geno.only=genonly,
                                                     Kin=kinship, glmm.control=glmm.contr)
                                             }, BPPARAM=BPPARAM,
                                         Xmodel=Xmodel, Zmodel=Zmodel, Y=Y, off.sets=off.sets,
                                         randlevels=randlevels, disper=disper, genonly=genonly,
                                         kin.ship=kin.ship, glmm.cont=glmm.cont, reml=reml, error.fail=fail.on.error)
                                }) # need to handle this output which is a bplist_error object

            # parse the bplist_error object
            if(all(bpok(bp.list))){
                model.list <- as(bp.list, "list")
            } else{
                model.list <- list()
                # set the failed results to all NA
                for(x in seq_along(bp.list)){
                    if(!bpok(bp.list)[x]){
                        bperr <- attr(bp.list[[x]], "traceback")
                        if(isTRUE(error.fail)){
                            stop(bperr)
                        }

                        model.list[[x]] <- list("FE"=NA, "RE"=NA, "Sigma"=NA,
                                                "converged"=FALSE, "Iters"=NA, "Dispersion"=NA,
                                                "Hessian"=NA, "SE"=NA, "t"=NA, "PSVAR"=NA,
                                                "COEFF"=NA, "P"=NA, "Vpartial"=NA, "Ginv"=NA,
                                                "Vsinv"=NA, "Winv"=NA, "VCOV"=NA, "LOGLIHOOD"=NA,
                                                "DF"=NA, "PVALS"=NA,
                                                "ERROR"=bperr)
                    }else{
                        model.list[[x]] <- bp.list[[x]]
                    }
                }
            }
            return(model.list)
        }


        if(!is.null(kinship)){
            if(isTRUE(geno.only)){
                message("Running genetic model with ", nrow(kinship), " individuals")
            } else{
                message("Running genetic model with ", nrow(z.model), " observations")
            }

            if(geno.only){
                fit <- glmmWrapper(Y=dge$counts, disper = 1/dispersion, Xmodel=x.model, Zmodel=z.model,
                                   off.sets=offsets, randlevels=rand.levels, reml=REML, glmm.contr = glmm.cont,
                                   genonly = geno.only, kin.ship=kinship,
                                   BPPARAM=BPPARAM, error.fail=fail.on.error)
            } else{
                fit <- glmmWrapper(Y=dge$counts, disper = 1/dispersion, Xmodel=x.model, Zmodel=z.model,
                                   off.sets=offsets, randlevels=rand.levels, reml=REML, glmm.contr = glmm.cont,
                                   genonly = geno.only, kin.ship=kinship,
                                   BPPARAM=BPPARAM, error.fail=fail.on.error)
            }

        } else{

            fit <- glmmWrapper(Y=dge$counts, disper = 1/dispersion, Xmodel=x.model, Zmodel=z.model,
                               off.sets=offsets, randlevels=rand.levels, reml=REML, glmm.contr = glmm.cont,
                               genonly = geno.only, kin.ship=kinship,
                               BPPARAM=BPPARAM)
        }

        # give warning about how many neighborhoods didn't converge and error is all nhoods failed
        if (sum(!(unlist(lapply(fit, `[[`, "converged"))), na.rm = TRUE)/length(unlist(lapply(fit, `[[`, "converged"))) > 0){
            if(all(is.na(unlist(lapply(fit, `[[`, "logFC"))))){
                err.list <- unique(unlist(lapply(fit, `[[`, "ERROR")))
                n.err <- length(err.list)
                stop("Lowest traceback returned: ", paste(err.list[1:3], err.list[c((n.err-3):(n.err))], collapse="\n")) # the first and last 3 traceback steps
            } else{
                warning(paste(sum(!unlist(lapply(fit, `[[`, "converged")), na.rm = TRUE), "out of", length(unlist(lapply(fit, `[[`, "converged"))),
                              "neighborhoods did not converge; increase number of iterations?"))
            }

        }

        # res has to reflect output from glmQLFit - express variance as a proportion as well.
        # this only reports the final fixed effect parameter
        ret.beta <- ncol(x.model)

        res <- cbind.data.frame("logFC" = unlist(lapply(lapply(fit, `[[`, "FE"), function(x) x[ret.beta])),
                                "logCPM"=log2((rowMeans(nhoodCounts(x)[keep.nh, ]/colSums2(nhoodCounts(x))))*1e6),
                                "SE"= unlist(lapply(lapply(fit, `[[`, "SE"), function(x) x[ret.beta])),
                                "tvalue" = unlist(lapply(lapply(fit, `[[`, "t"), function(x) x[ret.beta])),
                                "PValue" = unlist(lapply(lapply(fit, `[[`, "PVALS"), function(x) x[ret.beta])),
                                matrix(unlist(lapply(fit, `[[`, "Sigma")), ncol=length(rand.levels), byrow=TRUE),
                                "Converged"=unlist(lapply(fit, `[[`, "converged")), "Dispersion" = unlist(lapply(fit, `[[`, "Dispersion")),
                                "Logliklihood"=unlist(lapply(fit, `[[`, "LOGLIHOOD")))

        rownames(res) <- 1:length(fit)
        colnames(res)[6:(6+length(rand.levels)-1)] <- paste(names(rand.levels), "variance")
    } else {

        fit <- glmQLFit(dge, x.model, robust=robust)
        if(!is.null(model.contrasts)){
            mod.constrast <- makeContrasts(contrasts=model.contrasts, levels=x.model)
            res <- as.data.frame(topTags(glmQLFTest(fit, contrast=mod.constrast),
                                         sort.by='none', n=Inf))
        } else{
            n.coef <- ncol(x.model)
            res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))
        }
    }

    res$Nhood <- as.numeric(rownames(res))
    message("Performing spatial FDR correction with ", fdr.weighting[1], " weighting")
    # res1 <- na.omit(res)
    mod.spatialfdr <- graphSpatialFDR(x.nhoods=nhoods(x),
                                      graph=graph(x),
                                      weighting=fdr.weighting,
                                      k=x@.k,
                                      pvalues=res[order(res$Nhood), ]$PValue,
                                      indices=nhoodIndex(x),
                                      distances=nhoodDistances(x),
                                      reduced.dimensions=reducedDim(x, reduced.dim))
    res$SpatialFDR[order(res$Nhood)] <- mod.spatialfdr
    res
}
