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
#' @param genotypes (optional) An n X g \code{matrix} containing individual genotypes over g
#' SNPs/SNVs/other genetic variants. Rownames should correspond to the column names
#' of \code{nhoods(x)} and rownames of \code{design.df}. Genotypes should be in additive
#' format, i.e. AA=0, Aa=1, aa=2.
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
#'
#' @details
#' This function wraps up several steps of differential abundance testing using
#' the \code{edgeR} functions. These could be performed separately for users
#' who want to exercise more contol over their DA testing. By default this
#' function sets the \code{lib.sizes} to the colSums(x), and uses the
#' Quasi-Likelihood F-test in \code{glmQLFTest} for DA testing. FDR correction
#' is performed separately as the default multiple-testing correction is
#' inappropriate for neighbourhoods with overlapping cells.
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
#' @importFrom stats model.matrix
#' @importFrom Matrix colSums rowMeans
#' @importFrom stats dist median
#' @importFrom limma makeContrasts
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags calcNormFactors
testNhoods <- function(x, design, design.df, genotypes=NULL,
                       fdr.weighting=c("k-distance", "neighbour-distance", "max", "graph-overlap", "none"),
                       min.mean=0, model.contrasts=NULL, robust=TRUE, reduced.dim="PCA", REML=TRUE,
                       norm.method=c("TMM", "RLE", "logMS"), max.iters = 50, max.tol = 1e-5){
    is.lmm <- FALSE
    geno.only <- FALSE
    if(is(design, "formula")){
        # parse to find random and fixed effects
        parse <- unlist(strsplit(gsub(design, pattern="~", replacement=""), split= " + "))
        find_re <- any(grepl(parse, pattern="1*\\|"))
        if(find_re | !is.null(genotypes)){
            message("Random effects found")

            is.lmm <- TRUE
            if(find_re | is.null(genotypes)){
                # make model matrices for fixed and random effects
                z.model <- .parse_formula(design, design.df, vtype="re")
                rownames(z.model) <- rownames(design.df)
            } else if(find_re | !is.null(genotypes)){
                if(!all(rownames(genotypes) == rownames(design.df))){
                    stop("Genotype rownames do not match design.df rownames")
                }

                z.model <- .parse_formula(design, design.df, vtype="re")
                rownames(z.model) <- rownames(design.df)

                # rescale genotypes by 1/sqrt(m), where m = number of SNPs
                genotypes <- (genotypes - 1)/(sqrt(ncol(genotypes)))
            } else if(!find_re | !is.null(genotypes)){
                genotypes <- (genotypes - 1)/(sqrt(ncol(genotypes)))
                z.model <- diag(nrow(genotypes))
                colnames(z.model) <- paste0("Genetic", seq_len(nrow(genotypes)))
                rownames(z.model) <- rownames(design.df)
                geno.only <- TRUE
            }

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
            keep.nh <- rowMeans(nhoodCounts(x)[, rownames(x.model)]) >= min.mean
        } else{
            keep.nh <- rowMeans(nhoodCounts(x)) >= min.mean
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

    # estimate disperions _before_ all models
    dge <- estimateDisp(dge, x.model)

    if (is.lmm) {
        message("Running GLMM model - this may take a few minutes")
        
        if(isFALSE(geno.only)){
            rand.levels <- lapply(seq_along(colnames(z.model)), FUN=function(X) unique(z.model[, X]))
            names(rand.levels) <- colnames(z.model)
        } else{
            rand.levels <- list("Genetic"=colnames(z.model))
        }

        # extract tagwise dispersion for glmm
        dispersion <- dge$tagwise.dispersion
        offsets <- dge$samples$norm.factors
        glmm.cont <- list(theta.tol=max.tol, max.iter=max.iters)
        
        #wrapper function is the same for all analyses
        glmmWrapper <- function(y, dispersion, x.model, z.model, offsets, rand.levels, REML, glmm.control, geno.only=FALSE, kin=NULL){
            model.list <- NULL
            for (i in 1:nrow(y)) {
                model.list[[i]] <- tryCatch({fitGLMM(X=x.model, Z=z.model, y=y[i,], offsets=offsets, random.levels=rand.levels, REML = TRUE,
                                  dispersion=dispersion[i], geno.only=geno.only, Kin=kin, glmm.control=glmm.cont)},
                         error = function(err){
                             #message(paste("Neighborhood", i, "failed:", err))
                             return(list("FE"=NA, "RE"=NA, "Sigma"=NA,
                                         "converged"=NA, "Iters"=NA, "Dispersion"=NA,
                                         "Hessian"=NA, "SE"=NA, "t"=NA,
                                         "COEFF"=NA, "P"=NA, "Vpartial"=NA, "Ginv"=NA,
                                         "Vsinv"=NA, "Winv"=NA, "VCOV"=NA, "DF"=NA, "PVALS"=NA))
                         })
            }
            return(model.list)
        }

        if(!is.null(genotypes)){
            if(isTRUE(geno.only)){
                message("Running genetic model with ", nrow(genotypes), " genetic variants")
            } else{
                message("Running genetic model with ", nrow(z.model), " observations and ", ncol(genotypes), " genetic variants")
            }

            Kin <- genotypes %*% genotypes

            if(geno.only){
                fit <- glmmWrapper(y=dge$counts, dispersion = 1/dispersion, x.model, z.model, offsets, rand.levels, REML, glmm.control = glmm.cont, geno.only = geno.only, Kin=kin.matrix)    
            } else{
                fit <- glmmWrapper(y=dge$counts, dispersion = 1/dispersion, x.model, z.model, offsets, rand.levels, REML, glmm.control = glmm.cont, genotypes=genotypes, Kin=kin.matrix)    
            }
   
        } else{
            
            fit <- glmmWrapper(y=dge$counts, dispersion = 1/dispersion, x.model, z.model, offsets, rand.levels, REML, glmm.control = glmm.cont)    

            # give warning about how many neighborhoods didn't converge
            if (sum(!(unlist(lapply(fit, `[[`, "converged"))), na.rm = TRUE)/length(unlist(lapply(fit, `[[`, "converged"))) > 0){
                warning(paste(sum(!unlist(lapply(fit, `[[`, "converged")), na.rm = TRUE), "out of", length(unlist(lapply(fit, `[[`, "converged"))),"neighborhoods did not converge; increase number of iterations?"))

            }
            
            # res has to reflect output from glmQLFit
            res <- cbind.data.frame("Estimate" = unlist(lapply(lapply(fit, `[[`, "FE"), function(x) tail(x,1))), "Std. Error"= unlist(lapply(lapply(fit, `[[`, "SE"), function(x) tail(x,1))),
                                    "t value" = unlist(lapply(lapply(fit, `[[`, "t"), function(x) tail(x,1))), #"Df" = unlist(lapply(fit, `[[`, 11)),
                                    "PValue" = unlist(lapply(lapply(fit, `[[`, "PVALS"), function(x) tail(x,1))), matrix(unlist(lapply(fit, `[[`, "Sigma")), ncol=length(rand.levels), byrow=T),
                                    "Converged"=unlist(lapply(fit, `[[`, "converged")), "Dispersion" = unlist(lapply(fit, `[[`, "Dispersion")))
            rownames(res) <- 1:length(fit)
            colnames(res)[5:(5+length(rand.levels)-1)] <- paste(names(rand.levels), "variance")
        }

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
    res1 <- na.omit(res)
    mod.spatialfdr <- graphSpatialFDR(x.nhoods=nhoods(x),
                                      graph=graph(x),
                                      weighting=fdr.weighting,
                                      k=x@.k,
                                      pvalues=res1[order(res1$Nhood), ]$PValue,
                                      indices=nhoodIndex(x),
                                      distances=nhoodDistances(x),
                                      reduced.dimensions=reducedDim(x, reduced.dim))
    message(length(res$Estimate[order(res$Nhood)][!is.na(res$Estimate)]))
    message(length(mod.spatialfdr))
    res$SpatialFDR[order(res$Nhood)][!is.na(res$Estimate)] <- mod.spatialfdr
    res
}
