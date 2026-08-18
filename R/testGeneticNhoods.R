#' Genome-wide cell state QTL testing on neighbourhoods
#'
#' Test for association between genetic variants and neighbourhood cell
#' abundance across a set of variants, using a negative binomial generalised
#' linear mixed model with a genetic relationship matrix.
#'
#' This function is intended for genome-wide analyses, where the same set of
#' neighbourhoods is tested against many thousands to millions of variants. It
#' is deliberately separate from \code{\link{testNhoods}}: \code{testNhoods}
#' fits one model per neighbourhood for a single variable of interest and
#' returns Wald statistics, whereas \code{testGeneticNhoods} fits a null model
#' once per neighbourhood and then evaluates a score test for every variant
#' against that null. Use \code{\link{testNhoods}} for designs with a single
#' variable of interest such as age or disease status.
#'
#' @param x A \code{\linkS4class{Milo}} object with a non-empty
#' \code{nhoodCounts} slot.
#' @param design A \code{formula} describing the null model, i.e. the fixed and
#' random effects \emph{excluding} the genetic variants. Random effects are
#' given in \code{lme4} syntax, e.g. \code{~ sex + age + (1|batch)}.
#' @param design.df A \code{data.frame} containing the model variables, with
#' rownames matching the columns of \code{nhoodCounts(x)}.
#' @param genotypes A \code{matrix} of genotypes with one row per sample and one
#' column per variant. Rownames must match the rownames of \code{design.df}.
#' Genotypes are usually allele dosages coded 0, 1 and 2.
#' @param kinship A \code{matrix} genetic relationship matrix with rownames and
#' colnames matching the rownames of \code{design.df}.
#' @param reduced.dim A \code{character} scalar naming the reduced dimensional
#' slot used to build the graph. Used for the spatial FDR weighting.
#' @param fdr.weighting A \code{character} scalar giving the spatial FDR
#' weighting scheme, passed to \code{\link{graphSpatialFDR}}. One of
#' \code{"graph-overlap"}, \code{"k-distance"}, \code{"neighbour-distance"},
#' \code{"max"} or \code{"none"}.
#' @param norm.method A \code{character} scalar giving the normalisation used to
#' compute model offsets. One of \code{"TMM"}, \code{"RLE"} or \code{"logMS"}.
#' @param cell.sizes A named \code{numeric} vector of column sizes used for
#' normalisation. Defaults to the neighbourhood count column sums.
#' @param min.count A \code{numeric} scalar giving the minimum average count for
#' a neighbourhood to be tested.
#' @param REML A \code{logical} scalar controlling whether the variance
#' components are estimated by REML.
#' @param block.size A \code{numeric} scalar giving the number of variants
#' evaluated per matrix multiplication. Larger blocks are faster but use more
#' memory; the working set is roughly \code{n * block.size} doubles.
#' @param max.iters A \code{numeric} scalar giving the maximum number of PQL
#' iterations for each null model fit.
#' @param theta.conv A \code{numeric} scalar giving the convergence tolerance
#' for the null model parameter estimates.
#' @param sigma.boundary A \code{numeric} scalar. Neighbourhoods whose genetic
#' variance component falls below this value are flagged in the output as
#' having hit the boundary of the parameter space.
#' @param BPPARAM A \code{\linkS4class{BiocParallelParam}} object controlling
#' parallelisation across neighbourhoods.
#'
#' @details The model fitted for each neighbourhood is the same negative
#' binomial GLMM used elsewhere in Milo. Because the variance components do not
#' depend on the variant being tested, the null model is fitted once and the
#' REML projection \eqn{P} and the vector \eqn{P y^*} are reused for every
#' variant. For a variant with genotype vector \eqn{g} the score statistic is
#' \eqn{U = g' P y^*} with variance \eqn{I = g' P g}, giving
#' \eqn{\chi^2_1 = U^2/I} and a one-step effect estimate \eqn{U/I}. This is the
#' approach used by genome-wide mixed model methods such as EMMAX and fastGWA.
#'
#' Note that the genetic variance component is only weakly identified when the
#' genetic relationship matrix is close to the identity, because in that case
#' \eqn{\sigma_g I} is confounded with the negative binomial dispersion. This is
#' the expected situation for a cohort of unrelated donors. The
#' \code{SigmaGenetic} and \code{SigmaBoundary} columns of the output report the
#' fitted value and flag neighbourhoods where it has collapsed, so that the
#' behaviour is visible rather than silent.
#'
#' Variants passing a screening threshold should be refitted exactly with
#' \code{\link{refineGeneticHits}} to obtain final effect size estimates, since
#' the score test holds the variance components at their null values.
#'
#' @return A \code{data.frame} with one row per tested neighbourhood and variant
#' combination, containing:
#' \describe{
#' \item{\code{Nhood}:}{the neighbourhood index.}
#' \item{\code{SNP}:}{the variant identifier, taken from the genotype colnames.}
#' \item{\code{logFC}:}{the one-step log fold change estimate.}
#' \item{\code{SE}:}{the standard error of the log fold change.}
#' \item{\code{Chisq}:}{the score test statistic on 1 degree of freedom.}
#' \item{\code{PValue}:}{the unadjusted p-value.}
#' \item{\code{SpatialFDR}:}{the spatial FDR, computed across neighbourhoods separately for each variant.}
#' \item{\code{SigmaGenetic}:}{the genetic variance component of the null model.}
#' \item{\code{SigmaBoundary}:}{\code{TRUE} where the genetic variance component collapsed to the boundary.}
#' \item{\code{Dispersion}:}{the null model dispersion estimate.}
#' \item{\code{NullConverged}:}{whether the null model met the convergence tolerance.}
#' }
#'
#' @author Mike Morgan
#'
#' @examples
#' library(SingleCellExperiment)
#' ux <- matrix(rpois(12000, 5), ncol=400)
#' vx <- log2(ux + 1)
#' pca <- prcomp(t(vx))
#' sce <- SingleCellExperiment(assays=list(counts=ux, logcounts=vx),
#'                             reducedDims=SimpleList(PCA=pca$x))
#' milo <- Milo(sce)
#' milo <- buildGraph(milo, k=20, d=10, transposed=TRUE)
#' milo <- makeNhoods(milo, k=20, d=10, prop=0.3)
#' milo <- calcNhoodDistance(milo, d=10)
#' meta.df <- data.frame(Sample=rep(paste0("S", seq_len(10)), each=40),
#'                       Sex=rep(c(0, 1), 200))
#' milo <- countCells(milo, meta.data=meta.df, samples="Sample")
#' dd <- data.frame(Sex=rep(c(0, 1), 5))
#' rownames(dd) <- paste0("S", seq_len(10))
#' geno <- matrix(rbinom(20, 2, 0.4), nrow=10,
#'                dimnames=list(rownames(dd), c("rs1", "rs2")))
#' kin <- diag(10)
#' dimnames(kin) <- list(rownames(dd), rownames(dd))
#' res <- testGeneticNhoods(milo, design=~Sex, design.df=dd,
#'                          genotypes=geno, kinship=kin,
#'                          fdr.weighting="graph-overlap")
#' head(res)
#'
#' @name testGeneticNhoods
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @importFrom stats model.matrix pchisq as.formula terms
#' @importFrom Matrix rowMeans
#' @importFrom SingleCellExperiment reducedDim
#' @export
testGeneticNhoods <- function(x, design, design.df, genotypes, kinship,
                              reduced.dim="PCA", fdr.weighting="graph-overlap",
                              norm.method="TMM", cell.sizes=NULL, min.count=1,
                              REML=TRUE, block.size=1000, max.iters=100,
                              theta.conv=1e-6, sigma.boundary=1e-6,
                              BPPARAM=BiocParallel::SerialParam()){

    if(!is(x, "Milo")){
        stop("Unrecognised input type - must be of class Milo")
    }

    if(ncol(nhoodCounts(x)) == 1 & nrow(nhoodCounts(x)) == 1){
        stop("Neighbourhood counts missing - please run countCells first")
    }

    if(!is(design, "formula")){
        stop("design must be a formula")
    }

    if(is.null(rownames(design.df))){
        stop("design.df must have rownames matching the columns of nhoodCounts(x)")
    }

    genotypes <- as.matrix(genotypes)
    if(is.null(colnames(genotypes))){
        colnames(genotypes) <- paste0("SNP", seq_len(ncol(genotypes)))
    }

    if(is.null(rownames(genotypes))){
        stop("genotypes must have rownames matching the rownames of design.df")
    }

    kinship <- as.matrix(kinship)
    if(nrow(kinship) != ncol(kinship)){
        stop("kinship must be a square matrix")
    }

    # ---- align samples across every input --------------------------------
    keep.samps <- intersect(colnames(nhoodCounts(x)), rownames(design.df))
    keep.samps <- intersect(keep.samps, rownames(genotypes))
    keep.samps <- intersect(keep.samps, rownames(kinship))

    if(length(keep.samps) < 2){
        stop("Fewer than 2 samples shared between nhoodCounts, design.df, genotypes and kinship")
    }

    design.df <- design.df[keep.samps, , drop=FALSE]
    genotypes <- genotypes[keep.samps, , drop=FALSE]
    kinship <- kinship[keep.samps, keep.samps, drop=FALSE]
    n.samps <- length(keep.samps)

    # ---- split the formula into fixed and random parts --------------------
    has.re <- any(grepl("\\|", attr(terms(design), "term.labels")))
    fixed.form <- .dropRandomTerms(design)
    x.model <- model.matrix(fixed.form, data=design.df)
    rownames(x.model) <- rownames(design.df)

    if(isTRUE(has.re)){
        z.model <- as.matrix(.parse_formula(design, design.df, vtype="re"))
        rownames(z.model) <- rownames(design.df)
        z.model <- z.model[keep.samps, , drop=FALSE]
        z.full <- .expandRandomDesign(z.model)
        u.indices <- attr(z.full, "u_indices")
    } else {
        z.full <- matrix(0, nrow=n.samps, ncol=0)
        u.indices <- list()
    }

    n.re <- length(u.indices) + 1L # genetic component is always last

    # ---- offsets ----------------------------------------------------------
    counts <- nhoodCounts(x)[, keep.samps, drop=FALSE]
    keep.nh <- Matrix::rowMeans(counts) >= min.count
    if(sum(keep.nh) == 0){
        stop("No neighbourhoods pass min.count - lower the threshold")
    }
    counts <- counts[keep.nh, , drop=FALSE]

    if(is.null(cell.sizes)){
        cell.sizes <- colSums(nhoodCounts(x)[, keep.samps, drop=FALSE])
    }

    dge <- DGEList(counts=counts, lib.size=cell.sizes)
    if(norm.method %in% c("TMM", "RLE")){
        dge <- calcNormFactors(dge, method=norm.method)
    }
    offsets <- log(dge$samples$lib.size * dge$samples$norm.factors)
    dge <- estimateDisp(dge, x.model)
    init.disp <- max(1e-2, 1 / mean(dge$tagwise.dispersion, na.rm=TRUE))

    # ---- K inverse is invariant across every neighbourhood and variant ----
    kin.inv <- tryCatch(chol2inv(chol(kinship)),
                        error=function(e){
                            warning("Kinship matrix is not positive definite - using a generalised inverse")
                            ek <- eigen(kinship, symmetric=TRUE)
                            pos <- ek$values > max(ek$values) * 1e-9
                            ek$vectors[, pos, drop=FALSE] %*%
                                diag(1 / ek$values[pos], sum(pos)) %*%
                                t(ek$vectors[, pos, drop=FALSE])
                        })

    # unname: which() carries names through, and a named length-1 column makes
    # data.frame() warn about row names from a short variable
    nhood.ids <- unname(which(keep.nh))
    n.snp <- ncol(genotypes)
    blocks <- split(seq_len(n.snp), ceiling(seq_len(n.snp) / block.size))

    message("Fitting ", length(nhood.ids), " null models over ", n.snp, " variants")

    # ---- one null model per neighbourhood, then score every variant -------
    per.nhood <- bplapply(seq_along(nhood.ids), BPPARAM=BPPARAM, FUN=function(i){
        yi <- as.numeric(counts[i, ])
        mu0 <- rep(max(mean(yi), 1e-2), n.samps)

        nullfit <- tryCatch(
            fitGeneticNullGlmm(Z=z.full, X=x.model, K=kinship,
                               muvec=mu0, offsets=offsets,
                               curr_beta=c(log(max(mean(yi), 1e-2)) - mean(offsets),
                                           rep(0, ncol(x.model) - 1)),
                               curr_u=rep(0, ncol(z.full) + n.samps),
                               curr_sigma=rep(0.5, n.re),
                               y=yi, u_indices=u.indices,
                               theta_conv=theta.conv, curr_disp=init.disp,
                               REML=REML, maxit=max.iters, Kinv_=kin.inv,
                               return_projection=TRUE),
            error=function(e) NULL)

        if(is.null(nullfit)){
            return(NULL)
        }

        sigmas <- as.numeric(nullfit$Sigma)
        sigma.g <- sigmas[length(sigmas)]

        stats <- lapply(blocks, FUN=function(bl){
            gb <- genotypes[, bl, drop=FALSE]
            scoreTestGeneticSNPs(P=nullfit$P,
                                 Pystar=as.numeric(nullfit$Pystar),
                                 G=gb)
        })

        data.frame(Nhood=nhood.ids[i],
                   SNP=colnames(genotypes),
                   logFC=unlist(lapply(stats, function(s) as.numeric(s$Beta))),
                   SE=unlist(lapply(stats, function(s) as.numeric(s$SE))),
                   Chisq=unlist(lapply(stats, function(s) as.numeric(s$Chisq))),
                   SigmaGenetic=sigma.g,
                   SigmaBoundary=sigma.g <= sigma.boundary,
                   Dispersion=as.numeric(nullfit$Dispersion),
                   NullConverged=as.logical(nullfit$converged),
                   stringsAsFactors=FALSE)
    })

    failed <- vapply(per.nhood, is.null, logical(1))
    if(all(failed)){
        stop("All null models failed to fit - check the design and inputs")
    }
    if(any(failed)){
        warning(sum(failed), " of ", length(failed), " null models failed to fit")
    }

    res <- do.call(rbind.data.frame, per.nhood[!failed])
    res$PValue <- pchisq(res$Chisq, df=1, lower.tail=FALSE)

    n.boundary <- length(unique(res$Nhood[res$SigmaBoundary]))
    if(n.boundary > 0){
        message(n.boundary, " of ", length(unique(res$Nhood)),
                " neighbourhoods have a genetic variance component at the boundary.",
                " This is expected when the kinship matrix is close to the identity,",
                " where the genetic variance is confounded with the dispersion.")
    }

    # ---- spatial FDR, computed across neighbourhoods within each variant ---
    res$SpatialFDR <- NA_real_
    for(s in unique(res$SNP)){
        idx <- which(res$SNP == s)
        ord <- idx[order(res$Nhood[idx])]
        sfdr <- graphSpatialFDR(x.nhoods=nhoods(x)[, res$Nhood[ord], drop=FALSE],
                                graph=graph(x),
                                weighting=fdr.weighting,
                                k=x@.k,
                                pvalues=res$PValue[ord],
                                indices=nhoodIndex(x)[res$Nhood[ord]],
                                distances=nhoodDistances(x),
                                reduced.dimensions=reducedDim(x, reduced.dim))
        res$SpatialFDR[ord] <- sfdr
    }

    rownames(res) <- NULL
    res
}


#' Drop random effect terms from a model formula
#'
#' Internal helper returning the fixed effect part of a mixed model formula.
#'
#' @param f A \code{formula} possibly containing \code{lme4}-style random
#' effect terms.
#'
#' @return A \code{formula} containing only the fixed effect terms.
#'
#' @author Mike Morgan
#'
#' @examples
#' miloR:::.dropRandomTerms(~ age + (1 | batch))
#'
#' @name dot-dropRandomTerms
#' @importFrom stats as.formula terms
.dropRandomTerms <- function(f){
    tl <- attr(terms(f), "term.labels")
    keep <- tl[!grepl("\\|", tl)]
    if(length(keep) == 0){
        return(as.formula("~ 1"))
    }
    as.formula(paste("~", paste(keep, collapse=" + ")))
}


#' Expand a random effect design into indicator columns
#'
#' Internal helper turning a matrix of random effect grouping variables into a
#' binary indicator design matrix, recording which columns belong to which
#' random effect.
#'
#' @param z A \code{matrix} or \code{data.frame} of grouping variables, one
#' column per random effect.
#'
#' @return A \code{matrix} of indicator columns with a \code{u_indices}
#' attribute giving the 1-based column indices belonging to each random effect.
#'
#' @author Mike Morgan
#'
#' @examples
#' zz <- matrix(c(1, 1, 2, 2), ncol=1, dimnames=list(NULL, "batch"))
#' miloR:::.expandRandomDesign(zz)
#'
#' @name dot-expandRandomDesign
.expandRandomDesign <- function(z){
    z <- as.matrix(z)
    blocks <- list()
    idx <- list()
    offset <- 0L
    for(j in seq_len(ncol(z))){
        lv <- sort(unique(z[, j]))
        zj <- vapply(lv, FUN=function(l) as.numeric(z[, j] == l),
                     FUN.VALUE=numeric(nrow(z)))
        zj <- matrix(zj, nrow=nrow(z))
        colnames(zj) <- paste0(colnames(z)[j], lv)
        blocks[[j]] <- zj
        idx[[j]] <- seq_len(ncol(zj)) + offset
        offset <- offset + ncol(zj)
    }
    out <- do.call(cbind, blocks)
    names(idx) <- colnames(z)
    attr(out, "u_indices") <- idx
    out
}
