# Columns referred to by non-standard evaluation inside ggplot2 aesthetics.
# Declaring them keeps R CMD check from reporting them as undefined globals.
utils::globalVariables(c("Expected", "Observed"))

#' Refit variants passing a screening threshold with the exact model
#'
#' Take the results of \code{\link{testGeneticNhoods}}, select the
#' neighbourhood and variant combinations that pass a lenient screening
#' threshold, and refit the full negative binomial GLMM with the variant in the
#' design to obtain exact Wald statistics.
#'
#' The score test used for the genome-wide scan holds the variance components
#' at the values estimated under the null. That is the standard approximation
#' for genome-wide mixed models and is asymptotically equivalent to the Wald
#' test, but it is not exact for large effects. This function performs the exact
#' refit for the small number of combinations that survive screening, so that
#' reported effect sizes and p-values do not rest on the approximation.
#'
#' @param x A \code{\linkS4class{Milo}} object with a non-empty
#' \code{nhoodCounts} slot.
#' @param results A \code{data.frame} returned by
#' \code{\link{testGeneticNhoods}}.
#' @param design A \code{formula} describing the null model, as passed to
#' \code{\link{testGeneticNhoods}}.
#' @param design.df A \code{data.frame} containing the model variables.
#' @param genotypes A \code{matrix} of genotypes, samples in rows and variants
#' in columns.
#' @param kinship A \code{matrix} genetic relationship matrix.
#' @param screen.threshold A \code{numeric} scalar. Combinations with a
#' \code{PValue} at or below this value are refitted. The default of 1e-4 is
#' deliberately lenient relative to a genome-wide threshold, so that the exact
#' refit rather than the screen decides significance.
#' @param screen.on A \code{character} scalar naming the column used for
#' screening, either \code{"PValue"} or \code{"SpatialFDR"}.
#' @param norm.method A \code{character} scalar giving the normalisation used to
#' compute offsets. One of \code{"TMM"}, \code{"RLE"} or \code{"logMS"}.
#' @param cell.sizes A named \code{numeric} vector of column sizes used for
#' normalisation.
#' @param REML A \code{logical} scalar controlling REML estimation of the
#' variance components.
#' @param max.iters A \code{numeric} scalar giving the maximum number of PQL
#' iterations per refit.
#' @param theta.conv A \code{numeric} scalar giving the convergence tolerance.
#' @param BPPARAM A \code{\linkS4class{BiocParallelParam}} object.
#'
#' @details Each refit is warm-started from the null model estimates for that
#' neighbourhood, which are recomputed once per neighbourhood and shared across
#' every variant screened in it.
#'
#' @return A \code{data.frame} with one row per refitted combination:
#' \describe{
#' \item{\code{Nhood}:}{the neighbourhood index.}
#' \item{\code{SNP}:}{the variant identifier.}
#' \item{\code{logFC}:}{the exact log fold change from the full refit.}
#' \item{\code{SE}:}{the standard error of the log fold change.}
#' \item{\code{tvalue}:}{the Wald statistic.}
#' \item{\code{PValue}:}{the exact Wald p-value.}
#' \item{\code{ScorePValue}:}{the screening p-value from the score test.}
#' \item{\code{ScoreLogFC}:}{the one-step effect estimate from the score test.}
#' \item{\code{Converged}:}{whether the refit met the convergence tolerance.}
#' }
#'
#' @author Mike Morgan
#'
#' @examples
#' NULL
#'
#' @name refineGeneticHits
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom stats model.matrix pchisq terms
#' @export
refineGeneticHits <- function(x, results, design, design.df, genotypes, kinship,
                              screen.threshold=1e-4, screen.on="PValue",
                              norm.method="TMM", cell.sizes=NULL, REML=TRUE,
                              max.iters=100, theta.conv=1e-6,
                              BPPARAM=BiocParallel::SerialParam()){

    if(!is(x, "Milo")){
        stop("Unrecognised input type - must be of class Milo")
    }
    screen.on <- match.arg(screen.on, c("PValue", "SpatialFDR"))
    if(!screen.on %in% colnames(results)){
        stop(screen.on, " not found in results - was this produced by testGeneticNhoods?")
    }

    hits <- results[!is.na(results[[screen.on]]) &
                        results[[screen.on]] <= screen.threshold, , drop=FALSE]
    if(nrow(hits) == 0){
        message("No variants pass the screening threshold of ", screen.threshold)
        return(hits[0, , drop=FALSE])
    }
    message(nrow(hits), " neighbourhood-variant combinations pass the ",
            screen.on, " screen at ", screen.threshold,
            " across ", length(unique(hits$Nhood)), " neighbourhoods")

    genotypes <- as.matrix(genotypes)
    kinship <- as.matrix(kinship)

    keep.samps <- intersect(colnames(nhoodCounts(x)), rownames(design.df))
    keep.samps <- intersect(keep.samps, rownames(genotypes))
    keep.samps <- intersect(keep.samps, rownames(kinship))

    design.df <- design.df[keep.samps, , drop=FALSE]
    genotypes <- genotypes[keep.samps, , drop=FALSE]
    kinship <- kinship[keep.samps, keep.samps, drop=FALSE]
    n.samps <- length(keep.samps)

    has.re <- any(grepl("\\|", attr(terms(design), "term.labels")))
    x.model <- model.matrix(.dropRandomTerms(design), data=design.df)
    rownames(x.model) <- rownames(design.df)

    if(isTRUE(has.re)){
        z.model <- as.matrix(.parse_formula(design, design.df, vtype="re"))
        rownames(z.model) <- rownames(design.df)
        z.full <- .expandRandomDesign(z.model[keep.samps, , drop=FALSE])
        u.indices <- attr(z.full, "u_indices")
    } else {
        z.full <- matrix(0, nrow=n.samps, ncol=0)
        u.indices <- list()
    }
    n.re <- length(u.indices) + 1L

    counts <- nhoodCounts(x)[, keep.samps, drop=FALSE]
    if(is.null(cell.sizes)){
        cell.sizes <- colSums(counts)
    }
    dge <- DGEList(counts=counts, lib.size=cell.sizes)
    if(norm.method %in% c("TMM", "RLE")){
        dge <- calcNormFactors(dge, method=norm.method)
    }
    offsets <- log(dge$samples$lib.size * dge$samples$norm.factors)

    kin.inv <- tryCatch(chol2inv(chol(kinship)), error=function(e) NULL)
    by.nhood <- split(hits, hits$Nhood)

    refit <- bplapply(names(by.nhood), BPPARAM=BPPARAM, FUN=function(nm){
        hn <- by.nhood[[nm]]
        nh <- as.integer(nm)
        yi <- as.numeric(counts[nh, ])
        mu0 <- rep(max(mean(yi), 1e-2), n.samps)
        b0 <- c(log(max(mean(yi), 1e-2)) - mean(offsets), rep(0, ncol(x.model) - 1))

        # one null fit per neighbourhood, shared by every variant screened in it
        nullfit <- tryCatch(
            fitGeneticNullGlmm(Z=z.full, X=x.model, K=kinship, muvec=mu0,
                               offsets=offsets, curr_beta=b0,
                               curr_u=rep(0, ncol(z.full) + n.samps),
                               curr_sigma=rep(0.5, n.re), y=yi,
                               u_indices=u.indices, theta_conv=theta.conv,
                               curr_disp=1, REML=REML, maxit=max.iters,
                               Kinv_=kin.inv, return_projection=FALSE),
            error=function(e) NULL)

        null.sigma <- if(is.null(nullfit)) NULL else as.numeric(nullfit$Sigma)
        null.disp <- if(is.null(nullfit)) -1 else nullfit$Dispersion
        null.beta <- if(is.null(nullfit)) NULL else as.numeric(nullfit$FE)

        rows <- lapply(seq_len(nrow(hn)), FUN=function(j){
            snp <- hn$SNP[j]
            xa <- cbind(x.model, SNP=genotypes[, snp])
            alt <- tryCatch(
                fitGeneticNullGlmm(Z=z.full, X=xa, K=kinship, muvec=mu0,
                                   offsets=offsets, curr_beta=c(b0, 0),
                                   curr_u=rep(0, ncol(z.full) + n.samps),
                                   curr_sigma=rep(0.5, n.re), y=yi,
                                   u_indices=u.indices, theta_conv=theta.conv,
                                   curr_disp=1, REML=REML, maxit=max.iters,
                                   Kinv_=kin.inv, null_sigma_=null.sigma,
                                   null_beta_=null.beta, null_disp=null.disp,
                                   return_projection=FALSE),
                error=function(e) NULL)
            if(is.null(alt)){
                return(NULL)
            }
            k <- ncol(xa)
            tv <- as.numeric(alt$t)[k]
            data.frame(Nhood=nh, SNP=snp,
                       logFC=as.numeric(alt$FE)[k],
                       SE=as.numeric(alt$SE)[k],
                       tvalue=tv,
                       PValue=pchisq(tv^2, df=1, lower.tail=FALSE),
                       ScorePValue=hn$PValue[j],
                       ScoreLogFC=hn$logFC[j],
                       Converged=as.logical(alt$converged),
                       stringsAsFactors=FALSE)
        })
        do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    })

    out <- do.call(rbind, refit[!vapply(refit, is.null, logical(1))])
    if(is.null(out) || nrow(out) == 0){
        warning("All exact refits failed")
        return(hits[0, , drop=FALSE])
    }
    rownames(out) <- NULL
    out[order(out$PValue), ]
}


#' Two-panel QQ plot with genomic inflation for genetic neighbourhood tests
#'
#' Produce a two-panel quantile-quantile plot for the output of
#' \code{\link{testGeneticNhoods}}. The left panel shows the raw p-values and
#' the right panel the spatial FDR corrected values, both against the uniform
#' expectation on the -log10 scale.
#'
#' @param results A \code{data.frame} returned by
#' \code{\link{testGeneticNhoods}}, containing \code{PValue} and
#' \code{SpatialFDR} columns.
#' @param facet.by An optional \code{character} scalar naming a column of
#' \code{results} to split the plot by, for example \code{"SNP"}.
#' @param title An optional \code{character} scalar giving the plot title.
#' @param point.size A \code{numeric} scalar giving the point size.
#'
#' @details The genomic inflation factor reported for the left panel is the
#' usual \eqn{\lambda_{GC}}, the median observed chi-squared statistic divided
#' by the median of the null chi-squared distribution on one degree of freedom.
#' A value near 1 indicates calibrated test statistics; values above 1 indicate
#' inflation, typically from uncontrolled population structure.
#'
#' The same quantity is reported for the right panel for comparability, but note
#' it is not a genomic inflation factor in the usual sense: spatial FDR values
#' are adjusted p-values and are not uniformly distributed under the null even
#' when the raw p-values are. The right panel should be read as showing the
#' effect of the correction rather than as a calibration diagnostic in its own
#' right.
#'
#' @return A \code{ggplot} object. The genomic inflation factors are attached as
#' the \code{lambda} attribute, a named \code{numeric} vector with elements
#' \code{PValue} and \code{SpatialFDR}.
#'
#' @author Mike Morgan
#'
#' @examples
#' res <- data.frame(Nhood=rep(seq_len(50), 2),
#'                   SNP=rep(c("rs1", "rs2"), each=50),
#'                   PValue=runif(100),
#'                   SpatialFDR=runif(100))
#' pl <- plotGeneticQQ(res)
#' attr(pl, "lambda")
#'
#' @name plotGeneticQQ
#' @importFrom ggplot2 ggplot aes geom_point geom_abline facet_wrap labs theme_bw
#' @importFrom ggplot2 theme element_text

#' @importFrom stats qchisq median ppoints
#' @export
plotGeneticQQ <- function(results, facet.by=NULL, title=NULL, point.size=0.8){

    if(!all(c("PValue", "SpatialFDR") %in% colnames(results))){
        stop("results must contain PValue and SpatialFDR columns - was this produced by testGeneticNhoods?")
    }

    .lambda <- function(p){
        p <- p[!is.na(p) & p > 0 & p <= 1]
        if(length(p) == 0){
            return(NA_real_)
        }
        stats::median(stats::qchisq(p, df=1, lower.tail=FALSE)) /
            stats::qchisq(0.5, df=1, lower.tail=FALSE)
    }

    .qqframe <- function(p, panel){
        keep <- !is.na(p) & p > 0 & p <= 1
        p <- sort(p[keep])
        if(length(p) == 0){
            return(NULL)
        }
        n <- length(p)
        data.frame(Expected=-log10(stats::ppoints(n)),
                   Observed=-log10(p),
                   Panel=panel,
                   stringsAsFactors=FALSE)
    }

    lam <- c(PValue=.lambda(results$PValue), SpatialFDR=.lambda(results$SpatialFDR))

    if(!is.null(facet.by)){
        if(!facet.by %in% colnames(results)){
            stop(facet.by, " not found in results")
        }
        grps <- split(results, results[[facet.by]])
        qq <- do.call(rbind, lapply(names(grps), FUN=function(g){
            gd <- grps[[g]]
            fr <- rbind(.qqframe(gd$PValue, "Raw p-value"),
                        .qqframe(gd$SpatialFDR, "Spatial FDR"))
            if(is.null(fr)){
                return(NULL)
            }
            fr$Group <- g
            fr
        }))
    } else {
        qq <- rbind(.qqframe(results$PValue, "Raw p-value"),
                    .qqframe(results$SpatialFDR, "Spatial FDR"))
    }

    if(is.null(qq) || nrow(qq) == 0){
        stop("No usable p-values to plot")
    }

    qq$Panel <- factor(qq$Panel, levels=c("Raw p-value", "Spatial FDR"))
    lab.raw <- sprintf("Raw p-value  (lambda = %.3f)", lam[["PValue"]])
    lab.fdr <- sprintf("Spatial FDR  (lambda = %.3f)", lam[["SpatialFDR"]])
    levels(qq$Panel) <- c(lab.raw, lab.fdr)

    if(is.null(title)){
        title <- "Neighbourhood csQTL quantile-quantile plot"
    }

    pl <- ggplot(qq, aes(x=Expected, y=Observed)) +
        geom_abline(slope=1, intercept=0, colour="grey50", linetype="dashed") +
        geom_point(size=point.size) +
        labs(x=expression(Expected~-log[10](p)),
             y=expression(Observed~-log[10](p)),
             title=title) +
        theme_bw() +
        theme(strip.text=element_text(size=10))

    if(!is.null(facet.by)){
        pl <- pl + facet_wrap(~ Panel + Group, scales="free")
    } else {
        pl <- pl + facet_wrap(~ Panel, scales="free")
    }

    attr(pl, "lambda") <- lam
    pl
}
