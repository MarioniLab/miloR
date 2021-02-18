#' Perform post-hoc differential gene expression analysis
#'
#' This function will perform differential gene expression analysis with
#' differentially abundant neighbourhoods, by first aggregating adjacent and
#' concordantly DA neighbourhoods, then comparing cells \emph{within} these
#' aggregated groups for differential gene experession using the input design. For
#' comparing \emph{between} DA neighbourhoods see \code{\link{findNhoodMarkers}}.
#'
#' @param x A \code{\linkS4class{Milo}} object containing single-cell gene expression
#' and neighbourhoods.
#' @param da.res A \code{data.frame} containing DA results, as expected from running
#' \code{testNhoods}.
#' @param design A \code{formula} or \code{model.matrix} object describing the
#' experimental design for differential gene expression testing. The last component
#' of the formula or last column of the model matrix are by default the test
#' variable. This behaviour can be overridden by setting the \code{model.contrasts}
#' argument. This should be the same as was used for DA testing.
#' @param meta.data A cell X variable \code{data.frame} containing single-cell meta-data
#' to which \code{design} refers. The order of rows (cells) must be the same as the
#' \code{\linkS4class{Milo}} object columns.
#' @param model.contrasts A string vector that defines the contrasts used to perform
#' DA testing. This should be the same as was used for DA testing.
#' @param assay A character scalar determining which \code{assays} slot to extract from the
#' \code{\linkS4class{Milo}} object to use for DGE testing.
#' @param subset.row A logical, integer or character vector indicating the rows
#' of \code{x} to use for sumamrizing over cells in neighbourhoods.
#' @param subset.nhoods A logical, integer or character vector indicating which neighbourhoods
#' to subset before aggregation and DGE testing (default: NULL).
#' @param gene.offset A logical scalar the determines whether a per-cell offset
#' is provided in the DGE GLM to adjust for the number of detected genes with
#' expression > 0.
#' @param n.coef A numeric scalar refering to the coefficient to select from the
#' DGE model. This is especially pertinent when passing an ordered variable and
#' only one specific type of effects are to be tested.
#' @param na.function A valid NA action function to apply, should be one of
#' \code{na.fail, na.omit, na.exclude, na.pass}.
#'
#' @details
#' Adjacent neighbourhoods are first merged based on two criteria: 1) they share at
#' least \code{overlap} number of cells, and 2) the DA log fold change sign is concordant.
#' This behaviour can be modulated by setting \code{overlap} to be more or less stringent.
#' Additionally, a threshold on the log fold-changes can be set, such that \code{lfc.threshold}
#' is required to merge adjacent neighbourhoods. Note: adjacent neighbourhoods will never be
#' merged with opposite signs unless \code{merge.discord=TRUE}.
#'
#' Within each aggregated group of cells differential gene expression testing is performed
#' using the single-cell log normalized gene expression with a GLM
#' (for details see \code{\link[limma]{limma-package}}), or the single-cell counts using a
#' negative binomial GLM (for details see \code{\link[edgeR]{edgeR-package}}). When using
#' single-cell data for DGE it is recommended to set \code{gene.offset=TRUE} as this
#' behaviour adjusts the model by the number of detected genes in each cell as a proxy for
#' differences in capture efficiency and cellular RNA content.
#'
#'
#' @return A \code{list} containing a \code{data.frame} of DGE results for each aggregated
#' group of neighbourhoods.
#'
#' @author Mike Morgan & Emma Dann
#'
#' @examples
#' data(sim_discrete)
#'
#' milo <- Milo(sim_discrete$SCE)
#' milo <- buildGraph(milo, k=20, d=10, transposed=TRUE)
#' milo <- makeNhoods(milo, k=20, d=10, prop=0.3)
#'
#' meta.df <- sim_discrete$meta
#' meta.df$SampID <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
#' milo <- countCells(milo, meta.data=meta.df, samples="SampID")
#'
#' test.meta <- data.frame("Condition"=c(rep("A", 3), rep("B", 3)), "Replicate"=rep(c("R1", "R2", "R3"), 2))
#' test.meta$Sample <- paste(test.meta$Condition, test.meta$Replicate, sep="_")
#' rownames(test.meta) <- test.meta$Sample
#' da.res <- testNhoods(milo, design=~Condition, design.df=test.meta[colnames(nhoodCounts(milo)), ])
#'
#' nhood.dge <- testDiffExp(milo, da.res, da.fdr=0.2, design=~Condition, meta.data=meta.df, overlap=1, compute.new=TRUE)
#' nhood.dge
#'
#' @name testDiffExp
#' @export
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment assayNames
#' @importFrom stats na.pass na.fail na.omit na.exclude
testDiffExp <- function(x, da.res, design, meta.data, model.contrasts=NULL,
                        assay="logcounts",
                        subset.nhoods=NULL,
                        subset.row=NULL, gene.offset=TRUE, n.coef=NULL,
                        na.function="na.pass"
                        ){

    if(!is(x, "Milo")){
        stop("Unrecognised input type - must be of class Milo")
    } else if(any(!assay %in% assayNames(x))){
        stop(paste0("Unrecognised assay slot: ", assay))
    }

    if(is.null(na.function)){
        warning("NULL passed to na.function, using na.pass")
        na.func <- get("na.pass")
    } else{
        tryCatch({
            na.func <- get(na.function)
        }, warning=function(warn){
            warning(warn)
        }, error=function(err){
            stop(paste0("NA function ", na.function, " not recognised"))
        }, finally={
        })
    }

    n.da <- sum(na.func(da.res$SpatialFDR < da.fdr))
    if(!is.na(n.da) & n.da == 0){
        stop("No DA neighbourhoods found")
    }

    if(any(is.na(da.res$SpatialFDR))){
        warning("NA values found in SpatialFDR vector")
    }

    # assign this group level information to the consituent cells using the input meta.data
    copy.meta <- meta.data # make a copy assume the order is the same, just check the rownames are the same

    if(!all(rownames(copy.meta) == colnames(x))){
        warning("Column names of x are not the same as meta-data rownames")
    }
    
    copy.meta$Nhood.Group <- NA
    nhs.da.gr <- da.res$NhoodGroup
    names(nhs.da.gr) <- da.res$Nhood
    
    if(!is.null(subset.nhoods)){
        nhs.da.gr <- nhs.da.gr[subset.nhoods]
    }
    
    nhood.gr <- unique(nhs.da.gr)

    for(i in seq_along(nhood.gr)){
        nhood.x <- nhs.da.gr %in% nhood.gr[i]
        # get the nhoods
        nhs <- nhoods(x)
        if(!is.null(subset.nhoods)){
            nhs <- nhs[,subset.nhoods]
        }
        nhood.gr.cells <- rowSums(nhs[, nhood.x, drop=FALSE]) > 0
        copy.meta[nhood.gr.cells,]$Nhood.Group <- nhood.gr[i]
    }

    # subset to non-NA group cells
    subset.dims <- !is.na(copy.meta$Nhood.Group)
    x <- x[, subset.dims]
    copy.meta <- copy.meta[subset.dims, ]

    if(is(design, "formula")){
        model <- model.matrix(design, data=copy.meta)
        rownames(model) <- rownames(copy.meta)
    } else if(is.matrix(design)){
        model <- design
        if(nrow(model) != nrow(copy.meta)){
            message("Subsetting input design matrix to DA neighbourhood cells")
            if(length(subset.dims) == nrow(model)){
                model <- model[subset.dims, ]
            } else{
                stop(paste0("Cannot subset model matrix, subsetting vector is wrong length:", length(subset.dims)))
            }
        }
        if(any(rownames(model) != rownames(copy.meta))){
            warning("Design matrix and meta-data dimnames are not the same")
        }
    }

    if(ncol(x) != nrow(model)){
        stop(paste0("Design matrix (", nrow(model), ") and milo objects (",
                    ncol(x), ") are not the same dimension"))
    }

    if(!is.null(subset.row)){
        x <- x[subset.row, , drop=FALSE]
    }

    sink(file="/dev/null")
    gc()
    sink(file=NULL)

    # perform DGE _within_ each group of cells using the input design matrix
    dge.list <- list()
    for(i in seq_along(nhood.gr)){
        i.meta <- copy.meta[copy.meta$Nhood.Group == nhood.gr[i], ,drop=FALSE]
        i.exprs <- assay(x[, copy.meta$Nhood.Group == nhood.gr[i],drop=FALSE], assay)
        i.model <- model[copy.meta$Nhood.Group == nhood.gr[i], ,drop=FALSE]

        if(!is.null(ncol(i.exprs))){
            if(ncol(i.exprs) > (ncol(i.model) + 1) & nrow(i.meta) > (ncol(i.model) + 1)){
                if(assay == "logcounts"){
                    i.res <- .perform_lognormal_dge(i.exprs, i.model,
                                                    model.contrasts=model.contrasts,
                                                    gene.offset=gene.offset, n.coef=n.coef)
                } else if(assay == "counts"){
                    i.res <- .perform_counts_dge(i.exprs, i.model,
                                                 model.contrasts=model.contrasts,
                                                 gene.offset=gene.offset, n.coef=n.coef)
                } else{
                    warning("Assay type is not counts or logcounts - assuming (log)-normal distribution. Use these results at your peril")
                    i.res <- .perform_lognormal_dge(i.exprs, i.model,
                                                    model.contrasts=model.contrasts,
                                                    gene.offset=gene.offset, n.coef=n.coef)
                }
                i.res$adj.P.Val[is.na(i.res$adj.P.Val)] <- 1
                i.res$logFC[is.infinite(i.res$logFC)] <- 0
                i.res$Nhood.Group <- nhood.gr[i]

                dge.list[[paste0(nhood.gr[i])]] <- i.res
        }else if(ncol(i.exprs) <= ncol(i.model)){
            warning("Not enough cells to perform DE testing in this neighbourhood")
        }
    } else{
        warning("Not enough cells to perform DE testing in this neighbourhood")
    }
    }

    return(dge.list)
 }


########################################
## utility functions for testDiffExp ###
########################################

#' @importFrom limma makeContrasts lmFit topTreat eBayes contrasts.fit
.perform_lognormal_dge <- function(exprs.data, test.model, gene.offset=gene.offset,
                                   model.contrasts=NULL, n.coef=NULL){

    if(isTRUE(gene.offset)){
        n.gene <- apply(exprs.data, 2, function(X) sum(X > 0))
        old.col <- colnames(test.model)
        if(all(test.model[, 1] == 1)){
            test.model <- cbind(test.model[, 1], n.gene, test.model[, c(2:ncol(test.model))])
            colnames(test.model) <- c(old.col[1], "NGenes", old.col[c(2:length(old.col))])
        } else{
            test.model <- cbind(n.gene, test.model)
            colnames(test.model) <- c("NGenes", old.col)
        }
    }

    i.fit <- lmFit(exprs.data, test.model)
    if(!is.null(model.contrasts)){
        mod.constrast <- makeContrasts(contrasts=model.contrasts, levels=test.model)
        i.fit <- contrasts.fit(i.fit, contrasts=mod.constrast)
        i.fit <- eBayes(i.fit, trend=TRUE)
        i.res <- as.data.frame(topTreat(i.fit, number = Inf, sort.by = "p", p.value = 1))
    } else{
        i.fit <- eBayes(i.fit, trend=TRUE)
        if(is.null(n.coef)){
            n.coef <- ncol(test.model)
        }

        i.res <- as.data.frame(topTreat(i.fit, coef=ncol(test.model), number = Inf, sort.by = "p", p.value = 1))
    }
    return(i.res)
}


#' @importFrom limma makeContrasts
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags
.perform_counts_dge <- function(exprs.data, test.model, gene.offset=gene.offset,
                                model.contrasts=NULL, n.coef=NULL){

    i.dge <- DGEList(counts=exprs.data,
                     lib.size=log(colSums(exprs.data)))

    if(isTRUE(gene.offset)){
        n.gene <- apply(exprs.data, 2, function(X) sum(X > 0))
        if(ncol(test.model) == 2){
            test.model <- cbind(test.model, n.gene)
            colnames(test.model) <- c(colnames(test.model)[1:2], "NGenes")
        } else if (ncol(test.model) > 2){
            test.model <- cbind(test.model[, 1], n.gene, test.model[, c(2:ncol(test.model))])
            colnames(test.model) <- c(colnames(test.model)[1], "NGenes", colnames(test.model[, c(2:ncol(test.model))]))
        } else{
            if(ncol(test.model) < 2){
                warning("Only one column in model matrix - must have at least 2. gene.offset forced to  FALSE")
                gene.offset <- FALSE
            }
        }
    }

    i.dge <- estimateDisp(i.dge, test.model)
    i.fit <- glmQLFit(i.dge, test.model, robust=TRUE)

    if(!is.null(model.contrasts)){
        mod.constrast <- makeContrasts(contrasts=model.contrasts, levels=test.model)
        i.res <- as.data.frame(topTags(glmQLFTest(i.fit, contrast=mod.constrast),
                                       sort.by='none', n=Inf))
    } else{
        if(is.null(n.coef)){
            n.coef <- ncol(test.model)
        }
        i.res <- as.data.frame(topTags(glmQLFTest(i.fit, coef=n.coef), sort.by='none', n=Inf))
    }
    return(i.res)
}
