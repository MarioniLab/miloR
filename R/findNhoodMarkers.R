#' Identify post-hoc neighbourhood marker genes
#'
#' This function will perform differential gene expression analysis on
#' differentially abundant neighbourhoods, by first aggregating adjacent and
#' concordantly DA neighbourhoods, then comparing cells \emph{between} these
#' aggregated groups. For differential gene experession based on an input design
#' \emph{within} DA neighbourhoods see \code{\link{testDiffExp}}.
#'
#' @param x A \code{\linkS4class{Milo}} object containing single-cell gene expression
#' and neighbourhoods.
#' @param da.res A \code{data.frame} containing DA results, as expected from running
#' \code{testNhoods}.
#' @param da.fdr A numeric scalar that determines at what FDR neighbourhoods are declared
#' DA for the purposes of aggregating across concorantly DA neighbourhoods.
#' @param assay A character scalar determining which \code{assays} slot to extract from the
#' \code{\linkS4class{Milo}} object to use for DGE testing.
#' @param overlap A scalar integer that determines the number of cells that must
#' overlap between adjacent neighbourhoods for merging.
#' @param lfc.threshold A scalar that determines the absolute log fold change above
#' which neighbourhoods should be considerd 'DA' for merging. Default=NULL
#' @param merge.discord A logical scalar that overrides the default behaviour and allows
#' adjacent neighbourhoods to be merged if they have discordant log fold change signs. Using
#' this argument is generally discouraged, but may be useful for constructing an empirical null
#' group of cells, regardless of DA sign.
#' @param subset.row A logical, integer or character vector indicating the rows
#' of \code{x} to use for sumamrizing over cells in neighbourhoods.
#' @param gene.offset A logical scalar the determines whether a per-cell offset
#' is provided in the DGE GLM to adjust for the number of detected genes with
#' expression > 0.
#' @param return.groups A logical scalar that returns a \code{\link{data.frame}} of the
#' aggregated groups per single-cell. Cells that are members of non-DA neighbourhoods contain
#' \code{NA} values.
#' @param subset.nhoods A logical, integer or character vector indicating which neighbourhoods
#' to subset before aggregation and DGE testing.
#' @param na.function A valid NA action function to apply, should be one of
#' \code{na.fail, na.omit, na.exclude, na.pass}.
#'
#' @details
#' Adjacent neighbourhoods are first merged based on two criteria: 1) they share at
#' least \code{overlap} number of cells, and 2) the DA log fold change sign is concordant.
#' This behaviour can be modulated by setting \code{overlap} to be more or less stringent.
#' Additionally, a threshold on the log fold-changes can be set, such that \code{lfc.threshold}
#' is required to merge adjacent neighbourhoods. Note: adjacent neighbourhoods will never be
#' merged with opposite signs.
#'
#' Using a one vs. all approach, each aggregated group of cells is compared to all others
#' using the single-cell log normalized gene expression with a GLM
#' (for details see \code{\link[limma]{limma-package}}), or the single-cell counts using a
#' negative binomial GLM (for details see \code{\link[edgeR]{edgeR-package}}). When using
#' the latter it is recommended to set \code{gene.offset=TRUE} as this behaviour adjusts
#' the model offsets by the number of detected genes in each cell.
#'
#'
#' @return A \code{data.frame} of DGE results containing a log fold change and adjusted
#' p-value for each aggregated group of neighbourhoods. If \code{return.groups} then
#' the return value is a list with the slots \code{groups} and \code{dge} containing the
#' aggregated neighbourhood groups per single-cell and marker gene results, respectively.
#'
#' @author Mike Morgan & Emma Dann
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
#' colnames(sce) <- paste0("Cell", 1:ncol(sce))
#' milo <- Milo(sce)
#' milo <- buildGraph(milo, k=20, d=10, transposed=TRUE)
#' milo <- makeNhoods(milo, k=20, d=10, prop=0.3)
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
#' da.res <- testNhoods(milo, design=~Condition, design.df=test.meta[colnames(nhoodCounts(milo)), ])
#'
#' nhood.dge <- findNhoodMarkers(milo, da.res, overlap=15)
#' nhood.dge
#'
#' @name findNhoodMarkers
NULL


#' @export
#' @importFrom stats model.matrix as.formula
#' @importFrom Matrix colSums
findNhoodMarkers <- function(x, da.res, da.fdr=0.1, assay="logcounts",
                             overlap=1, lfc.threshold=NULL, merge.discord=FALSE,
                             subset.row=NULL, gene.offset=TRUE,
                             return.groups=FALSE, subset.nhoods=NULL,
                             na.function="na.pass"){

    if(class(x) != "Milo"){
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

    message(paste0("Found ", n.da, " DA neighbourhoods at FDR ", da.fdr*100, "%"))

    nhs.da.gr <- .group_nhoods_by_overlap(nhoods(x),
                                          da.res=da.res,
                                          is.da=da.res$SpatialFDR < da.fdr,
                                          overlap=overlap,
                                          subset.nhoods=subset.nhoods) # returns a vector group values for each nhood
    nhood.gr <- unique(nhs.da.gr)
    # perform DGE _within_ each group of cells using the input design matrix
    message(paste0("Nhoods aggregated into ", length(nhood.gr), " groups"))

    fake.meta <- data.frame("CellID"=colnames(x), "Nhood.Group"=rep(NA, ncol(x)))
    rownames(fake.meta) <- fake.meta$CellID

    for(i in seq_along(nhood.gr)){
        nhood.x <- nhs.da.gr == nhood.gr[i]
        fake.meta[unlist(nhoods(x)[da.res$SpatialFDR < da.fdr][nhood.x]),]$Nhood.Group <- nhood.gr[i]
    }

    # only compare against the other DA neighbourhoods
    x <- x[, !is.na(fake.meta$Nhood.Group)]
    fake.meta <- fake.meta[!is.na(fake.meta$Nhood.Group), ]

    if(!is.null(subset.row)){
        x <- x[subset.row, , drop=FALSE]
    }

    exprs <- assay(x, assay)

    marker.list <- list()
    i.contrast <- c("TestTest - TestRef") # always use contrasts for this

    for(i in seq_along(nhood.gr)){
        i.meta <- fake.meta
        i.meta$Test <- "Ref"
        i.meta$Test[fake.meta$Nhood.Group == nhood.gr[i]] <- "Test"

        if(ncol(exprs) > 1 & nrow(i.meta) > 1){
            i.design <- as.formula(" ~ 0 + Test")
            i.model <- model.matrix(i.design, data=i.meta)
            rownames(i.model) <- rownames(i.meta)
        }

        sink(file="/dev/null")
        gc()
        sink(file=NULL)

        if(assay == "logcounts"){
            i.res <- .perform_lognormal_dge(exprs, i.model, model.contrasts=i.contrast,
                                            gene.offset=gene.offset)
        } else if(assay == "counts"){
            i.res <- .perform_counts_dge(exprs, i.model, model.contrasts=i.contrast,
                                         gene.offset=gene.offset)
        } else{
            warning("Assay type is not counts or logcounts - assuming (log)-normal distribution. Use these results at your peril")
            i.res <- .perform_lognormal_dge(exprs, i.model,
                                            model.contrasts=i.contrast,
                                            gene.offset=gene.offset)
        }

        i.res$adj.P.Val[is.na(i.res$adj.P.Val)] <- 1
        i.res$logFC[is.infinite(i.res$logFC)] <- 0

        i.res <- i.res[, c("logFC", "adj.P.Val")]
        colnames(i.res) <- paste(colnames(i.res), nhood.gr[i], sep="_")
        marker.list[[paste0(nhood.gr[i])]] <- i.res
    }

    marker.df <- do.call(cbind.data.frame, marker.list)
    colnames(marker.df) <- gsub(colnames(marker.df), pattern="^[0-9]+\\.", replacement="")
    marker.df$GeneID <- rownames(i.res)

    if(isTRUE(return.groups)){
        out.list <- list("groups"=fake.meta, "dge"=marker.df)
        return(out.list)
    }else{
        return(marker.df)
    }
}






