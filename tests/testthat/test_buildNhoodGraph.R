context("Testing buildNhoodGraph function")
library(miloR)

### Set up a mock data set using simulated data
library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(patchwork)

data("sim_trajectory", package = "miloR")

## Extract SingleCellExperiment object
traj_sce <- sim_trajectory[['SCE']]

## Extract sample metadata to use for testing
traj_meta <- sim_trajectory[["meta"]]

## Add metadata to colData slot
colData(traj_sce) <- DataFrame(traj_meta)
logcounts(traj_sce) <- log(counts(traj_sce) + 1)
traj_sce <- runPCA(traj_sce, ncomponents=30)
traj_sce <- runUMAP(traj_sce)
traj_milo <- Milo(traj_sce)
reducedDim(traj_milo, "UMAP") <- reducedDim(traj_sce, "UMAP")

traj_milo <- buildGraph(traj_milo, k = 10, d = 30)

test_that("Expected errors on missing input", {
    expect_error(buildNhoodGraph(traj_milo),
                 "No neighbourhoods found")

    expect_error(buildNhoodGraph(traj_sce),
                 "Not a valid Milo object")
})

test_that("Code produced identical output", {
    library(igraph)
    library(gtools)
    set.seed(42)
    overlap <- 1
    traj_milo <- suppressWarnings(makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE))
    real.graph <- nhoodGraph(buildNhoodGraph(traj_milo, overlap=overlap))

    ## Build adjacency matrix for nhoods
    nhoods <- nhoods(traj_milo)

    ## Make igraph object
    nms <- permutations(n = ncol(nhoods), v = colnames(nhoods), r = 2, repeats.allowed = TRUE)
    # keep_pairs <- sapply( 1:nrow(nms) , function(x) any(nhoods[[nms[x,1]]] %in% nhoods[[ nms[x,2] ]]))
    nh_intersect_mat <- Matrix::crossprod(nhoods)
    nh_intersect_mat[nh_intersect_mat < overlap] <- 0

    message("Calculating nhood adjacency....")
    rownames(nh_intersect_mat) <- colnames(nhoods)
    colnames(nh_intersect_mat) <- colnames(nhoods)

    ig <- graph_from_adjacency_matrix(nh_intersect_mat, mode="undirected", weighted=TRUE)
    nhood_sizes <- sapply(nhoods(traj_milo), length)
    ig <- set_vertex_attr(ig, name = 'size', value = nhood_sizes[vertex.attributes(ig)$name])

    expect_true(is_igraph(real.graph))
    expect_true(is_igraph(ig))

    expect_equal(length(V(ig)), length(V(real.graph)))
    expect_equal(length(E(ig)), length(E(real.graph)))

    expect_equal(length(V(intersection(ig, real.graph))), length(V(real.graph)))
    expect_equal(length(E(difference(ig, real.graph))), 0)
})

test_that("Overlap reduces edge connections", {
    set.seed(42)
    data("sim_discrete", package = "miloR")

    ## Extract SingleCellExperiment object
    traj_sce <- sim_trajectory[['SCE']]

    ## Extract sample metadata to use for testing
    traj_meta <- sim_trajectory[["meta"]]

    ## Add metadata to colData slot
    colData(traj_sce) <- DataFrame(traj_meta)
    logcounts(traj_sce) <- log(counts(traj_sce) + 1)
    traj_sce <- runPCA(traj_sce, ncomponents=30)
    traj_sce <- runUMAP(traj_sce)
    traj_milo <- Milo(traj_sce)
    reducedDim(traj_milo, "UMAP") <- reducedDim(traj_sce, "UMAP")

    traj_milo <- buildGraph(traj_milo, k = 10, d = 30)
    traj_milo <- suppressWarnings(makeNhoods(traj_milo, prop = 0.3, k = 10, d=30, refined = TRUE))

    overlap <- 1
    full.graph <- nhoodGraph(buildNhoodGraph(traj_milo, overlap=overlap))
    sub.graph <- nhoodGraph(buildNhoodGraph(traj_milo, overlap=150))

    expect_true(length(E(full.graph)) > length(E(sub.graph)))
    expect_equal(length(V(full.graph)), length(V(sub.graph)))
})


