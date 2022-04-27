context("Testing testNhoods function")
library(miloR)

### Set up a mock data set using simulated data
library(SingleCellExperiment)
library(scran)
library(scater)
library(irlba)
library(MASS)
library(mvtnorm)

set.seed(42)
r.n <- 1000
n.dim <- 50
block1.cells <- 500
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block1.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block1.eigens <- block1.eigens[order(block1.eigens)]
block1.p <- qr.Q(qr(matrix(rnorm(block1.cells^2, mean=4, sd=0.01), block1.cells)))
block1.sigma <- crossprod(block1.p, block1.p*block1.eigens)
block1.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block1.cells, mean=2, sd=0.01), sigma=block1.sigma))


block2.cells <- 500
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block2.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block2.eigens <- block2.eigens[order(block2.eigens)]
block2.p <- qr.Q(qr(matrix(rnorm(block2.cells^2, mean=4, sd=0.01), block2.cells)))
block2.sigma <- crossprod(block2.p, block2.p*block2.eigens)
block2.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block2.cells, mean=4, sd=0.01), sigma=block2.sigma))


block3.cells <- 650
# select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
block3.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
block3.eigens <- block3.eigens[order(block3.eigens)]
block3.p <- qr.Q(qr(matrix(rnorm(block3.cells^2, mean=4, sd=0.01), block3.cells)))
block3.sigma <- crossprod(block3.p, block3.p*block3.eigens)
block3.gex <- abs(rmvnorm(n=r.n, mean=rnorm(n=block3.cells, mean=5, sd=0.01), sigma=block3.sigma))

sim1.gex <- do.call(cbind, list("b1"=block1.gex, "b2"=block2.gex, "b3"=block3.gex))
colnames(sim1.gex) <- paste0("Cell", 1:ncol(sim1.gex))
rownames(sim1.gex) <- paste0("Gene", 1:nrow(sim1.gex))
sim1.pca <- prcomp_irlba(t(sim1.gex), n=50, scale.=TRUE, center=TRUE)

set.seed(42)
block1.cond <- rep("A", block1.cells)
block1.a <- sample(1:block1.cells, size=floor(block1.cells*0.9))
block1.b <- setdiff(1:block1.cells, block1.a)
block1.cond[block1.b] <- "B"

block2.cond <- rep("A", block2.cells)
block2.a <- sample(1:block2.cells, size=floor(block2.cells*0.05))
block2.b <- setdiff(1:block2.cells, block2.a)
block2.cond[block2.b] <- "B"

block3.cond <- rep("A", block3.cells)
block3.a <- sample(1:block3.cells, size=floor(block3.cells*0.5))
block3.b <- setdiff(1:block3.cells, block3.a)
block3.cond[block3.b] <- "B"

meta.df <- data.frame("Block"=c(rep("B1", block1.cells), rep("B2", block2.cells), rep("B3", block3.cells)),
                      "Condition"=c(block1.cond, block2.cond, block3.cond),
                      "Replicate"=c(rep("R1", floor(block1.cells*0.33)), rep("R2", floor(block1.cells*0.33)),
                                    rep("R3", block1.cells-(2*floor(block1.cells*0.33))),
                                    rep("R1", floor(block2.cells*0.33)), rep("R2", floor(block2.cells*0.33)),
                                    rep("R3", block2.cells-(2*floor(block2.cells*0.33))),
                                    rep("R1", floor(block3.cells*0.33)), rep("R2", floor(block3.cells*0.33)),
                                    rep("R3", block3.cells-(2*floor(block3.cells*0.33)))))
colnames(meta.df) <- c("Block", "Condition", "Replicate")
# define a "sample" as teh combination of condition and replicate
meta.df$Sample <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
meta.df$Vertex <- c(1:nrow(meta.df))

sim1.sce <- SingleCellExperiment(assays=list(logcounts=sim1.gex),
                                 reducedDims=list("PCA"=sim1.pca$x))

sim1.mylo <- Milo(sim1.sce)

# build a graph - this can take a while for large graphs - will need to play
# around with the parallelisation options
sim1.mylo <- buildGraph(sim1.mylo, k=21, d=30)

# define neighbourhoods - this is slow for large data sets
# how can this be sped up? There are probably some parallelisable steps
sim1.mylo <- makeNhoods(sim1.mylo, k=21, prop=0.1, refined=TRUE,
                                d=30,
                                reduced_dims="PCA")
sim1.mylo <- calcNhoodDistance(sim1.mylo, d=30)

sim1.meta <- data.frame("Condition"=c(rep("A", 3), rep("B", 3)),
                        "Replicate"=rep(c("R1", "R2", "R3"), 2))
sim1.meta$Sample <- paste(sim1.meta$Condition, sim1.meta$Replicate, sep="_")
rownames(sim1.meta) <- sim1.meta$Sample

test_that("Wrong input gives errors", {
    # missing neighbourhood counts
    expect_error(testNhoods(sim1.mylo, design=~Condition,
                                    design.df=sim1.meta),
                 "Neighbourhood counts missing - please run countCells first")

    # count cells
    sim1.mylo <- countCells(sim1.mylo, samples="Sample", meta.data=meta.df)
    expect_error(testNhoods(nhoodCounts(sim1.mylo), design=~Condition,
                            design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]),
                 "Unrecognised input type - must be of class Milo")

    # missing reduced dimension slot
    expect_error(testNhoods(sim1.mylo, design=~Condition,
                            design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ],
                            reduced.dim="blah"),
                 "is not found in reducedDimNames")

})

sim1.mylo <- countCells(sim1.mylo, samples="Sample", meta.data=meta.df)

test_that("Discordant dimension names gives an error", {
    design.matrix <- model.matrix(~Condition, data=sim1.meta)
    rownames(design.matrix) <- c(1:nrow(design.matrix))
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                             design.df=sim1.meta)),
                                  "Design matrix and model matrix rownames are not a subset")

    # warning if only a subset are present, or order is wrong
    design.matrix <- model.matrix(~Condition, data=sim1.meta)
    rownames(design.matrix) <- rownames(sim1.meta)

    set.seed(42)
    design.matrix <- design.matrix[sample(rownames(design.matrix)), ]
    expect_warning(testNhoods(sim1.mylo, design=design.matrix,
                              design.df=sim1.meta),
                   "Sample names in design matrix and nhood counts are not matched. Reordering")

})

test_that("Discordant dimensions between input and design gives an error", {
    add.meta <- sim1.meta[c(1:5), ]
    rownames(add.meta) <- paste0(rownames(add.meta), "_add")
    big.meta <- rbind.data.frame(sim1.meta, add.meta)
    design.matrix <- model.matrix(~Condition, data=big.meta)

    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                                     design.df=sim1.meta)),
                 "Design matrix and model matrix are not the same dimensionality")
})


test_that("Concordant dimensions between input and output", {
    in.rows <- nrow(nhoodCounts(sim1.mylo))
    out.rows <- nrow(testNhoods(sim1.mylo, design=~Condition,
                                        design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    expect_identical(in.rows, out.rows)
})


test_that("Identical results are produced with identical input", {
    # run for each weighting scheme
    kd.ref1 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    kd.ref2 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(kd.ref1, kd.ref2)

    # mean nhood distance
    nd.ref1 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="neighbour-distance",
                          design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    nd.ref2 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="neighbour-distance",
                          design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(nd.ref1, nd.ref2)
})

test_that("testNhoods produces reproducible results with equivalent input", {
    # same input, but different format should give the same results
    design.matrix <- model.matrix(~Condition, data=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    mod.ref <- suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                                   design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    form.ref <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition,
                                                    design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    # test for each column of output
    expect_equal(mod.ref$Pvalue, form.ref$Pvalue)
    expect_equal(mod.ref$SpatialFDR, form.ref$SpatialFDR)
    expect_equal(mod.ref$logFC, form.ref$logFC)
    expect_equal(mod.ref$`F`, form.ref$`F`)
    expect_equal(mod.ref$FDR, form.ref$FDR)
})

test_that("Filtering nhoods provides reproducible results", {
    require(Matrix)
    exp.nh <- sum(Matrix::rowMeans(nhoodCounts(sim1.mylo)) >= 5)
    out.da <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                 min.mean=5,
                                 design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(nrow(out.da), exp.nh)

    kd.ref1 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  min.mean=5,
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    kd.ref2 <- testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                  min.mean=5,
                                  design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    expect_identical(kd.ref1, kd.ref2)
})

test_that("Model contrasts provide expected results", {
    design.matrix <- model.matrix(~0 + Condition, data=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ])
    cont.ref <- suppressWarnings(testNhoods(sim1.mylo, design=design.matrix,
                                                    model.contrasts=c("ConditionB-ConditionA"),
                                                    design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))
    form.ref <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition,
                                                    design.df=sim1.meta[colnames(nhoodCounts(sim1.mylo)), ]))

    # test for each column of output
    expect_equal(cont.ref$Pvalue, form.ref$Pvalue)
    expect_equal(cont.ref$SpatialFDR, form.ref$SpatialFDR)
    expect_equal(cont.ref$logFC, form.ref$logFC)
    expect_equal(cont.ref$`F`, form.ref$`F`)
    expect_equal(cont.ref$FDR, form.ref$FDR)
})

test_that("Providing a subset model.matrix is reproducible", {
    require(Matrix)
    set.seed(42)
    subset.samples <- sample(rownames(sim1.meta))
    exp.nh <- sum(Matrix::rowMeans(nhoodCounts(sim1.mylo)[, subset.samples]) >= 1)
    out.da <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                          min.mean=1,
                                          design.df=sim1.meta[subset.samples, ]))
    expect_identical(nrow(out.da), exp.nh)

    kd.ref1 <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                           min.mean=1,
                                           design.df=sim1.meta[subset.samples, ]))
    kd.ref2 <- suppressWarnings(testNhoods(sim1.mylo, design=~Condition, fdr.weighting="k-distance",
                                           min.mean=1,
                                           design.df=sim1.meta[subset.samples, ]))
    expect_identical(kd.ref1, kd.ref2)
})

sim1.meta$Condition_num <- paste0("Condition_num", c(1, 1, 1, 0, 0, 0))
sim1.meta$Replicate_num <- paste0("Replicate_num", c(1, 2, 3, 1, 2, 3))
sim1.meta$Replicate2 <- paste0("Replicate2", c(1, 2, 1, 2, 1, 2))

test_that("Singular Hessians are detectable and fail appropriately", {

    # having a singular Hessian depends on some of the staring values <- this test needs to
    # be reproducible and not depend on setting a specific seed. The easiest way might be to have
    # a variance component that is effectively 0.

    # collinear fixed and random effects
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=~Condition + (1|Condition),
                            design.df=sim1.meta)),
                 "Zero eigenvalues in D")

    # more variables than observations
    set.seed(42)
    expect_error(suppressWarnings(testNhoods(sim1.mylo, design=~Condition_num + (1|Replicate_num) + (1|Replicate2),
                                             design.df=sim1.meta)),
                 "Zero eigenvalues in D")
})

initializeFullZsim <- function(Z, cluster_levels, stand.cols=FALSE){
    # construct the full Z with all random effect levels
    n.cols <- ncol(Z)
    col.classes <- apply(Z, 2, class)
    i.z.list <- list()
    for(i in seq_len(n.cols)){
        i.class <- col.classes[i]
        if(i.class %in% c("factor")){ # treat as factors
            i.levels <- levels(Z[, i, drop=FALSE])
            i.levels <- as.factor(paste(sort(as.integer(i.levels))))
            i.z <- sapply(i.levels, FUN=function(X) (Z[, i] == X) + 0, simplify=TRUE)
        } else if(i.class %in% c("character")){
            i.levels <- unique(Z[, i, drop=FALSE])
            i.levels <- as.factor(paste(sort(as.integer(i.levels))))
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
            i.levels <- (unique(Z[, i]))
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

SimulateXZ <- function(N, n.fe, n.re, re.levels, fe.levels){

    # create a per-level mean effect for each FE
    if(length(fe.levels) != n.fe){
        stop("List entries need to match number of input fixed effects")
    }

    if(length(re.levels) != n.re){
        stop("List entries need to match number of input random effects")
    }

    # create the design matrices
    X <- matrix(0L, ncol=n.fe+1, nrow=N)
    X[, 1] <- 1
    colnames(X) <- c("Intercept", names(fe.levels))

    Z <- matrix(0L, ncol=n.re, nrow=N)

    for(i in seq_len(n.fe)){
        if(fe.levels[[i]] == 1){
            X[, i+1] <- sapply(seq_len(N), FUN=function(B){
                rnorm(1, mean=0, sd=1)
            })
        } else if(fe.levels[[i]] == 2){
            X[, i+1] <- sapply(seq_len(N), FUN=function(B){
                sample(c(0, 1), 1)
            })
            X[, i+1] <- as.factor(X[, i+1])
        }else{
            X[, i+1] <- sapply(seq_len(N), FUN=function(B){
                sample(seq_len(fe.levels[[i]]), 1)
            })
            X[, i+1] <- as.factor(X[, i+1])
        }
    }

    # Make categorical effects 0 or 1 (not 1 or 2)
    X[,2] <- X[,2] - 1

    for(j in seq_len(n.re)){
        if(re.levels[[j]] == 1){
            Z[, j] <- sapply(seq_len, FUN=function(R){
                rnorm(1, mean=1, sd=1)
            })
        } else{
            Z[, j] <- sapply(seq_len(N), FUN=function(R){
                sample(seq_len(re.levels[[j]]), 1)
            })
            Z[, j] <- factor(Z[, j], levels=c(1:re.levels[[j]]))
        }
    }
    colnames(Z) <- names(re.levels)

    sim.data <- do.call(cbind.data.frame, list(X, Z))
    return(sim.data)
}


SimulateY <- function(N, X, Z, fe.betas, re.sigmas,
                      dispersion, grand.mean, n.fe, n.re,
                      re.levels,
                      fe.levels){

    # create a per-level mean effect for each FE
    if(length(fe.levels) != n.fe){
        stop("List entries need to match number of input fixed effects")
    }

    if(length(re.levels) != n.re){
        stop("List entries need to match number of input random effects")
    }

    # construct the full Z
    random.levels <- sapply(seq_len(length(re.levels)), FUN=function(RX) {
        rx.name <- names(re.levels)[RX]
        paste(rx.name, seq_len(re.levels[[rx.name]]), sep="_")
    }, simplify=FALSE)
    names(random.levels) <- names(re.levels)

    full.Z <- initializeFullZsim(Z, random.levels)

    # get a combination over random effects
    # and sample each level from the same ~Normal(0, sigma)
    # note that the variance would be G if we also had random slopes
    re.thetas <- list()
    for(i in seq_len(length(re.levels))){
        i.re <- names(random.levels[i])
        i.levels <- length(random.levels[[i.re]])
        i.re.means <- rnorm(n=i.levels, 0, sd=sqrt(re.sigmas[[i.re]])) # sample a random effect value
        i.re.list <- sapply(seq_len(i.levels), FUN=function(X) i.re.means[X])
        names(i.re.list) <- random.levels[[i.re]]
        re.thetas[[i.re]] <- i.re.list
    }

    B <- full.Z %*% unlist(re.thetas)
    # map the fixed effects to mean values
    betas <- c(grand.mean, unlist(fe.betas))
    Beta <- X %*% betas

    i.error <- matrix(data = rnorm(N, mean=0, sd=0.001), ncol = 1)

    # construct the y.means equation, depending on desired distribution and FE/RE
    y.means <- exp(Beta + B)
    y.means <- y.means + i.error

    y.counts <- rnbinom(N, mu = y.means, size = dispersion)

    sim.data <- data.frame("Mean.Count"=y.counts)
    sim.data <- do.call(cbind.data.frame, list(sim.data, X, Z))

    return(sim.data)
}

N=150
fe.levels <- list("FE1"=2)
re.levels <- list("RE1"=10)
design.sim <- SimulateXZ(N=N, n.fe=length(fe.levels), n.re=length(re.levels), re.levels=re.levels, fe.levels=fe.levels)

n <- 10 # number of neighborhoods
sim.list <- c()

for (i in 1:n) {
    r.dispersion <- runif(1, min = 2, max = 2.5)
    fe.betas=list("FE1"=runif(1, min = 0.15, max = 0.25))
    re.sigmas=list("RE1"=runif(1, min = 0.05, max = 0.1))
    grand.mean=runif(1, min = 1, max = 2)
    sim.list[[i]]  <- SimulateY(N=N, X=sapply(design.sim[,1:2], as.numeric), Z=design.sim[,3, drop=FALSE], fe.betas=fe.betas, re.sigmas=re.sigmas, dispersion=r.dispersion, grand.mean=grand.mean, n.fe=length(fe.betas), n.re=length(re.sigmas), re.levels=re.levels, fe.levels=fe.levels)
}

names(sim.list) <- 1:n

y <- lapply(sim.list, `[[`, 1)
y_matrix <- matrix(unlist(y), nrow = 150, ncol = n)
colnames(y_matrix) <- paste("n", 1:n, sep = "_")

y_counts <- Matrix(t(y_matrix), sparse = T)
rownames(y_counts) <- 1:nrow(y_counts)
colnames(y_counts) <- 1:ncol(y_counts)
sim1.mylo@nhoodCounts <- y_counts

test_that("Providing a subset model.matrix is reproducible for glmm", {
    X <- sapply(design.sim[,2, drop = F], as.numeric)
    Z <- design.sim[,3, drop = F]
    design.df <- cbind(X, Z)
    colnames(design.df) <- c("ConditionB", "RE")
    rownames(design.df) <- 1:nrow(design.df)

    require(Matrix)
    subset.samples <- sample(rownames(design.df))
    exp.nh <- sum(Matrix::rowMeans(nhoodCounts(sim1.mylo)[, subset.samples]) >= 1)
    out.da <- suppressWarnings(testNhoods(sim1.mylo, design=~ConditionB + (1|RE), fdr.weighting="k-distance",
                                          min.mean=1,
                                          design.df=design.df[subset.samples, ]))
    expect_identical(nrow(out.da), exp.nh)

    set.seed(42)
    kd.ref1 <- suppressWarnings(testNhoods(sim1.mylo, design=~ConditionB + (1|RE), fdr.weighting="k-distance",
                                           min.mean=1,
                                           design.df=design.df[subset.samples, ]))
    kd.ref2 <- suppressWarnings(testNhoods(sim1.mylo, design=~ConditionB + (1|RE), fdr.weighting="k-distance",
                                           min.mean=1,
                                           design.df=design.df[subset.samples, ]))
    expect_equal(kd.ref1, kd.ref2, tolerance=1e-6)
})

