# 2.0.X (2024-10-30)
+ Warning on GLMM if glmm.solver not set
+ Bug fix in model contrasts vignette with multiple contrasts

# 2.0.1 (2024-04-30)
+ Introduce NB-GLMM into Milo 2.0 for random effect variables and modelling dependencies between observations
+ Diagnostic function for checking model separation for experimental variables, i.e. splitting zero from non-zero counts perfectly
+ Vignette describing basic usage of GLMM functions in `testNhoods`

# 1.7.1 (2023-02-15)
+ Patch to fix `NA` plotting in `plotDABeeswarm`

# 1.60 (2022-11-02)
+ Vignette describing the use of contrasts in `testNhoods`

# 1.5.0 (2022-04-27)
+ Introduce plotting function to visualise neighbourhood count distributions for nhoods interest: `plotNhoodCounts`. Implemented by Nick Hirschmüller

# 1.3.1 (2022-01-07)
+ Fix bug in `findNhoodGroupMarkers` to merge on gene IDs explicitly
+ Fix bug in `makeNhoods` to include index cell in nhoods() matrix
+ Introduce graph-based neighbourhood definition - allows full compatibility with graph-only batch correction and graphs constructed by third-party tools
+ Introduce graph-based spatial FDR correction to obviate the need for any distance calculations
+ Add vignette to describe the use and application of contrasts in `testNhoods`
+ Patch to correct SpatialFDR with sparse nhoods where density is ~0

# 1.1.0 (2021-10-12)
+ Fix bug in testNhoods to use user-specific reduced dimensions
+ Vignettes now include set rownames() to avoid confusion
+ Numerous doc-string typo fixes

# 0.99.1 (2021-03-13)
+ Fix model normalisation bug - now using TMM normalisation by default. Log(M_s) offset can be used by passing `norm.method="logMS"` to `testNhoods`.

# 0.99.0 (2021-03-04)
+ Submitted to Bioconductor

