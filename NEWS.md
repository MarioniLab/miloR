# 2.9.2 (2026-08-19)
+ **Behaviour change:** the negative binomial overdispersion in `fitGeneticNullGlmm()` and `testGeneticNhoods()` is now estimated as an additional variance component on the REML objective (`disp_as_vc=TRUE`, the default). Previously the dispersion was estimated by a golden-section search on the conditional NB likelihood evaluated at the current fitted means. That is circular, because those means already contain the random effect BLUPs: the random effects absorb the overdispersion before it is estimated, so the size parameter is driven to the Poisson limit (order 1e4 against a true value of 5 in simulation) and the random effect variances inflate to compensate (genetic component 0.372 against a true 0.1225). The legacy behaviour is incorrect and is retained only for reproducing earlier results - set `disp_as_vc=FALSE` to use it.
+ Fix `nbLogLik()`, which did not evaluate the negative binomial log-likelihood: it omitted `lgamma(y + r)`, used `r * (1 - mu/(mu + r))` in place of `r * log(r / (r + mu))`, and carried `lgamma(y + 1)` with the wrong sign. The first two terms depend on `r`, so `phiGoldenSearch()` was maximising the wrong function of the dispersion - it peaked at r = 0.785 against a true value of 5. This affects every GLMM fit, including `fitPLGlmm()` and `fitGeneticPLGlmm()`, since the working weights are `W = diag(1/phi + 1/mu)`. Together with the change above, null p-value calibration improves from lambda_GC = 0.26 to 1.02.
+ Introduce `testGeneticNhoods()` for genome-wide cell state QTL testing: the null model is fitted once per neighbourhood and every variant is scored against it, rather than refitting the full model per neighbourhood and variant. Use `testNhoods()` for designs with a single variable of interest such as age or disease status.
+ Introduce `refineGeneticHits()` to refit variants passing a screening threshold with the exact model, and `plotGeneticQQ()` for two-panel QQ plots with genomic inflation factors.

# 2.3.X (2024-11-20)
+ Allow groupNhoods() to retain original behaviour or force intuitive grouping such that no nhood group has discordant LFCs

# 2.2.0 (2024-10-30)
+ Warning on GLMM if glmm.solver not set
+ Bug fix in model contrasts vignette with multiple contrasts
+ testNhoods will error if N<60 and using GLMM - introduce force=TRUE to override (with a warning)
+ DA nhoods can be emphasised in plotNhoodGraphDA with `highlight.da`

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

