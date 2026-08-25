# 2.9.2 (2026-08-23)
+ **New feature:** Introduce `testGeneticNhoods()` for genome-wide cell state QTL testing: the null model is fitted once per neighbourhood and every variant is scored against it, rather than refitting the full model per neighbourhood and variant. Use `testNhoods()` for designs with a single variable of interest such as age or disease status.
+ **New feature:** Introduce `refineGeneticHits()` to refit variants passing a screening threshold with the exact model, and `plotGeneticQQ()` for two-panel QQ plots with genomic inflation factors.
+ **New feature:** `testNhoods()` now rejects a mixed model design whose variance components cannot be identified, before fitting anything. A component is identifiable only if its design contributes something the fixed effects do not: `rank([X Z_j])` is always at least `rank(X)`, and equality means the levels of `Z_j` are fully determined by the fixed effect design. Rank deficiency in `X` itself is rejected on the same pass. Previously a design like `~Condition + (1|Condition)` surfaced downstream as a singular Hessian, which depended on the random starting values. With the dispersion estimated as a variance component the pseudo-variance is better conditioned and such a fit converges instead, leaving the degeneracy visible only as Satterthwaite degrees of freedom of order 1e-28 and p-values of exactly 1, which is easy to miss.
+ **New feature:** `testNhoods()` rejects a design that asks for more random effect levels than there are observations. A relatedness matrix gives the genetic effect one level per observation, which is identified on its own because its covariance is `sigma_g * K` with `K` fixed - one variance parameter, not n - but a second random effect on top of it is not, since there is no residual left to distinguish it from.
+ **Behaviour change:** the negative binomial overdispersion in `fitGLMM()`, `testNhoods()`, `fitGeneticNullGlmm()` and `testGeneticNhoods()` is now estimated as an additional variance component on the REML objective (`disp_as_vc=TRUE`, the default). Previously the dispersion was estimated by a golden-section search on the conditional NB likelihood evaluated at the current fitted means. That is circular, because those means already contain the random effect BLUPs: the random effects absorb the overdispersion before it is estimated, so the size parameter is driven to the Poisson limit and the random effect variances inflate to compensate. Under `disp.as.vc` the dispersion enters as sigma_0 with `dV*/dsigma_0 = I`, so it is estimated by the same Fisher scoring, the same step-halving at the boundary and the same convergence test as the other components, rather than off-line. The legacy behaviour is incorrect and is retained only for reproducing earlier results - set `disp_as_vc=FALSE` to use it.
+ **Behaviour change:** The genetic variance component is fitted on `K = I + E`. With both `sigma_0 I` and `sigma_g K` in the pseudo-variance the two partial derivatives are `I` and `K`, which for a relatedness matrix of nominally unrelated individuals are numerically the same matrix - the REML information then has two identical columns and is exactly rank 1 when `K = I`. Writing `E = K - I` and collecting terms gives `V* = W + tau I + gamma E` with `tau = sigma_0 + sigma_g` and `gamma = sigma_g`. This is a linear change of basis with Jacobian `[[1,1],[0,1]]`, so the information transforms as `J' I J` and its rank is preserved. Verified behaviour-preserving: fits agree with the previous parameterisation to seven or eight significant figures on every seed tested, with identical iteration counts. `testNhoods()` checks that a supplied relatedness matrix can actually separate the genetic variance from the negative binomial overdispersion, and refuses the fit when it cannot. Writing `K = I + E`, the two components enter the pseudo-variance as `sigma_g (I + E)` and `sigma_0 I`, so they differ only through the off-diagonal structure; when `E` is negligible the partial derivatives with respect to the two parameters are the same matrix and the REML information has two identical columns. The statistic is `||E||_F / ||K||_F`, the share of the matrix norm carried off the diagonal, which scales with cohort size for fixed per-pair relatedness - `||E||_F^2` grows as `n(n-1)` against a diagonal contributing only `n` - so one threshold serves any `n`. Fits are refused at or below 0.1 and warned between 0.1 and 0.35; `force=TRUE` downgrades the refusal to a warning.
+ **Bug fix:** `invertPseudoVar()` computed `Z * G` and `I + Z' W^-1 * (Z * G)` in two `omp parallel sections`. The second reads the product the first is still writing, so the block was a data race rather than a speed-up; the two steps should be sequential.
+ + **Bug fix:** Fix error in `nbLogLik()`, which omitted `lgamma(y + r)`, used `r * (1 - mu/(mu + r))` in place of `r * log(r / (r + mu))`, and carried `lgamma(y + 1)` with the wrong sign. The first two terms depend on `r`. Together with the change above, null p-value calibration improves from lambda_GC = 0.26 to 1.02.
+ **Bug fix:** the GLMM path of `testNhoods()` failed with `could not find function "bpstopOnError"` for anyone who had not separately attached `BiocParallel`, which is the normal case for `library(miloR)`. `glmmWrapper()` passed the stop-on-error setting through `BPOPTIONS`, and `bplapply` resolves that option from the calling frame rather than from its own namespace. 
+ **Bug fix:** The variance component update is relaxed against period-2 limit cycles. The pseudolikelihood iteration alternates variance components -> fixed effects and BLUPs -> fitted means -> variance components, and on strongly overdispersed data that composite map can have a stable 2-cycle rather than a fixed point. A relaxed update `sigma <- sigma + alpha (update - sigma)` collapses the cycle for `alpha < 1` and cannot move the estimate, since at a fixed point the update equals the current value for every `alpha`; the feasible set is convex, so the `tau > gamma` constraint survives a convex combination too. `alpha` starts at 1 and halves only when a cycle is actually detected. The cycle test is deliberately strict - the iterate must return to within a millionth of its position two steps ago while still moving appreciably. A loose threshold will activate on ordinary damped oscillation while a genuine limit cycle sits at a ratio of order 1e-15, leaving nine orders of margin.
+ **Bug fix:** the REML score for the variance components in `fitPLGlmm()` and `fitGeneticPLGlmm()` mixed two bases. The score is now taken entirely in the `P` basis, and the function no longer accepts a second basis to get wrong.
+ **Bug fix:** the model offset was left in the working response used for the mixed model equations and for the variance component estimation. The offset is part of the linear predictor but is not a column of `X`, so it is neither annihilated by `P` nor held fixed by the mixed model equations - the intercept absorbed it and then fed it back into the linear predictor on the next iteration, growing by `mean(offset)` per iteration without ever converging. The offset is now removed before the working response is used.
+ **Bug fix:** the golden-section bracket for the legacy dispersion was `[phi/2, phi]`, whose upper end is the current estimate, so the search could only ever return a smaller value. The dispersion halved on every iteration until the convergence check stopped it, regardless of the data - from a starting value of 4 it reached 0.0155 = 4 x 2^-8 on data with a true size of 5. The bracket is now two-sided, `[phi/2, 2 phi]`, capped at 1e4 as in `fitGeneticNullGlmm()`. This is only relevant for legacy models.
+ **Bug fix:** the dispersion update was skipped once two consecutive searches agreed to within 1e-2. The search now runs while the model is still iterating, and its bracket tolerance is resolved relative to the bracket rather than at a flat 1e-2, which was coarse enough to prevent the outer fixed point from converging.
+ **Bug fix:** `nnlsSolve()` initialised the Lawson-Hanson passive and active sets on the wrong side of zero. A zero starting vector - which is what both the `HE-NNLS` solver and the negative-variance fallback supply - put every index in the passive set and left the active set empty, and an empty active set is itself a termination condition, so the routine returned the zero vector it was handed without doing any work. Every `HE-NNLS` fit reported variance components of exactly zero. The pivot search is now restricted to the active set, and the loops are bounded.
+ **Bug fix:** the residual column of the vectorised Haseman-Elston design was `vec(I)`. Under REML the moment being regressed is `E[P y* y*' P] = P W P + sum_j sigma_j P Z_j Z_j' P`, so that column must be `vec(P W P)`; `W = diag(1/phi + 1/mu_i)` is not a multiple of the identity, so the residual column could not absorb the working weights and the difference was taken out of the random effect instead, driving it negative. The ML variants use `vec(W)` for the same reason.
+ **Bug fix:** the genetic random effect applied the relatedness matrix twice. `fitGLMM()` sets the genetic block of the design matrix to `t(chol(Kin))`, i.e. the Cholesky factor L, so that `Z_g Z_g' = L L' = K`. The corresponding block of G should then be `sigma_g I` but the C++ `initialiseG_G()`, which rebuilds G on every iteration, used `sigma_g Kin` instead. The pseudo-variance therefore contained `sigma_g L K L'` while every partial derivative of V* with respect to `sigma_g` remained `K`. G and its inverse are now built in the whitened parameterisation, matching the design matrix and the derivatives; internal consistency of estimates is now resolved.
+ **Bug fix:**F ix a shadowed declaration in `fitGeneticPLGlmm()` that redeclared `VP_partial` inside the ML branch, so the `Vpartial` returned to R was empty for ML with Fisher scoring and `varCovar()` ran on an empty list.
+ **Bug fix:** Negative variance component updates under Fisher scoring now retreat along the ascent direction by step-halving, as in `fitGeneticNullGlmm()`, rather than switching solver and clamping to the boundary, which pinned the component there for every later iteration.
+ **Deprecation:** Retire `initialiseG_G()`, `invGmat_G()`, `subMatG()` and `broadcastInverseMatrix()`. The whitened parameterisation makes G and its inverse diagonal, so the genetic path no longer inverts the relatedness matrix at all - one dense O(n^3) inverse per fit.
  
# 2.3.1 (2024-11-20)
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

