# 1.3.1 (2022-01-07)
+ Fix bug in `findNhoodGroupMarkers` to merge on gene IDs explicitly
+ Introduce graph-based neighbourhood definition - allows full compatibility with graph-only batch correction and graphs constructed by third-party tools
+ Introduce graph-based spatial FDR correction to obviate the need for any distance calculations

# 1.1.0 (2021-10-12)
+ Fix bug in testNhoods to use user-specific reduced dimensions
+ Vignettes now include set rownames() to avoid confusion
+ Numerous doc-string typo fixes

# 0.99.1 (2021-03-13)
+ Fix model normalisation bug - now using TMM normalisation by default. Log(M_s) offset can be used by passing `norm.method="logMS"` to `testNhoods`.

# 0.99.0 (2021-03-04)
+ Submitted to Bioconductor

