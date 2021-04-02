# miloR
_Milo_ is a method for differential abundance analysis on KNN graph from single-cell datasets. For more details, read [our preprint](https://www.biorxiv.org/content/10.1101/2020.11.23.393769v1). 

<p align="center">
  <img src="docs/milo_schematic.png" width="500">
</p>

[![Build Status](https://travis-ci.com/MarioniLab/miloR.svg?branch=master)](https://travis-ci.com/MarioniLab/miloR)
[![Coverage](https://codecov.io/gh/MarioniLab/miloR/branch/master/graph/badge.svg)](https://codecov.io/gh/MarioniLab/miloR)
[![R-CMD-check](https://github.com/MarioniLab/miloR/actions/workflows/RCMD_check.yml/badge.svg)](https://github.com/MarioniLab/miloR/actions/workflows/RCMD_check.yml)

### Installation

```
## Install development version
devtools::install_github("MarioniLab/miloR") 
```

### Tutorials

1. [Basic Milo example on simulated dataset](https://rawcdn.githack.com/MarioniLab/miloR/3646391023f600bae00efd9d940b888503d7a536/docs/articles/milo_demo.html)
2. [Milo example on mouse gastrulation dataset](https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#5_Finding_markers_of_DA_populations): this includes a demo for downstream analysis functions.
3. [Integrating Milo in scanpy/anndata workflow](https://github.com/MarioniLab/milo_analysis_2020/blob/main/notebooks/milo_in_python.ipynb)

### Example work flow
An example of the `Milo` work flow to get started:

```{r}
data(sim_trajectory)
milo.meta <- sim_trajectory$meta
milo.obj <- Milo(sim_trajectory$SCE)
milo.obj
```

Build a graph and neighbourhoods.

```{r}
milo.obj <- buildGraph(milo.obj, k=20, d=30)
milo.obj <- makeNhoods(milo.obj, k=20, d=30, refined=TRUE, prop=0.2)
```

Calculate distances, count cells according to an experimental design and perform DA testing.

```{r}
milo.obj <- calcNhoodDistances(milo.obj, d=30)
milo.obj <- countCells(milo.obj, samples="Sample", meta.data=milo.meta)

milo.design <- as.data.frame(xtabs(~ Condition + Sample, data=milo.meta))
milo.design <- milo.design[milo.design$Freq > 0, ]

milo.res <- testNhoods(milo.obj, design=~Condition, design.df=milo.design)
head(milo.res)
```

### Support

For any question, feature request or bug report please create a new issue in this repository.






