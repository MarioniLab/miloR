# miloR
_Milo_ is a method for differential abundance analysis on KNN graph from single-cell datasets.

<p align="center">
  <img src="docs/milo_schematic.png" width="500">
</p>

[![Build Status](https://travis-ci.com/MikeDMorgan/miloR.svg?branch=master)](https://travis-ci.com/MikeDMorgan/miloR)
[![Coverage](https://codecov.io/gh/MikeDMorgan/miloR/branch/master/graph/badge.svg)](https://codecov.io/gh/MikeDMorgan/miloR)


### Installation

```
## Install development version
devtools::install_github("MarioniLab/miloR") 
```

### Tutorials

1. [Basic Milo example on simulated dataset](https://raw.githack.com/MarioniLab/miloR/use_pkgdown/docs/articles/milo_demo.html)
2. [Milo example on mouse gastrulation dataset](https://raw.githack.com/MarioniLab/miloR/use_pkgdown/docs/articles/milo_gastrulation.html): this includes a demo for downstream analysis functions.
3. [Integrating Milo in scanpy/anndata workflow](https://github.com/MarioniLab/milo/blob/master/notebooks/milo_in_python.ipynb)

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






