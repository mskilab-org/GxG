[![Build Status](https://travis-ci.org/mskilab/GxG.svg?branch=master)](https://travis-ci.org/mskilab/GxG)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/GxG.svg)](https://codecov.io/github/mskilab/GxG?branch=master)

# GxG

  GenomicRanges-based classes for representing, querying, and manipulating
  matrices of genomic data (ie a map of pairs of genomic coordinates to a
  numeric value).  Applications include visualization and analysis of Hi-C
  contact maps, barcode overlap in 10X, microhomology heatmaps, pairwise LD, and
  epistasis. 
  
## Install

1. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
install.packages('devtools')
install.packages('testthat')
```


2. Install dependent mskilab R packages

```{r}
devtools::install_github('mskilab/gUtils')
devtools::install_github('mskilab/gTrack')
```

3. Install GxG

```{r}
devtools::install_github('mskilab/GxG)
```

Documentation 
------------

[GxG Tutorial)](http://mskilab.com/GxG/tutorial.html)

[![alttext](https://github.com/mskilab/gGnome/raw/master/docs/gGnome.png)](http://mskilab.com/GxG/tutorial.html))

<!---
[GxG Developer Reference](docs/reference.md)
-->

<div id="attributions"/>

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill Cornell Medicine
> Core Member, New York Genome Center.

Funding sources
------------

<img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819f02b6a28750f79597c/1524111879079/DDCF.jpeg?format=1500w"
height="150" class ="center"> <img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819b8aa4a996c2d584594/1524111841815/BWF.png?format=500w"
height="150" class ="center">




```
