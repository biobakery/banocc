---
title: "Introduction to BAnOCC (Bayesian Analaysis of Compositional Covariance)"
author: "Emma Schwager"
date: "2016-06-14"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{BAnOCC (Bayesian Analysis of Compositional Covariance)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

BAnOCC is a package for analyzing compositional covariance. Briefly, the model assumes that the unobserved counts are log-normally distributed and then infers the correlation matrix of the log-basis. 

## How To Install
There are three options for installing BAnOCC:

* Within R
* Using compressed file from bitbucket
* Directly from bitbucket

### To install from within R 
**This is not yet available**

```r
source("<bioconductor URL>")
biocLite("BAnOCC")
```

### To install using a compressed file from bitbucket
**This is not yet available**

### To install directly from bitbucket
```bash
git clone https://<your-user-name>@bitbucket.org/biobakery/banocc.git
R CMD INSTALL banocc
```

## To Run

```r
# This code is not run
library(banocc)
data(banocc_data)
bayesModel   <- rstan::stan_model(model_code=bayesStanModel)
bayes_output <- banocc::run_banocc(banocc_data)
```
