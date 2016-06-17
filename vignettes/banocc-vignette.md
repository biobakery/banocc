Introduction to BAnOCC (Bayesian Analaysis Of Compositional Covariance)
================
Emma Schwager
2016-06-17

-   [How To Install](#how-to-install)
    -   [To install from within R](#to-install-from-within-r)
    -   [To install using a compressed file from bitbucket](#to-install-using-a-compressed-file-from-bitbucket)
    -   [To install directly from bitbucket](#to-install-directly-from-bitbucket)
-   [How To Run](#how-to-run)
    -   [The Model](#the-model)

Compositional data occur in many disciplines: geology, nutrition, economics, and ecology, to name a few. Data are compositional when each sample is sum-constrained. For example, mineral compositions describe a mineral in terms of the weight percentage coming from various elements; or taxonomic compositions break down a community by the fraction of individuals that come from a particular species. In ecology in particular, the covariance between features is often of interest to determine which species possible interact with each other. However, the sum constraint of compositional data makes naive measures inappropriate.

BAnOCC is a package for analyzing compositional covariance while accounting for the compositional structure. Briefly, the model assumes that the unobserved counts are log-normally distributed and then infers the correlation matrix of the log-basis.

How To Install
--------------

There are three options for installing BAnOCC:

-   Within R
-   Using compressed file from bitbucket
-   Directly from bitbucket

### To install from within R

**This is not yet available**

### To install using a compressed file from bitbucket

**This is not yet available**

### To install directly from bitbucket

``` bash
git clone https://<your-user-name>@bitbucket.org/biobakery/banocc.git
R CMD INSTALL banocc
```

How To Run
----------

The BAnOCC package contains three things:

-   `bayesStanModel`, which is the BAnOCC model in the `rstan` format
-   `run_banocc`, a wrapper function for `rstan::sampling` that samples from the model and returns a list with various useful elements
-   `banocc_data`, a small test dataset

``` r
# This code is not run
library(banocc)
data(banocc_data)
bayesModel   <- rstan::stan_model(model_code=bayesStanModel)
bayes_output <- banocc::run_banocc(banocc_data, bayes_model=bayesModel)
```

For a full description of the possible parameters for `run_banocc`, their default values, and the output, see

``` r
?run_banocc 
```
