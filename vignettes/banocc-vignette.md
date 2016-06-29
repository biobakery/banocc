Introduction to BAnOCC (Bayesian Analaysis Of Compositional Covariance)
================
Emma Schwager
2016-06-29

-   [How To Install](#how-to-install)
    -   [To install from within R](#to-install-from-within-r)
    -   [To install using a compressed file from bitbucket](#to-install-using-a-compressed-file-from-bitbucket)
    -   [To install directly from bitbucket](#to-install-directly-from-bitbucket)
-   [How To Run](#how-to-run)
    -   [The Model](#the-model)
-   [References](#references)

Compositional data occur in many disciplines: geology, nutrition, economics, and ecology, to name a few. Data are compositional when each sample is sum-constrained. For example, mineral compositions describe a mineral in terms of the weight percentage coming from various elements; or taxonomic compositions break down a community by the fraction of individuals that come from a particular species. In ecology in particular, the covariance between features is often of interest to determine which species possibly interact with each other. However, the sum constraint of compositional data makes naive measures inappropriate.

BAnOCC is a package for analyzing compositional covariance while accounting for the compositional structure. Briefly, the model assumes that the unobserved counts are log-normally distributed and then infers the correlation matrix of the log-basis. The inference is made using No U-Turn Sampling for Hamiltonian Monte Carlo (Hoffman and Gelman 2014) as implemented in the `rstan` R package (Stan Development Team 2015).

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

We first need to load the package:

``` r
library(banocc)
```

    ## Loading required package: rstan

    ## Loading required package: ggplot2

    ## Loading required package: StanHeaders

    ## rstan (Version 2.10.1, packaged: 2016-06-24 13:22:16 UTC, GitRev: 85f7a56811da)

    ## For execution on a local, multicore CPU with excess RAM we recommend calling
    ## rstan_options(auto_write = TRUE)
    ## options(mc.cores = parallel::detectCores())

The BAnOCC package contains three things:

-   `bayesStanModel`, which is the BAnOCC model in the `rstan` format
-   `run_banocc`, a wrapper function for `rstan::sampling` that samples from the model and returns a list with various useful elements
-   `banocc_data`, a small test dataset

The simplest way to run the model is to load the test dataset, compile the model, and sample from it:

``` r
# This code is not run 
data(banocc_data) 
banocc_model  <- rstan::stan_model(model_code=bayesStanModel) 
banocc_output <- banocc::run_banocc(banocc_data, bayes_model=banocc_model) 
```

For a full description of the possible parameters for `run_banocc`, their default values, and the output, see

``` r
?run_banocc 
```

### The Model

A pictoral representation of the model is shown below. Briefly, the basis (or unobserved, unrestricted counts) for each sample is assumed to be a lognormal distribution with parameters ***m*** and ***S***. The prior on ***m*** is a normal distribution parametrized by mean ***n*** and variance-covariance matrix ***L***. Since we are interested in the correlation structure, we break ***S*** into a correlation matrix ***W*** and a vector of standard deviations ***s***. The prior on ***W*** is an LKJ distribution (Lewandowski, Kurowicka, and Joe 2009) with shrinkage parameter *η*, while the prior on eash *s*<sub>*j*</sub> is a gamma prior with shape *a*<sub>*j*</sub> and rate *b*<sub>*j*</sub>.

![plate-diagram](Figure3.png)

If we print the model, we can actually see the code. It is written in the format required by the `rstan` package, since `banocc` uses this package to sample from the model.

``` r
# This code is not run
cat(bayesStanModel)
```

References
----------

Hoffman, Matthew D., and Andrew Gelman. 2014. “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.” *J. Mach. Learn. Res.* 15 (1). JMLR.org: 1593–1623. <http://dl.acm.org/citation.cfm?id=2627435.2638586>.

Lewandowski, Daniel, Dorota Kurowicka, and Harry Joe. 2009. “Generating Random Correlation Matrices Based on Vines and Extended Onion Method.” *Journal of Multivariate Analysis* 100 (9): 1989–2001. doi:[http://dx.doi.org/10.1016/j.jmva.2009.04.008](https://doi.org/http://dx.doi.org/10.1016/j.jmva.2009.04.008).

Stan Development Team. 2015. “Stan: A C++ Library for Probability and Sampling, Version 2.10.0.” <http://mc-stan.org/>.
