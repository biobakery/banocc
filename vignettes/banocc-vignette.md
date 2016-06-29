Introduction to BAnOCC (Bayesian Analaysis Of Compositional Covariance)
================
Emma Schwager
2016-06-29

-   [Introduction](#introduction)
-   [How To Install](#how-to-install)
    -   [To install from within R](#to-install-from-within-r)
    -   [To install using a compressed file from bitbucket](#to-install-using-a-compressed-file-from-bitbucket)
    -   [To install directly from bitbucket](#to-install-directly-from-bitbucket)
-   [How To Run](#how-to-run)
    -   [Data and Prior Input](#data-and-prior-input)
    -   [Sampling Control](#sampling-control)
    -   [Output Control](#output-control)
-   [Assessing Convergence](#assessing-convergence)
-   [Choosing Priors](#choosing-priors)
-   [The Model](#the-model)
-   [References](#references)

Introduction
------------

Compositional data occur in many disciplines: geology, nutrition, economics, and ecology, to name a few. Data are compositional when each sample is sum-constrained. For example, mineral compositions describe a mineral in terms of the weight percentage coming from various elements; or taxonomic compositions break down a community by the fraction of individuals that come from a particular species. In ecology in particular, the covariance between features is often of interest to determine which species possibly interact with each other. However, the sum constraint of compositional data makes naive measures inappropriate.

BAnOCC is a package for analyzing compositional covariance while accounting for the compositional structure. Briefly, the model assumes that the unobserved counts are log-normally distributed and then infers the correlation matrix of the log-basis (see [The Model](#the-model) section for a more detailed explanation). The inference is made using No U-Turn Sampling for Hamiltonian Monte Carlo (Hoffman and Gelman 2014) as implemented in the `rstan` R package (Stan Development Team 2015).

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
banocc_model  <- rstan::stan_model(model_code=banocc_model) 
banocc_output <- banocc::run_banocc(C = banocc_data, banocc_model=banocc_model) 
```

For a full and complete description of the possible parameters for `run_banocc`, their default values, and the output, see

``` r
?run_banocc 
```

### Data and Prior Input

The only required input to `run_banocc` is the dataset, `C` and the compiled stan model, `banocc_model`. The dataset is assumed to be *n* × *p*, with *n* samples and *p* features. The row sums are therefore required to be less than 1 for all samples.

The hyperparameter values can be specified as input. Their names correspond to the parameters in the plate diagram figure (see section [The Model](#the-model)). For example,

``` r
# This code is not run
p <- ncol(banocc_data)
b_hp <- banocc::run_banocc(C = banocc_data, 
                           banocc_model = banocc_model,
                           n = rep(0, p),
                           L = 10 * diag(p), 
                           a = rep(1,   p),
                           b = rep(0.5, p),
                           eta = 1)
```

The hyperparameter values for ***s*** can be alternatively specified as the means and variances of the prior using the parameters `sd_mean` and `sd_var`. For example,

``` r
# This code is not run
b_sd_mean <- banocc::run_banocc(C = banocc_data,
                                banocc_model = banocc_model,
                                sd_mean = rep(2, p),
                                sd_var  = rep(4, p))
```

### Sampling Control

There are several options to control the behavior of the HMC sampler. This is simply a call to `rstan::sampling`, and so many of the parameters are the same. The number of chains, iterations, and warmup iterations as well as the rate of thinning can be specified using the same parameters. For example, the following code gives a total of three iterations from each of two chains. These parameters are used only for brevity and are NOT recommended in practice.

``` r
# This code is not run
b_sampling <- banocc::run_banocc(C = banocc_data,
                                 banocc_model = banocc_model,
                                 chains = 2,
                                 iter = 11,
                                 warmup = 5,
                                 thin = 2)
```

The number of cores used for sampling on a multi-processor machine can also be specified, which allows chains to run in parallel and therefore decreases computation time. Since its purpose is running chains in parallel, computation time will decrease as cores are added up to when the number of cores and the number of chains are equal.

``` r
# This code is not run
b_cores <- banocc::run_banocc(C = banocc_data,
                              banocc_model = banocc_model,
                              chains = 2,
                              cores = 2)
```

The initial values are sampled from the priors by default because the positive definite constraint on the covariance matrix makes it difficult for `rstan::sampling` to generate initial values automatically. They can also be set to a particular value by using a list. `WChol` is the cholesky decomposition of ***W*** that is used by the sampler for computational reasons.

``` r
# This code is not run
init <- list(list(m = rep(0, p),
                  WChol = t(chol(diag(p))),
          s = rep(1, p)),
             list(m = runif(p),
                  WChol = t(chol(2 * diag(p))),
                  s = runif(p, 0.1, 2)))
b_init <- banocc::run_banocc(C = banocc_data,
                             banocc_model = banocc_model,
                             chains = 2,
                             init = init)
```

More specific control of the sampler's behavior comes from the `control` argument to `rstan::sampling`. Details about this argument can be found in the help for the `rstan::stan` function:

``` r
?stan
```

### Output Control

There are several parameters that control the type of output which is returned.

The width of the returned credible intervals is controlled by `conf_alpha`. A 100%\*(1 − *α*<sub>conf</sub>) credible interval is returned:

``` r
# This code is not run

# Get 90% credible intervals
b_90 <- banocc::run_banocc(C = banocc_data,
                           banocc_model = banocc_model,
                           conf_alpha = 0.1)
# Get 99% credible intervals
b_99 <- banocc::run_banocc(C = banocc_data,
                           banocc_model = banocc_model,
                           conf_alpha = 0.01)
```

Two types of output can be requested for each correlation that are not included by default:

1.  The smallest credible interval width that includes zero
2.  The scaled neighborhood criterion, or SNC (Li and Lin 2010)

``` r
# This code is not run

# Get the smallest credible interval width that includes zero
b_min_width <- banocc::run_banocc(C = banocc_data,
                                  banocc_model = banocc_model,
                                  get_min_width = TRUE)

# Get the scaled neighborhood criterion
b_snc <- banocc::run_banocc(C = banocc_data,
                            banocc_model = banocc_model,
                            calc_snc = TRUE)
```

Detailed statements about the function's execution can also be printed using the `verbose` argument. The relative indentation of the verbose output indicates the nesting level of the function. The starting indentation can be set with `num_level`.

Assessing Convergence
---------------------

**Coming Soon!**

Choosing Priors
---------------

**Coming Soon!**

The Model
---------

A pictoral representation of the model is shown below. Briefly, the basis (or unobserved, unrestricted counts) for each sample is assumed to be a lognormal distribution with parameters ***m*** and ***S***. The prior on ***m*** is a normal distribution parametrized by mean ***n*** and variance-covariance matrix ***L***. Since we are interested in the correlation structure, we break ***S*** into a correlation matrix ***W*** and a vector of standard deviations ***s***. The prior on ***W*** is an LKJ distribution (Lewandowski, Kurowicka, and Joe 2009) with shrinkage parameter *η*, while the prior on eash *s*<sub>*j*</sub> is a gamma prior with shape *a*<sub>*j*</sub> and rate *b*<sub>*j*</sub>.

![plate-diagram](Figure3.png)

If we print the model, we can actually see the code. It is written in the format required by the `rstan` package, since `banocc` uses this package to sample from the model.

``` r
# This code is not run
cat(banocc_model)
```

References
----------

Hoffman, Matthew D., and Andrew Gelman. 2014. “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.” *J. Mach. Learn. Res.* 15 (1). JMLR.org: 1593–1623. <http://dl.acm.org/citation.cfm?id=2627435.2638586>.

Lewandowski, Daniel, Dorota Kurowicka, and Harry Joe. 2009. “Generating Random Correlation Matrices Based on Vines and Extended Onion Method.” *Journal of Multivariate Analysis* 100 (9): 1989–2001. <http://dx.doi.org/10.1016/j.jmva.2009.04.008>.

Li, Qing, and Nan Lin. 2010. “The Bayesian Elastic Net.” *Bayesian Anal.* 5 (1). International Society for Bayesian Analysis: 151–70. <http://dx.doi.org/10.1214/10-BA506>.

Stan Development Team. 2015. “Stan: A C++ Library for Probability and Sampling, Version 2.10.0.” <http://mc-stan.org/>.
