context("Step 3 - testing get_posterior_quantiles")

load("testthat_objects/sample_stan_data.RData")
load("testthat_objects/sample_stan_fit.RData")
posterior_samples <- rstan::extract(Fit)

test_that("get_posterior_quantiles returns a list", {
    test_is_list <- function(probs, list, p.names){
        get_posterior_quantiles(posterior_samples=posterior_samples,
                                probs=probs, list=list,
                                parameter.names=p.names)
    }
    expect_is(test_is_list(0.5, FALSE, "mu"),             "list")
    expect_is(test_is_list(0.5, TRUE, "mu"),              "list")
    expect_is(test_is_list(0.5, FALSE, c("mu", "Sigma")), "list")
    expect_is(test_is_list(0.5, TRUE, c("mu", "Sigma")),  "list")
    expect_is(test_is_list(0, FALSE, "mu"),               "list")
    expect_is(test_is_list(0, TRUE, "mu"),                "list")
    expect_is(test_is_list(1, FALSE, "mu"),               "list")
    expect_is(test_is_list(1, TRUE, "mu"),                "list")
    expect_is(test_is_list(c(0, 0.5, 1), FALSE, "mu"),    "list")
    expect_is(test_is_list(c(0, 0.5, 1), TRUE, "mu"),     "list")
})

test_that("get_posterior_quantiles returns with list parameter names", {
    test_names <- function(probs, list, idx){
        sort(names(get_posterior_quantiles(
            posterior_samples=posterior_samples,
            probs=probs, list=list,
            parameter.names=sort(names(posterior_samples))[idx])))
    }
    expected_names <- sort(names(posterior_samples))
    idx <- seq(1, length(posterior_samples))
    expect_equal(test_names(c(0, 0.5, 1), TRUE, 1),    expected_names[1])
    expect_equal(test_names(c(0, 0.5, 1), FALSE, 1),   expected_names[1])
    expect_equal(test_names(c(0, 0.5, 1), TRUE, idx),  expected_names)
    expect_equal(test_names(c(0, 0.5, 1), FALSE, idx), expected_names)
})

test_that("get_posterior_quantiles elt dimensions are correct when list=FALSE", {
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=c(0, 0.5, 1), list=FALSE,
                                  parameter.names=c("mu", "Sigma"))
    expect_equal(dim(pq$mu),    c(3, Data$P))
    expect_equal(dim(pq$Sigma), c(3, Data$P, Data$P))
})

test_that("get_posterior_quantiles elt lengths are correct when list=TRUE", {
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=c(0, 0.5, 1), list=TRUE,
                                  parameter.names=c("mu", "Sigma"))
    expect_equal(length(pq$mu), 3)
    expect_equal(length(pq$Sigma), 3)
})

test_that("get_posterior_quantiles quantile dimensions are correct when list=TRUE", {
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=c(0, 0.5, 1), list=TRUE,
                                  parameter.names=c("mu", "Sigma"))
    expect_equal(dim(pq$mu[[1]]),    NULL)
    expect_equal(length(pq$mu[[1]]), Data$P)
    expect_equal(dim(pq$Sigma[[1]]), c(Data$P, Data$P))
})

test_that("get_posterior_quantiles matches eltwise calcn for vectors when list=FALSE",{
    probs <- c(0, 0.5, 1)
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=probs, list=FALSE,
                                  parameter.names="mu")
    for (i in seq_len(Data$P)){
        expect_equal(pq$mu[, i],
                     quantile(posterior_samples$mu[, i], probs=probs))
    }
})

test_that("get_posterior_quantiles matches eltwise calcn for vectors when list=TRUE", {
    probs <- c(0, 0.5, 1)
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=probs, list=TRUE,
                                  parameter.names="mu")
    for (i in seq_len(Data$P)){
        q <- unname(quantile(posterior_samples$mu[, i], probs=probs))
        for (k in seq_along(probs)){
            expect_equal(pq$mu[[k]][i], q[k])
        }
    }
})

test_that("get_posterior_quantiles matches eltwise calcn for matrices when list=FALSE", {
    probs <- c(0, 0.5, 1)
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=probs, list=FALSE,
                                  parameter.names="Sigma")
    for (i in seq_len(Data$P)){
        for (k in seq_len(Data$P)){
            q <- quantile(posterior_samples$Sigma[, i, k], probs=probs)
            expect_equal(pq$Sigma[, i, k], q)
        }
    }
})

test_that("get_posterior_quantiles matches eltwise calcn for matrices when list=TRUE", {
    probs <- c(0, 0.5, 1)
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=probs, list=TRUE,
                                  parameter.names="Sigma")
    for (i in seq_len(Data$P)){
        for (k in seq_len(Data$P)){
            q <- unname(quantile(posterior_samples$Sigma[, i, k],
                                 probs=probs))
            for (l in seq_along(probs)){
                expect_equal(pq$Sigma[[l]][i, k], q[l])
            }
        }
    }
})
