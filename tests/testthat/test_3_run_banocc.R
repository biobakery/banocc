context("Step 3 - running BAnOCC")

test_that("get_alpha_beta gives error if alpha and beta not provided", {
    err_string <- "provide both 'alpha' and 'beta'"
    expect_error(get_alpha_beta(alpha=c(1, 2), beta=NULL, p=2),  err_string)
    expect_error(get_alpha_beta(alpha=NULL, beta=c(1, 2), p=2),  err_string)
    expect_error(get_alpha_beta(alpha=NULL, beta=NULL, p=2),     err_string)
})

test_that("get_alpha_beta gives error if sd_mean and sd_var not both provided", {
    err_string <- "provide both 'sd_mean' and 'sd_var'"
    test_get_alpha_beta <- function(x, y){
        get_alpha_beta(alpha=NULL, beta=NULL, sd_mean=x, sd_var=y, p=2)
    }
    expect_error(test_get_alpha_beta(NULL, c(2, 1)), err_string)
    expect_error(test_get_alpha_beta(c(1, 2), NULL), err_string)
    expect_error(test_get_alpha_beta(NULL, NULL),    "provide both[a-z ']*OR")
})

test_that("get_alpha_beta returns a list", {
    expect_that(get_alpha_beta(alpha=c(1, 2), beta=c(1, 2), p=2), is_a("list"))
    expect_that(get_alpha_beta(alpha=NULL, beta=NULL, sd_mean=c(1, 2),
                               sd_var=c(1, 2), p=2), is_a("list"))
})

test_that("get_alpha_beta return value has names alpha and beta", {
    test_alpha_names <- function(x, y, p){
        sort(names(get_alpha_beta(alpha=x, beta=y, p=p)))
    }
    test_sd_names <- function(x, y, p){
        sort(names(get_alpha_beta(alpha=NULL, beta=NULL, sd_mean=x, sd_var=y,
                                  p=p)))
    }
    expected_names <- c("alpha", "beta")
    expect_equal(test_alpha_names(c(1, 2), c(1, 2), 2), expected_names)
    expect_equal(test_sd_names(c(1, 2), c(1, 2), 2),    expected_names)
})

test_that("get_alpha_beta returns alpha and beta", {
    a <- c(1, 2)
    b  <- c(2, 3)
    expect_equal(get_alpha_beta(alpha=a, beta=b, p=2)$alpha, a)
    expect_equal(get_alpha_beta(alpha=a, beta=b, p=2)$beta,  b)
})

test_that("get_alpha_beta returns correct transformation of sd_mean, sd_var", {
    gamma_mean <- function(alpha, beta) alpha/beta
    gamma_var  <- function(alpha, beta) alpha/(beta^2)
    mean_val <- c(1, 2)
    var_val  <- c(2, 3)
    ab <- get_alpha_beta(alpha=NULL, beta=NULL, sd_mean=mean_val,
                         sd_var=var_val, p=2)
    expect_equal(gamma_mean(ab$alpha, ab$beta), mean_val)
    expect_equal(gamma_var(ab$alpha, ab$beta),  var_val)
})

test_that("get_eta gives error if eta < 1", {
    err_string <- "must be >= 1"
    expect_error(get_eta(0),  err_string)
    expect_error(get_eta(-1), err_string)
})

test_that("get_eta returns eta=1", {
    expect_equal(get_eta(1), 1)
})

test_that("get_eta returns eta if eta > 1", {
    expect_equal(get_eta(pi), pi)
    expect_equal(get_eta(2), 2)
})

