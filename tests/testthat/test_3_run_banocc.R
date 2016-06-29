context("Step 3 - running BAnOCC")
load("testthat_objects/banocc_model_test.RData")
load("testthat_objects/sample_stan_data.RData")

test_that("get_a_b gives error if a and b not provided", {
    err_string <- "provide both 'a' and 'b'"
    expect_error(get_a_b(a=c(1, 2), b=NULL, p=2),  err_string)
    expect_error(get_a_b(a=NULL, b=c(1, 2), p=2),  err_string)
    expect_error(get_a_b(a=NULL, b=NULL, p=2),     err_string)
})

test_that("get_a_b gives error if sd_mean and sd_var not both provided", {
    err_string <- "provide both 'sd_mean' and 'sd_var'"
    test_get_a_b <- function(x, y){
        get_a_b(a=NULL, b=NULL, sd_mean=x, sd_var=y, p=2)
    }
    expect_error(test_get_a_b(NULL, c(2, 1)), err_string)
    expect_error(test_get_a_b(c(1, 2), NULL), err_string)
    expect_error(test_get_a_b(NULL, NULL),    "provide both[a-z ']*OR")
})

test_that("get_a_b returns a list", {
    expect_that(get_a_b(a=c(1, 2), b=c(1, 2), p=2), is_a("list"))
    expect_that(get_a_b(a=NULL, b=NULL, sd_mean=c(1, 2),
                               sd_var=c(1, 2), p=2), is_a("list"))
})

test_that("get_a_b return value has names a and b", {
    test_a_names <- function(x, y, p){
        sort(names(get_a_b(a=x, b=y, p=p)))
    }
    test_sd_names <- function(x, y, p){
        sort(names(get_a_b(a=NULL, b=NULL, sd_mean=x, sd_var=y,
                                  p=p)))
    }
    expected_names <- c("a", "b")
    expect_equal(test_a_names(c(1, 2), c(1, 2), 2), expected_names)
    expect_equal(test_sd_names(c(1, 2), c(1, 2), 2),    expected_names)
})

test_that("get_a_b returns a and b", {
    a <- c(1, 2)
    b  <- c(2, 3)
    expect_equal(get_a_b(a=a, b=b, p=2)$a, a)
    expect_equal(get_a_b(a=a, b=b, p=2)$b,  b)
})

test_that("get_a_b returns correct transformation of sd_mean, sd_var", {
    gamma_mean <- function(a, b) a/b
    gamma_var  <- function(a, b) a/(b^2)
    mean_val <- c(1, 2)
    var_val  <- c(2, 3)
    ab <- get_a_b(a=NULL, b=NULL, sd_mean=mean_val,
                         sd_var=var_val, p=2)
    expect_equal(gamma_mean(ab$a, ab$b), mean_val)
    expect_equal(gamma_var(ab$a, ab$b),  var_val)
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

