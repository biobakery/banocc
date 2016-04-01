context("Checking input")

test_that("check_vector of non-numeric vector gives error", {
    test_check_vector <- function(x){
        check_vector(parm.name="test_vector", parm=x, p=length(x),
                     verbose=FALSE)
    }
    err_string <- "must be a vector"
    expect_error(test_check_vector("string"),             err_string)
    expect_error(test_check_vector(matrix(1:10, nrow=5)), err_string)
    expect_error(test_check_vector(list(12, 13)),         err_string)
})

test_that("check_vector gives warning with mismatched lengths",{
    test_check_vector <- function(x, p){
        check_vector(parm.name="test_vector", parm=x, p=p, verbose=FALSE)
    }
    expect_warning(test_check_vector(1:10, 5), "using first p elements")
    expect_warning(test_check_vector(1:10, 12), "recycling")
})

test_that("check_vector recycles with too-short lengths", {
    test_check_vector <- function(x, p){
        suppressWarnings(check_vector(parm.name="test_vector", parm=x, p=p,
                                      verbose=FALSE))
    }
    expect_equal(test_check_vector(1:10, 5), 1:5)
    expect_equal(test_check_vector(1:10, 12), c(1:10, 1:2))
})

test_that("check_nu works with all lengths", {
    test_check_nu <- function(x, p){
        suppressWarnings(check_nu(nu=x, p=p))
    }
    expect_equal(test_check_nu(1:10, 5), 1:5)
    expect_equal(test_check_nu(1:10, 10), 1:10)
    expect_equal(test_check_nu(1:10, 12), c(1:10, 1:2))
})

test_that("check_Lambda fails if Lambda not numeric", {
    err_string <- "must be numeric"
    expect_error(check_Lambda(Lambda="string",   p=1), err_string)
    expect_error(check_Lambda(Lambda=list(1, 2), p=2), err_string)
})

test_that("check_Lambda fails if Lambda not pxp", {
    err_string <- "must be a square matrix with the same number of columns as C"
    expect_error(check_Lambda(Lambda=matrix(1:10, nrow=5), p=5), err_string)
    expect_error(check_Lambda(Lambda=matrix(1:16, nrow=4), p=5), err_string)
})

test_that("check_Lambda fails if Lambda not positive definite", {
    err_string <- "not positive definite"
    expect_error(check_Lambda(Lambda=matrix(-16:-1, nrow=4), p=4), err_string)
})

test_that("check_Lambda fails if Lambda not symmetric", {
    err_string <- "must be symmetric"
    Lambda_test <- diag(3)
    Lambda_test[1, 2] <- 0.1
    expect_error(check_Lambda(Lambda=Lambda_test, p=3), err_string)
})

test_that("check_Lambda warns if Lambda is a vector with mismatched length", {
    expect_warning(check_Lambda(Lambda=1:10, p=5),  "only first p elements")
    expect_warning(check_Lambda(Lambda=1:10, p=12), "recycled")
})

test_that("check_Lamba returns a matrix", {
    test_check_Lambda <- function(x, p) {
        suppressWarnings(check_Lambda(Lambda=x, p=p))
    }
    lambda_vec <- 1:5
    expected_class <- is_a("matrix")
    expect_that(test_check_Lambda(diag(lambda_vec), 5),  expected_class)
    expect_that(test_check_Lambda(lambda_vec, 5),        expected_class)
    expect_that(test_check_Lambda(lambda_vec, 3),        expected_class)
    expect_that(test_check_Lambda(lambda_vec, 7),        expected_class)

})

test_that("check_Lambda works overall", {
    test_check_Lambda <- function(x, p){
        suppressWarnings(check_Lambda(Lambda=x, p=p))
    }
    l <- 1:5
    l_long <- c(l, l[1:2])
    l_short <- l[1:3]
    expect_equal(test_check_Lambda(diag(l), length(l)), diag(l))
    expect_equal(test_check_Lambda(l, length(l)),       diag(l))
    expect_equal(test_check_Lambda(l, length(l_short)), diag(l_short))
    expect_equal(test_check_Lambda(l, length(l_long)),  diag(l_long))
})

test_that("check_alpha_beta fails if alpha and beta have unequal length", {
    err_string <- "must be of equal length"
    expect_error(check_alpha_beta(alpha=1, beta=c(1, 2), p=1), err_string)
    expect_error(check_alpha_beta(alpha=c(1, 2), beta=1, p=1), err_string)
    expect_error(check_alpha_beta(alpha=1, beta=c(1, 2), p=2), err_string)
    expect_error(check_alpha_beta(alpha=c(1, 2), beta=1, p=2), err_string)
})

test_that("check_alpha_beta fails if alpha or beta have negative values", {
    err_string <- "values must be positive"
    expect_error(check_alpha_beta(alpha=c(1, -1), beta=1:2, p=2), err_string)
    expect_error(check_alpha_beta(alpha=c(-1, 1), beta=1:2, p=2), err_string)
    expect_error(check_alpha_beta(alpha=1:2, beta=c(1, -1), p=2), err_string)
    expect_error(check_alpha_beta(alpha=1:2, beta=c(-1, 1), p=2), err_string)
})

test_that("check_alpha_beta fails if alpha or beta have zero values", {
    err_string <- "values must be positive"
    expect_error(check_alpha_beta(alpha=c(0, 1), beta=1:2, p=2),  err_string)
    expect_error(check_alpha_beta(alpha=c(1, 0), beta=1:2, p=2),  err_string)
    expect_error(check_alpha_beta(alpha=1:2, beta=c(0, 1), p=2),  err_string)
    expect_error(check_alpha_beta(alpha=1:2, beta=c(1, 0), p=2),  err_string)
})

test_that("check_sd_mean_var fails if sd_mean, sd_var of unequal length", {
    err_string <- "must be of equal length"
    expect_error(check_sd_mean_var(sd_mean=1, sd_var=c(1, 2), p=1), err_string)
    expect_error(check_sd_mean_var(sd_mean=c(1, 2), sd_var=1, p=1), err_string)
    expect_error(check_sd_mean_var(sd_mean=1, sd_var=c(1, 2), p=2), err_string)
    expect_error(check_sd_mean_var(sd_mean=c(1, 2), sd_var=1, p=2), err_string)
})

test_that("check_sd_mean_var fails if sd_mean or sd_var have negative values", {
    err_string <- "values must be positive"
    expect_error(check_sd_mean_var(sd_mean=c(1, -1), sd_var=1:2, p=2), err_string)
    expect_error(check_sd_mean_var(sd_mean=c(-1, 1), sd_var=1:2, p=2), err_string)
    expect_error(check_sd_mean_var(sd_mean=1:2, sd_var=c(1, -1), p=2), err_string)
    expect_error(check_sd_mean_var(sd_mean=1:2, sd_var=c(-1, 1), p=2), err_string)
})

test_that("check_sd_mean_var fails if sd_mean or sd_var have zero values", {
    err_string <- "values must be positive"
    expect_error(check_sd_mean_var(sd_mean=c(0, 1), sd_var=1:2, p=2),  err_string)
    expect_error(check_sd_mean_var(sd_mean=c(1, 0), sd_var=1:2, p=2),  err_string)
    expect_error(check_sd_mean_var(sd_mean=1:2, sd_var=c(0, 1), p=2),  err_string)
    expect_error(check_sd_mean_var(sd_mean=1:2, sd_var=c(1, 0), p=2),  err_string)
})
