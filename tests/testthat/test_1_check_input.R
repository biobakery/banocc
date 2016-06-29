context("Step 1 - checking input")

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

test_that("check_n works with all lengths", {
    test_check_n <- function(x, p){
        suppressWarnings(check_n(n=x, p=p))
    }
    expect_equal(test_check_n(1:10, 5), 1:5)
    expect_equal(test_check_n(1:10, 10), 1:10)
    expect_equal(test_check_n(1:10, 12), c(1:10, 1:2))
})

test_that("check_L fails if L not numeric", {
    err_string <- "must be numeric"
    expect_error(check_L(L="string",   p=1), err_string)
    expect_error(check_L(L=list(1, 2), p=2), err_string)
})

test_that("check_L fails if L not pxp", {
    err_string <- "must be a square matrix with the same number of columns as C"
    expect_error(check_L(L=matrix(1:10, nrow=5), p=5), err_string)
    expect_error(check_L(L=matrix(1:16, nrow=4), p=5), err_string)
})

test_that("check_L fails if L not positive definite", {
    err_string <- "not positive definite"
    expect_error(check_L(L=matrix(-16:-1, nrow=4), p=4), err_string)
})

test_that("check_L fails if L not symmetric", {
    err_string <- "must be symmetric"
    L_test <- diag(3)
    L_test[1, 2] <- 0.1
    expect_error(check_L(L=L_test, p=3), err_string)
})

test_that("check_L warns if L is a vector with mismatched length", {
    expect_warning(check_L(L=1:10, p=5),  "only first p elements")
    expect_warning(check_L(L=1:10, p=12), "recycled")
})

test_that("check_L returns a matrix", {
    test_check_L <- function(x, p) {
        suppressWarnings(check_L(L=x, p=p))
    }
    l_vec <- 1:5
    expected_class <- is_a("matrix")
    expect_that(test_check_L(diag(l_vec), 5),  expected_class)
    expect_that(test_check_L(l_vec, 5),        expected_class)
    expect_that(test_check_L(l_vec, 3),        expected_class)
    expect_that(test_check_L(l_vec, 7),        expected_class)

})

test_that("check_L works overall", {
    test_check_L <- function(x, p){
        suppressWarnings(check_L(L=x, p=p))
    }
    l <- 1:5
    l_long <- c(l, l[1:2])
    l_short <- l[1:3]
    expect_equal(test_check_L(diag(l), length(l)), diag(l))
    expect_equal(test_check_L(l, length(l)),       diag(l))
    expect_equal(test_check_L(l, length(l_short)), diag(l_short))
    expect_equal(test_check_L(l, length(l_long)),  diag(l_long))
})

test_that("check_a_b fails if a and b have unequal length", {
    err_string <- "must be of equal length"
    expect_error(check_a_b(a=1, b=c(1, 2), p=1), err_string)
    expect_error(check_a_b(a=c(1, 2), b=1, p=1), err_string)
    expect_error(check_a_b(a=1, b=c(1, 2), p=2), err_string)
    expect_error(check_a_b(a=c(1, 2), b=1, p=2), err_string)
})

test_that("check_a_b fails if a or b have negative values", {
    err_string <- "values must be positive"
    expect_error(check_a_b(a=c(1, -1), b=1:2, p=2), err_string)
    expect_error(check_a_b(a=c(-1, 1), b=1:2, p=2), err_string)
    expect_error(check_a_b(a=1:2, b=c(1, -1), p=2), err_string)
    expect_error(check_a_b(a=1:2, b=c(-1, 1), p=2), err_string)
})

test_that("check_a_b fails if a or b have zero values", {
    err_string <- "values must be positive"
    expect_error(check_a_b(a=c(0, 1), b=1:2, p=2),  err_string)
    expect_error(check_a_b(a=c(1, 0), b=1:2, p=2),  err_string)
    expect_error(check_a_b(a=1:2, b=c(0, 1), p=2),  err_string)
    expect_error(check_a_b(a=1:2, b=c(1, 0), p=2),  err_string)
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
