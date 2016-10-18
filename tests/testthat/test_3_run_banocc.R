context("Step 3 - running BAnOCC")
load("testthat_objects/banocc_model_test.RData")
load("testthat_objects/sample_stan_data.RData")


test_that("get_lambda gives error if lambda <= 0", {
    err_string <- "must be > 0"
    expect_error(get_lambda(0),  err_string)
    expect_error(get_lambda(-1), err_string)
})

test_that("get_lambda returns lambda if lambda > 0", {
    expect_equal(get_lambda(pi), pi)
    expect_equal(get_lambda(2), 2)
})

n  <- rep(0, Data$P)
L  <- 10 * diag(Data$P)

init1    <- list(list(m = n, O=diag(Data$P)))
init2    <- NULL
init_opt <- list(init1, init2)

conf_alpha_opt    <- list(0.05, 0.5)
get_min_width_opt <- list(TRUE, FALSE)
calc_snc_opt      <- list(TRUE, FALSE)
use_matrix        <- list(TRUE, FALSE)
eval_convergence  <- list(TRUE, FALSE)

opt_idx <- expand.grid(init          = seq_along(init_opt),
                       conf_alpha    = seq_along(conf_alpha_opt),
                       get_min_width = seq_along(get_min_width_opt),
                       calc_snc      = seq_along(calc_snc_opt),
                       use_matrix    = seq_along(use_matrix),
                       eval_convergence = seq_along(eval_convergence))

all_args <- lapply(seq_len(nrow(opt_idx) * 2), function(i){
    j <- ifelse(i %% 2 == 0, i / 2, (i + 1) / 2)
    idx <- c(opt_idx$init[j],
             opt_idx$conf_alpha[j],
             opt_idx$get_min_width[j],
             opt_idx$calc_snc[j],
             opt_idx$use_matrix[j],
             opt_idx$eval_convergence[j])
    args <- list(init=init_opt[[idx[1]]],
                 conf_alpha=conf_alpha_opt[[idx[2]]],
                 get_min_width=get_min_width_opt[[idx[3]]],
                 calc_snc=calc_snc_opt[[idx[4]]],
                 use_matrix=use_matrix[[idx[5]]],
                 eval_convergence=use_matrix[[idx[6]]])
    return(args)
})

sw_run_banocc <- function(conf_alpha, get_min_width, calc_snc,
                          eval_convergence, use_matrix=FALSE, init=NULL){
    sink("banocc.out", type="output")
    if (use_matrix){
        data <- as.matrix(Data$C)
    } else {
        data <- as.data.frame(Data$C)
    }
    rb <- suppressWarnings(run_banocc(
        compiled_banocc_model=compiled_banocc_model, C=data, lambda=0.02,
        n=n, L=L,
        chains=1, iter=4, warmup=2, init=init,
        conf_alpha=conf_alpha, get_min_width=get_min_width,
        calc_snc=calc_snc, eval_convergence=eval_convergence,
        verbose=FALSE, num_level=0
        ))
    sink()
    return(rb)
}

on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
on_bioc <- !identical(Sys.getenv("BBS_HOME"), "")
if (!on_cran && !on_bioc){
    rb <- lapply(seq_along(all_args), function(i){
        do.call(what=sw_run_banocc, args=all_args[[i]])
    })
}

test_that("run_banocc returns a list", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_is(rb[[i]], "list")
    }
})

test_that("run_banocc takes a matrix", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$use_matrix){
            expect_is(rb[[i]], "list")
        }
    }
})

test_that("run_banocc takes a data frame", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if(!all_args[[i]]$use_matrix){
            expect_is(rb[[i]], "list")
        }
    }
})

rb_names <- c("Data", "Fit", "CI.hpd", "Estimates.median")
extra_names <- c("Min.width", "SNC")

test_that("run_banocc returns list with correct names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        names_i <- rb_names
        if (all_args[[i]]$get_min_width) names_i <- c(names_i, extra_names[1])
        if (all_args[[i]]$calc_snc) names_i <- c(names_i, extra_names[2])
        expect_equal(sort(names(rb[[i]])), sort(names_i))
    }
})

test_that("run_banocc Data elt is list", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_is(rb[[i]]$Data, "list")
    }
})

data_names <- c("C", "N", "P", "n", "L", "lambda")
test_that("run_banocc Data elt has correct names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_equal(sort(names(rb[[i]]$Data)), sort(data_names))
    }
})

test_that("run_banocc Fit elt is stanfit object", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_is(rb[[i]]$Fit, "stanfit")
    }
})

test_that("run_banocc CI.hpd elt is list with correct names", {
    skip_on_cran()
    skip_on_bioc()
    ci_names <- c("lower", "upper")
    for (i in seq_along(rb)){
        expect_is(rb[[i]]$CI.hpd, "list")
        expect_equal(sort(names(rb[[i]]$CI.hpd)), sort(ci_names))
    }
})

p <- ncol(Data$C)
test_that("run_banocc CI.hpd elt elements are pxp matrices", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        for (k in seq_along(rb[[i]]$CI.hpd)){
            expect_is(rb[[i]]$CI.hpd[[k]], "matrix")
            expect_equal(dim(rb[[i]]$CI.hpd[[k]]), c(p, p))
        }
    }
})

test_that("run_banocc CI.hpd elt elements have col and row names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        for (k in seq_along(rb[[i]]$CI.hpd)){
            ci <- rb[[i]]$CI.hpd[[k]]
            expect_equal(colnames(ci), names(Data$C))
            expect_equal(rownames(ci), names(Data$C))
        }
    }
})

test_that("run_banocc CI.hpd elt elements are between -1 and 1", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        for (k in seq_along(rb[[i]]$CI.hpd)){
            if (!all_args[[i]]$eval_convergence){
                expect_true(all(rb[[i]]$CI.hpd[[k]] - 1 <= 1e-12))
                expect_true(all(rb[[i]]$CI.hpd[[k]] + 1 >= 1e-12))
            } else {
                expect_true(all(is.na(rb[[i]]$CI.hpd[[k]])))
            }
        }
    }
})

test_that("run_banocc Estimates.median elt is a pxp matrix", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_is(rb[[i]]$Estimates.median, "matrix")
        expect_equal(dim(rb[[i]]$Estimates.median), c(p, p))
    }
})

test_that("run_banocc Estimates.median has col and row names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        est.med <- rb[[i]]$Estimates.median
        expect_equal(colnames(est.med), names(Data$C))
        expect_equal(rownames(est.med), names(Data$C))
    }
})

test_that("run_banocc Estimates.median elts are between -1 and 1", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        est.med <- rb[[i]]$Estimates.median
        if (!all_args[[i]]$eval_convergence){
            expect_true(all(est.med + 1 >= 1e-12))
            expect_true(all(est.med - 1 <= 1e-12))
        } else {
            expect_true(all(is.na(est.med)))
        }
    }
})

test_that("run_banocc Min.width elt is a pxp matrix", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$get_min_width){
            expect_is(rb[[i]]$Min.width, "matrix")
            expect_equal(dim(rb[[i]]$Min.width), c(p, p))
        }
    }
})

test_that("run_banocc Min.width has col and row names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$get_min_width){
            mw <- rb[[i]]$Min.width
            expect_equal(colnames(mw), names(Data$C))
            expect_equal(rownames(mw), names(Data$C))
        }
    }
})

test_that("run_banocc Min.width elts are between 0 and 1", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$get_min_width){
            mw <- rb[[i]]$Min.width
            if (!all_args[[i]]$eval_convergence){
                expect_true(all(mw - 1<= 1e-12))
                expect_true(all(mw    >= 1e-12))
            } else {
                expect_true(all(is.na(mw)))
            }
        }
    }
})

test_that("run_banocc SNC elt is a pxp matrix", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$calc_snc){
            expect_is(rb[[i]]$SNC, "matrix")
            expect_equal(dim(rb[[i]]$SNC), c(p, p))
        }
    }
})

test_that("run_banocc SNC has col and row names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$calc_snc){
            snc <- rb[[i]]$SNC
            expect_equal(colnames(snc), names(Data$C))
            expect_equal(rownames(snc), names(Data$C))
        }
    }
})

test_that("run_banocc SNC elts are between 0 and 1", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$calc_snc){
            snc <- rb[[i]]$SNC
            if (!all_args[[i]]$eval_convergence){
                if (any(!is.na(snc))){
                    snc <- na.omit(as.vector(snc))
                    expect_true(all(snc - 1 <= 1e-12))
                    expect_true(all(snc    >= 1e-12))
                }
            } else {
                expect_true(all(is.na(snc)))
            }
        }
    }
})

## Need to test that banocc takes: 1) A data frame or 2) a matrix with 3) compositional rows
## Need to test that C can have zeros
