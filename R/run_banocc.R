#' Runs BAnOCC to fit the model and generate appropriate convergence metrics
#'   and inference.
#'
#' @param compiled_banocc_model The compiled stan model (as with
#'   \code{stan_model(model_code = banocc_model)}).
#' @param C The dataset as a data frame or matrix. This should be N by P
#'   with N samples as the rows and P features as the columns. 
#' @param n The prior mean for m; vectors of length less than P (the number
#'   of features/columns of \code{C}) will be recycled.
#' @param L The prior variance-covariance for m (must be
#'   positive-definite with dimension PxP where P=number of features/columns
#'   in \code{C}), or a vector of length p of variances for m. If a vector of
#'   length less than P is given, it will be recycled.
#' @param a The shape parameter of a gamma distribution (the prior on the
#'   shrinkage parameter lambda)
#' @param b The scale parameter of a gamma distribution (the prior on the
#'   shrinkage parameter lambda)
#' @param init The initial values as a list (see
#'   \code{\link[rstan]{sampling}} in the \code{rstan} package). Default
#'   value is NULL, which means that initial values are sampled from the
#'   priors for parameters m and lambda while O is set to the identity matrix.
#' @inheritParams rstan::sampling
#' @param get_min_width A boolean value: should the minimum CI width that
#'   includes zero be calculated?
#' @param conf_alpha The percentage of the posterior density outside the
#'   credible interval. That is, a \code{1-conf_alpha} * 100\% credible
#'   interval will be returned.
#' @param calc_snc Boolean: should the scaled neighborhood criterion be
#'   calculated?
#' @param eval_convergence Boolean: if `TRUE`, convergence will be evaluated
#'   using the Rhat statistic, and the fit output (estimates, credible
#'   intervals, etc.) will be missing if this statistic does not indicate
#'   convergence.
#' @param verbose Print informative statements as the function executes?
#' @param num_level The number of indentations to add to the output when
#'   \code{verbose = TRUE}.
#'
#' @importFrom rstan extract
#' @importFrom rstan sampling
#' @importFrom rstan summary
#' @importFrom stringr str_match
#'
#' @export
#'
#' @return
#' Returns a named list with the following elements:
#' \describe{
#'   \item{\emph{Data}}{The data formatted as a named list that includes the
#'     input data (\code{C}) and the prior parameters (\code{n}, \code{L},
#'     \code{a}, \code{b})}
#' 
#'   \item{\emph{Fit}}{The \code{stanfit} object returned by the call to
#'     \code{\link[rstan]{sampling}}}
#' 
#'   \item{\emph{CI}}{The \code{1-conf_alpha} * 100\% credible intervals}
#'
#'   \item{\emph{Estimates.median}}{The correlation estimates, which are the
#'     marginal posterior medians}
#' 
#'   \item{\emph{Min.width}}{Only present if the \code{get_min_width}
#'     argument is \code{TRUE}. The minimum CI width that includes zero for
#'     each correlation.}
#' 
#'   \item{\emph{SNC}}{Only present if the \code{calc_snc} argument is
#'     \code{TRUE}. The scaled neighborhood criterion for each correlation.}
#' }
#' 
#' @examples
#'   data(compositions_null)
#'   \dontrun{
#'     compiled_banocc_model <- rstan::stan_model(model_code=banocc_model)
#'     b_output <- run_banocc(C=compositions_null,
#'                            compiled_banocc_model=compiled_banocc_model)
#'   }
#'
#' @seealso \code{vignette("banocc-vignette")} for more examples.

run_banocc <- function(compiled_banocc_model, C, n = rep(0, ncol(C)),
                       L = 10*diag(ncol(C)), a=0.5, b=0.01,
                       cores = getOption("mc.cores", 1L),
                       chains = 4, iter = 50, warmup = floor(iter/2),
                       thin = 1, init = NULL, control=NULL,
                       sd_mean=NULL, sd_var=NULL, conf_alpha=0.05,
                       get_min_width=FALSE, calc_snc=FALSE,
                       eval_convergence=TRUE, verbose=FALSE, num_level=0){
    cat_v("Begin run_banocc\n", verbose, num_level=num_level)
    C <- check_C(C, verbose=verbose, num_level=num_level+1)
    Data <- list(C=C, N=nrow(C), P=ncol(C))
    
    Data$n <- check_n(n, Data$P, verbose, num_level=num_level+1)
    Data$L <- check_L(L, Data$P, verbose,
                      num_level=num_level+1)

    conf_alpha <- check_conf_alpha(conf_alpha, verbose,
                                   num_level=num_level+1)
    Data$a <- get_gamma_param(param=a, name="a")
    Data$b <- get_gamma_param(param=b, name="b")

    if (is.null(init)){
        init <- get_IVs(chains=chains, data=Data, verbose=verbose,
                        num_level=num_level + 1)
        test_output <- rstan::sampling(
            compiled_banocc_model, data=Data, init=init, chains=chains,
            iter=4, warmup=2, refresh=0, show_messages=FALSE)
        num_tests <- 1
        while (length(dimnames(test_output)) == 0 && num_tests < 10){
            init <- get_IVs(chains=chains, data=Data, verbose=verbose,
                            num_level=num_level + 1)
            test_output <- rstan::sampling(
                compiled_banocc_model, data=Data, init=init, chains=chains,
                iter=4, warmup=2, refresh=0, show_messages=FALSE)
            num_tests <- num_tests + 1
        }
        if (length(dimnames(test_output)) == 0){
            stop(paste0("Unable to generate workable starting values from",
                        " priors after 10 tries. Try specifying the values ",
                        "by hand."))
        }
    }

    cat_v("Begin fitting the model\n", verbose, num_level=num_level+1)
    refresh <- ifelse(verbose, max(iter/10, 1), 0)
    Fit <- rstan::sampling(compiled_banocc_model, data=Data,
                           chains=chains, iter=iter,
                           warmup=warmup, thin=thin,
                           init=init, cores=cores,
                           control=control,
                           show_messages=FALSE,
                           refresh=refresh)
    cat_v("End fitting the model\n", verbose, num_level=num_level+1)

    if (eval_convergence){
        cat_v("Begin evaluating convergence\n", verbose,
              num_level=num_level+1)
        rhat_stat <- rstan::summary(Fit)$summary[, "Rhat"]
        diag_elts <- grep("W.*\\[([0-9]*),[ ]?\\1\\]", names(rhat_stat))
        rhat_stat <- rhat_stat[-diag_elts]
        ## chol_ut_re <- "[Cc]hol\\[([0-9]+),[ ]?([0-9]+)"
        ## chol_ut_idx <- stringr::str_match(names(rhat_stat),
        ##                                   chol_ut_re)[, c(2, 3)]
                                                 
        ## chol_upper_tri_elts <- apply(chol_ut_idx, 1, function(r){
        ##     ifelse(all(is.na(r)), FALSE, which.max(as.numeric(r)) == 2)
        ## })
        ## rhat_stat <- rhat_stat[-which(chol_upper_tri_elts)]
        if (any(is.na(rhat_stat)) || max(rhat_stat) > 1.2){
            fit_converged <- FALSE
            warning(paste0("Fit has not converged as evaluated by the Rhat ",
                           "statistic. You might try a larger number of ",
                           "warmup iterations, different priors, or ",
                           "different initial values. See vignette for ",
                           "more on evaluating convergence."))
        } else {
            fit_converged <- TRUE
        }
        cat_v("End evaluating convergence\n", verbose, num_level=num_level+1)
    } else {
        fit_converged <- TRUE
    }

    post_samples_list <- rstan::extract(Fit)
    if (fit_converged){
        CI <- get_credible_intervals(posterior_samples=post_samples_list,
                                     list=TRUE, parameter.names=c("W"),
                                     conf=1-conf_alpha, type="marginal.hpd",
                                     verbose=verbose, num_level=num_level+1)
        
    } else {
        CI <- list(W=list(lower=matrix(NA, ncol=Data$P, nrow=Data$P),
                       upper=matrix(NA, ncol=Data$P, nrow=Data$P)))
    }
    dimnames(CI$W$lower) <- list(colnames(Data$C), colnames(Data$C))
    dimnames(CI$W$upper) <- list(colnames(Data$C), colnames(Data$C))
    CI <- CI$W

    if (fit_converged){
        Estimates <-
            get_posterior_estimates(posterior_samples=post_samples_list,
                                    estimate_method="median",
                                    parameter.names="W")
    } else {
        Estimates <- list(W=matrix(NA, ncol=Data$P, nrow=Data$P))
    }
    dimnames(Estimates$W) <- list(colnames(Data$C), colnames(Data$C))
    Estimates <- Estimates$W

    
    return_object <- list(Data=Data, Fit=Fit, 
                          CI.hpd=CI, Estimates.median=Estimates)

    if (get_min_width && fit_converged){
        min_width <- get_min_width(posterior_sample=post_samples_list,
                                   parameter.names=c("W"),
                                   null_value=0, type="marginal.hpd",
                                   precision=0.01, verbose=verbose,
                                   num_level=num_level + 1)
    } else if (get_min_width){
        min_width <- list(W=matrix(NA, ncol=Data$P, nrow=Data$P))
    }

    if (get_min_width){
        return_object$Min.width <- min_width$W
        colnames(return_object$Min.width) <- colnames(Data$C)
        rownames(return_object$Min.width) <- colnames(Data$C)
    }

    if (calc_snc && fit_converged){
        snc <- get_snc(posterior_samples=post_samples_list,
                       parameter.names=c("W"))
    } else if (calc_snc){
        snc <- list(W=matrix(NA, ncol=Data$P, nrow=Data$P))
    }

    if (calc_snc){
        return_object$SNC <- snc$W
        colnames(return_object$SNC) <- colnames(Data$C)
        rownames(return_object$SNC) <- colnames(Data$C)   
    }

    cat_v("End run_banocc\n", verbose, num_level=num_level)

    return(return_object)
}

check_n <- function(n, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_n\n", verbose, num_level=num_level)
    n <- as.numeric(n)
    n <- check_vector("n", n, p, verbose, num_level=num_level+1)
    cat_v("End check_n\n", verbose, num_level=num_level)
    return(n)
}


check_vector <- function(parm.name, parm, p, verbose=TRUE, num_level=0){
    cat_v("Begin check_vector...", verbose, num_level=num_level)
    if (!is.vector(parm) || mode(parm)!="numeric"){
        stop(paste0("'", parm.name, "' must be a vector"))
    }
    if ((length(parm) < p) && (length(parm) > 1)){
        warning(paste0("recycling '", parm.name, "'"))
        parm <- rep(parm, ceiling(p / length(parm)))[1:p]
    } else if (length(parm) == 1){
        parm <- rep(parm, p)
    } else if (length(parm) > p){
        warning("length of '", parm.name,
                "' is > p; only using first p elements")
        parm <- parm[1:p]
    }
    cat_v("Done.\n", verbose)
    return(parm)
}

check_L <- function(L, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_L...", verbose, num_level=num_level)
    if (!is.numeric(L)){
        stop("L must be numeric")
    }

    if (is.matrix(L)){
        if (any(dim(L) != p)){
            stop(paste0("L must be a square matrix with the same number of",
                        " columns as C"))
        }
        if (any(eigen(L)$values <= 0)){
            stop("'L' is not positive definite.")
        }
        if (!is.symmetric(L)){
            stop("'L' must be symmetric")
        }
    } else if (is.vector(L)){
        if ((length(L) < p) && (length(L) > 1)){
            warning("'L' is being recycled")
            num_rep <- ceiling(p / length(L))
            L <- diag(rep(L, num_rep)[1:p])
        } else if (length(L) == 1){
            L <- diag(rep(L, p))
        } else if (length(L) > p){
            warning("'L' has length > p; only first p elements will be ",
                    "used.")
            L <- diag(L[1:p])
        } else {
            L <- diag(L)
        }
    } else {
        stop("'L' must be either a matrix or a vector")
    }
    cat_v("Done.\n", verbose)
    return(L)
}

get_gamma_param <- function(param, name){
    if (param <= 0){
        stop(paste0("'", name, "' must be > 0"))
    } else {
        return(param)
    }

}

check_C <- function(C, zero_adj=0.0001, verbose=FALSE, num_level=0){
    cat_v("Begin checking input data matrix\n", verbose=verbose,
          num_level=num_level)
    if (any(C < 0)){
        stop("Some values of C are < 0")
    }
    if (any(C > 1)){
        stop("Some values of C are > 1")
    }
    if (any(rowSums(C) - 1 > 1e-8)){
        stop(paste0("Some row sums of C are not 1; perhaps you transposed ",
                    "features and samples?"))
    }
    if (!is.data.frame(C) && !is.matrix(C)){
        stop("C must be a data frame or matrix")
    }
    C <- as.matrix(C)
    if (!is.numeric(C)){
        stop("C must be numeric")
    }
    if (any(C == 0)){
        warning(paste0("Some values of C are zero. ",
                       "Since zero-inflation is not yet implemented, these ",
                       "will be changed to ", zero_adj, "."))
        C <- adjust_zeros(C, zero_adj=zero_adj, verbose=verbose,
                          num_level=num_level+1)
    }
    cat_v("End checking input data matrix\n", verbose=verbose,
          num_level=num_level)
    return(C)
}

# Check that conf_alpha is non-NULL and numeric
check_conf_alpha <- function(conf_alpha, verbose=FALSE, num_level=0){
    if (is.null(conf_alpha)){
        stop("conf_alpha must not be NULL")
    }
    conf_alpha_num <- suppressWarnings(as.numeric(conf_alpha))
    if (is.na(conf_alpha_num)){
        stop(paste0("conf_alpha must be coercible to numeric type. '",
                    conf_alpha, "' is not."))
    }
    return(conf_alpha_num)
}
