#' Takes a model fit from BAnOCC, evaluates convergence and generates
#'   appropriate convergence metrics and inference
#'
#' @param banoccfit Either a \code{stanfit} object (the \code{Fit} element
#'   returned by \code{run_banocc}), or the list returned by a call to
#'   \code{run_banocc}.
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
#' @inheritParams run_banocc
#'
#' @importFrom rstan extract
#' @importFrom rstan summary
#' @importFrom stringr str_to_lower
#' @importFrom stringr str_match
#' @export
#'
#' @return
#' Returns a named list with the following elements:
#' \describe{ 
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
#'
#'   \item{\emph{Fit}}{The \code{stanfit} object returned by the call to
#'     \code{run_banocc}.}
#' 
#'   \item{\emph{Data}}{Only present if the \code{banoccfit} argument is
#'     specified as the output of a call to \code{run_banocc}. It will be
#'     missing if \code{banoccfit} is specified as a \code{stanfit} object.}
#' }
#'
#' @examples
#' data(compositions_null)
#'   \dontrun{
#'     compiled_banocc_model <- rstan::stan_model(model_code=banocc_model)
#'     b_fit <- run_banocc(C=compositions_null,
#'                             compiled_banocc_model=compiled_banocc_model)
#'     b_output <- get_banocc_output(banoccfit=b_fit)
#'   }
#'
#' @seealso \code{vignette("banocc-vignette")} for more examples.

get_banocc_output <- function(banoccfit, conf_alpha=0.05, get_min_width=FALSE,
                              calc_snc=TRUE, eval_convergence=TRUE,
                              verbose=FALSE, num_level=0){
    cat_v("Begin get_banocc_output\n", verbose, num_level=num_level)
    b_stanfit <- get_stanfit(banoccfit)
    b_data    <- get_data(banoccfit)
    p <- ifelse(is.null(b_data),
                attr(b_stanfit, "par_dims")$m,
                b_data$P)
    conf_alpha <- check_conf_alpha(conf_alpha, verbose,
                                   num_level=num_level+1)


    if (eval_convergence){
        fit_converged <- evaluate_convergence(b_stanfit=b_stanfit,
                                              verbose=verbose,
                                              num_level=num_level + 1)
    } else {
        fit_converged <- TRUE
    }

    post_samples_list <- rstan::extract(b_stanfit)
    if (fit_converged){
        if (!('W' %in% names(post_samples_list))){
            post_samples_list$W <-
                aperm(array(apply(post_samples_list$O, 1, function(Oi){
                    cov2cor(solve(matrix(Oi, ncol=sqrt(length(Oi)))))
                }), dim=dim(post_samples_list$O)[c(3, 2, 1)]),
                      perm=c(3, 2, 1))
        }
        CI <- get_credible_intervals(posterior_samples=post_samples_list,
                                     list=TRUE, parameter.names=c("W"),
                                     conf=1-conf_alpha, type="marginal.hpd",
                                     verbose=verbose, num_level=num_level+1)
        
    } else {
        CI <- list(W=list(lower=matrix(NA, ncol=p, nrow=p),
                       upper=matrix(NA, ncol=p, nrow=p)))
    }
    dimnames(CI$W$lower) <- list(colnames(b_data$C), colnames(b_data$C))
    dimnames(CI$W$upper) <- list(colnames(b_data$C), colnames(b_data$C))
    CI <- CI$W

    if (fit_converged){
        Estimates <-
            get_posterior_estimates(posterior_samples=post_samples_list,
                                    estimate_method="median",
                                    parameter.names="W")
    } else {
        Estimates <- list(W=matrix(NA, ncol=p, nrow=p))
    }
    dimnames(Estimates$W) <- list(colnames(b_data$C), colnames(b_data$C))
    Estimates <- Estimates$W

    
    banocc_output <- list(Fit=b_stanfit, 
                          CI.hpd=CI, Estimates.median=Estimates)
    if (!is.null(b_data)) banocc_output$Data <- b_data

    if (get_min_width && fit_converged){
        min_width <- get_min_width(posterior_sample=post_samples_list,
                                   parameter.names=c("W"),
                                   null_value=0, type="marginal.hpd",
                                   precision=0.01, verbose=verbose,
                                   num_level=num_level + 1)
    } else if (get_min_width){
        min_width <- list(W=matrix(NA, ncol=p, nrow=p))
    }

    if (get_min_width){
        banocc_output$Min.width <- min_width$W
        colnames(banocc_output$Min.width) <- colnames(b_data$C)
        rownames(banocc_output$Min.width) <- colnames(b_data$C)
    }

    if (calc_snc && fit_converged){
        snc <- get_snc(posterior_samples=post_samples_list,
                       parameter.names=c("W"))
    } else if (calc_snc){
        snc <- list(W=matrix(NA, ncol=p, nrow=p))
    }

    if (calc_snc){
        banocc_output$SNC <- snc$W
        colnames(banocc_output$SNC) <- colnames(b_data$C)
        rownames(banocc_output$SNC) <- colnames(b_data$C)   
    }

    cat_v("End get_banocc_output\n", verbose, num_level=num_level)
    return(banocc_output)
}

get_stanfit <- function(banoccfit){
    if (methods::is(banoccfit, "stanfit")){
        return(banoccfit)
    } else if (methods::is(banoccfit, "list")){
        banoccfit_class <- unlist(lapply(banoccfit, methods::is, "stanfit"))
        if (any(banoccfit_class)){
            return(banoccfit[[which(banoccfit_class)]])
        } else {
            stop("No 'stanfit' object found in 'banoccfit' list.")
        }
    } else {
        stop("'banoccfit' must either have class 'stanfit' or 'list'")
    }
}

get_data <- function(banoccfit){
    if (methods::is(banoccfit, "list")){
        banoccfit_names <- stringr::str_to_lower(names(banoccfit))
        if ("data" %in% banoccfit_names){
            return(banoccfit[[which(banoccfit_names == "data")]])
        } else {
            warning("'banoccfit' given as a list, but no data ",
                    "element was found to return.")
            return(NULL)
        }
    } else {
        return(NULL)
    }
}

# Check that conf_alpha is non-NULL and numeric
check_conf_alpha <- function(conf_alpha, verbose=FALSE, num_level=0){
    if (is.null(conf_alpha)){
        stop("conf_alpha must not be NULL")
    }
    conf_alpha_num <- suppressWarnings(as.numeric(conf_alpha))
    if (is.na(conf_alpha_num)){
        stop("conf_alpha must be coercible to numeric type. '",
             conf_alpha, "' is not.")
    }
    return(conf_alpha_num)
}
