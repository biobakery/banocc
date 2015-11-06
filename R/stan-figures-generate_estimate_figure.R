#' Generate a heatmap figure comparing the posterior fits from a \code{stan}
#'   fit or optimization object.
#'
#' @param Fit Either a \code{stanfit} or a list with element \code{par} being
#'   a named vector of parameter estimates
#' @param A The dataset of true counts (absolute abundances)
#' @inheritParams get_posterior_quantiles
#' @inheritParams make_samples_list
#' @inheritParams sample_mu_prior
#' @name generate_estimate_figure
#' @importFrom rstan extract

#' @rdname generate_estimate_figure
generate_estimate_figure_MCMC <-
function(Fit, A, estimate_method="median", thin=1, verbose=FALSE){
    posterior.samples <- rstan::extract(Fit, permuted= FALSE)
    post.samples.list <- banocc::make_samples_list(posterior.samples, thin=thin,
                                                   concatenate.chains=TRUE)
    if (verbose) print(str(post.samples.list))
    estimates <-
            banocc::get_posterior_estimates(post.samples.list,
                                            estimate_method=estimate_method,
                                            parameter.names="ln_Rho")
    R.est <- estimates$ln_Rho

    categories <- c("Basis", paste0("Posterior ", estimate_method))
    banocc::compare_two_corrmat_figure(cor(A), R.est, categories)
}

#' @rdname generate_estimate_figure
generate_estimate_figure_optim <-
function(Fit, A, Data){
    R.est.vec  <- Fit$par[grep("^ln_Rho", names(Fit$par))]
    R.est      <- matrix(R.est.vec, ncol=Data$P)

    banocc::compare_two_corrmat_figure(cor(A), R.est,
                                       categories=c("Basis", "Posterior Mode"))
}
