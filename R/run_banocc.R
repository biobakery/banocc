#' Runs BAnOCC to fit the model and generate appropriate convergence metrics
#'   and inference.
#'
#' @param bayes_model The compiled stan model (as with
#'   \code{stan_model(model_code = bayesStanModel)}).
#' @param C The dataset (rows=samples, columns=features), which is nxp
#' @param nu The prior mean for mu; vectors of length less than p will be
#'   recycled.
#' @param Lambda The prior variance-covariance for mu (must be
#'   positive-definite with dimension pxp where p=number of features in C), or
#'   a vector of length p of variances for mu. If a vector of length less than p
#'   is given, it will be recycled.
#' @param alpha,beta Vectors of parameters for the gamma prior on the standard
#'   deviations. They must be of equal length. If a vector of length less than p
#'   is given, it will be recycled.
#' @param sd_mean,sd_var Vectors of means and variances for the gamma priors on
#'   the standard deviations. They must be of equal length. If a vector of length
#'   less than p is given, it will be recycled.
#' @param eta scale parameter of LKJ distribution; must be >= 1.
#' @param init The initial values as a list (see
#'   \code{\link[rstan]{sampling}}). Default value is NULL, which means that
#'   initial values are sampled from the priors.
#' @inheritParams rstan::sampling
#' @param get_min_width A boolean value: should the minimum CI width that
#'   includes zero be calculated?
#' @name run_banocc
#' @param verbose Print informative statements as the function executes?
#' @param num_level The number of the level (will determine the number of
#'   spaces to add to verbose output)
#'
#'
#' @export

#' @rdname run_banocc
run_banocc <- function(bayes_model, C, nu = rep(0, ncol(C)),
                       Lambda = 10*diag(ncol(C)),
                       alpha = rep(1, ncol(C)), beta = rep(0.5, ncol(C)),
                       eta = 1, cores = getOption("mc.cores", 1L),
                       chains = 4, iter = 2000, warmup = floor(iter/2),
                       thin = 1, init = NULL, control=NULL,
                       sd_mean=NULL, sd_var=NULL, conf_alpha=0.05,
                       get_min_width=FALSE, verbose=FALSE, num_level=0){
    cat_v("Begin run_banocc\n", verbose, num_level=num_level)
    Data <- list(C=C, N=nrow(C), P=ncol(C))
    
    Data$nu     <- banocc::check_nu(nu, Data$P, verbose, num_level=num_level+1)
    Data$Lambda <- banocc::check_Lambda(Lambda, Data$P, verbose,
                                        num_level=num_level+1)
    alpha_beta <- banocc::get_alpha_beta(alpha=alpha, beta=beta,
                                         sd_mean=sd_mean, sd_var=sd_var,
                                         p=Data$P,
                                         verbose=verbose, num_level=num_level+1)
    Data$alpha <- alpha_beta$alpha
    Data$beta <- alpha_beta$beta
    Data$eta <- banocc::get_eta(eta)

    if (is.null(init)){
        init <- banocc::get_IVs(chains=chains, data=Data, verbose=verbose,
                                num_level=num_level + 1)
    }

    cat_v("Begin fitting the model\n", verbose, num_level=num_level+1)
    Fit.all <- banocc::mycapture(rstan::sampling(bayes_model, data=Data,
                                                 chains=chains, iter=iter,
                                                 warmup=warmup, thin=thin,
                                                 init=init, cores=cores,
                                                 control=control))
    Fit <- Fit.all$output
    cat_v("End fitting the model\n", verbose, num_level=num_level+1)

    post.samples.list <- rstan::extract(Fit)
    CI <- banocc::get_credible_intervals(posterior_samples=post.samples.list,
                                         list=TRUE,
                                         parameter.names=c("ln_Rho"),
                                         conf=1-conf_alpha,
                                         type="marginal.hpd",
                                         verbose=verbose, num_level=num_level+1)

    dimnames(CI$ln_Rho$lower) <- list(colnames(Data$C), colnames(Data$C))
    dimnames(CI$ln_Rho$upper) <- list(colnames(Data$C), colnames(Data$C))
    CI <- CI$ln_Rho

    if (get_min_width){
        min_width <- banocc::get_min_width(posterior_sample=post.samples.list,
                                           parameter.names=c("ln_Rho"),
                                           null_value=0, type="marginal.hpd",
                                           precision=0.01, verbose=verbose,
                                           num_level=num_level + 1)
    } else {
        min_width <- list(ln_Rho=NULL)
    }
    min_width <- min_width$ln_Rho

    Estimates <-
        banocc::get_posterior_estimates(posterior_samples=post.samples.list,
                                        estimate_method="median",
                                        parameter.names="ln_Rho")
    dimnames(Estimates$ln_Rho) <- list(colnames(Data$C), colnames(Data$C))
    Estimates <- Estimates$ln_Rho

    
    return_object <- list(Data=Data, Fit=Fit, Fit.print=Fit.all$print.output,
                          CI.hpd=CI, Estimates.median=Estimates,
                          Min.width=min_width)

    cat_v("End run_banocc\n", verbose, num_level=num_level)

    return(return_object)
}

check_nu <- function(nu, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_nu\n", verbose, num_level=num_level)
    nu <- as.numeric(nu)
    nu <- banocc::check_vector("nu", nu, p, verbose, num_level=num_level+1)
    cat_v("End check_nu\n", verbose, num_level=num_level)
    return(nu)
}

check_alpha_beta <- function(alpha, beta, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_alpha_beta\n", verbose, num_level=num_level)
    if (length(alpha) != length(beta)){
        stop("'alpha' and 'beta' must be of equal length")
    }
    alpha <- as.numeric(alpha)
    alpha <- banocc::check_vector("alpha", alpha, p, verbose,
                                  num_level=num_level + 1)
    if (any(alpha <= 0)){
        stop("'alpha' values must be positive")
    }
    beta  <- as.numeric(beta)
    beta  <- banocc::check_vector("beta", beta, p, verbose,
                                  num_level=num_level + 1)
    if (any(beta <= 0)){
        stop("'beta' values must be positive")
    }
    cat_v("End check_alpha_beta\n", verbose, num_level=num_level)
    return(list(alpha=alpha, beta=beta))
}

check_sd_mean_var <- function(sd_mean, sd_var, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_sd_mean_var\n", verbose, num_level=num_level)
    if (length(sd_mean) != length(sd_var)){
        stop("'sd_mean' and 'sd_var' must be of equal length")
    }
    sd_mean <- as.numeric(sd_mean)
    sd_mean <- banocc::check_vector("sd_mean", sd_mean, p, verbose,
                                    num_level=num_level + 1)
    if (any(sd_mean <= 0)){
        stop("'sd_mean' values must be positive")
    }
    sd_var  <- as.numeric(sd_var)
    sd_var  <- banocc::check_vector("sd_var", sd_var, p, verbose,
                                    num_level = num_level + 1)
    if (any(sd_var <= 0)){
        stop("'sd_var' values must be positive")
    }
    cat_v("End check_sd_mean_var\n", verbose, num_level=num_level)
    return(list(sd_mean=sd_mean, sd_var=sd_var))
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

check_Lambda <- function(Lambda, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_Lambda...", verbose, num_level=num_level)
    if (!is.numeric(Lambda)){
        stop("Lambda must be numeric")
    }

    if (is.matrix(Lambda)){
        if (any(dim(Lambda) != p)){
            stop(paste0("Lambda must be a square matrix with the same number of",
                        " columns as C"))
        }
        if (any(eigen(Lambda)$values <= 0)){
            stop("'Lambda' is not positive definite.")
        }
        if (!is.symmetric(Lambda)){
            stop("'Lambda' must be symmetric")
        }
    } else if (is.vector(Lambda)){
        if ((length(Lambda) < p) && (length(Lambda) > 1)){
            warning("'Lambda' is being recycled")
            num_rep <- ceiling(p / length(Lambda))
            Lambda <- diag(rep(Lambda, num_rep)[1:p])
        } else if (length(Lambda) == 1){
            Lambda <- diag(rep(Lambda, p))
        } else if (length(Lambda) > p){
            warning("'Lambda' has length > p; only first p elements will be ",
                    "used.")
            Lambda <- diag(Lambda[1:p])
        } else {
            Lambda <- diag(Lambda)
        }
    } else {
        stop("'Lambda' must be either a matrix or a vector")
    }
    cat_v("Done.\n", verbose)
    return(Lambda)
}

get_alpha_beta <- function(alpha, beta, p, sd_mean=NULL, sd_var=NULL,
                           verbose=FALSE, num_level=0){
    cat_v("Begin get_alpha_beta...", verbose=verbose,
                  num_level=num_level)
    if (!is.null(alpha) || !is.null(beta)){
        if (is.null(alpha) || is.null(beta)){
            stop(paste0("Must provide both 'alpha' and 'beta' OR both 'sd_mean'",
                        " and 'sd_var'"))
        }
        alpha_beta <- banocc::check_alpha_beta(alpha, beta, p, verbose,
                                               num_level=num_level+1)
    } else if (!is.null(sd_var) || !is.null(sd_mean)){
        if (is.null(sd_var) || is.null(sd_mean)){
            stop(paste0("Must provide both 'sd_mean' and 'sd_var' OR both ",
                        "'alpha' and 'beta'"))
        }
        alpha_beta <- list()
        sd_mean_var <- banocc::check_sd_mean_var(sd_mean, sd_var, p,
                                                 num_level=num_level+1)
        alpha_beta$beta  <- sd_mean_var$sd_mean / sd_mean_var$sd_var
        alpha_beta$alpha <- sd_mean_var$sd_mean^2 / sd_mean_var$sd_var
    } else {
        stop(paste0("Must provide both 'alpha' and 'beta' OR both 'sd_mean'",
                    " and 'sd_var'"))
    }
    cat_v("Done.\n", verbose)
    return(alpha_beta)
}

get_eta <- function(eta){
    if (eta < 1){
        stop("'eta' must be >= 1")
    } else {
        return(eta)
    }

}
