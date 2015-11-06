#' Runs a stan workflow to fit the model and generate figures on SIMULATED DATA
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
#' @inheritParams rstan::sampling
#' @param Rho.true The true correlation matrix
#' @param workflow.name A string name of the workflow (to identify figures
#'   produced). If \code{NULL}, figures are returned but not output to a
#'   file.
#' @param n.prior The number of prior samples to use in plotting the prior
#'   densities.  If \code{NULL} or \code{0}, no priors are plotted.
#' @inheritParams cat_v
#' @name stan_workflow
#'
#' @importFrom rstan sampling
#' @importFrom rstan traceplot
#' @importFrom rstan extract

#' @rdname stan_workflow
stan_workflow <- function(bayes_model, C, nu, Lambda, alpha, beta,
                          eta,
                          cores = getOption("mc.cores", 1L),
                          chains = 4, iter = 2000, warmup = floor(iter/2),
                          thin = 1, init = 'random',
                          sd_mean=NULL, sd_var=NULL, n.prior=NULL,
                          workflow.name=NULL, verbose=FALSE, num_level=0,
                          conf_alpha=0.05){
    banocc::cat_v("Begin stan_workflow\n", verbose, num_level=num_level)
    Data <- list(C=C, N=nrow(C), P=ncol(C))
    
    Data$nu     <- banocc::check_nu(nu, Data$P, verbose, num_level=num_level+1)
    Data$Lambda <- banocc::check_Lambda(Lambda, Data$P, verbose,
                                        num_level=num_level+1)
    if (is.null(sd_mean) || is.null(sd_var)){
        if (is.null(alpha) || is.null(beta)){
            stop(paste0("Must provide both 'alpha' and 'beta' OR both 'sd_mean'",
                        " and 'sd_var'"))
        }
        alpha_beta <- banocc::check_alpha_beta(alpha, beta, Data$P, verbose,
                                               num_level=num_level+1)
        Data$alpha <- alpha_beta$alpha
        Data$beta  <- alpha_beta$beta
    }
    if (is.null(alpha) || is.null(beta)){
        if (is.null(sd_var) || is.null(sd_mean)){
            stop(paste0("Must provide both 'sd_mean' and 'sd_var' OR both ",
                        "'alpha' and 'beta'"))
        }
        sd_mean_var <- banocc::check_sd_mean_var(sd_mean, sd_var, Data$P,
                                                 num_level=num_level+1)
        Data$beta  <- sd_mean_var$sd_mean / sd_mean_var$sd_var
        Data$alpha <- sd_mean_var$sd_mean^2 / sd_mean_var$sd_var
    }
    if (eta < 1){
        stop("'eta' must be >= 1")
    }
    Data$eta <- eta

    if (!is.null(n.prior) && n.prior > 0) {
        banocc::cat_v("Begin getting prior plots.\n", verbose,
                      num_level=num_level + 1)

        prior.plots.normal <-
            banocc::get_LKJ_prior_plots(n.prior, eta=Data$eta,
                                        alpha_vec=Data$alpha,
                                        beta_vec=Data$beta, nu=Data$nu,
                                        Lambda=Data$Lambda)
        prior.plots.lognormal <-
            banocc::get_LKJ_prior_plots(n.prior, eta=Data$eta,
                                        alpha_vec=Data$alpha,
                                        beta_vec=Data$beta, nu=Data$nu,
                                        Lambda=Data$Lambda,
                                        transform="normal:lognormal")
        banocc::cat_v("End getting prior plots.\n", verbose,
                      num_level=num_level + 1)
        if (!is.null(workflow.name)){
            banocc::cat_v("Outputting prior plots...", verbose,
                          num_level=num_level + 1)

            pdf(paste0(workflow.name, "_normal_priors.pdf"))
            print(prior.plots.normal$mu.plot)
            print(prior.plots.normal$R.plot)
            print(prior.plots.normal$sigma.plot)
            dev.off()

            pdf(paste0(workflow.name, "_lognormal_priors.pdf"))
            print(prior.plots.lognormal$mu.plot)
            print(prior.plots.lognormal$R.plot)
            print(prior.plots.lognormal$sigma.plot)
            dev.off()
            banocc::cat_v("Done.\n")
        }
    }

    banocc::cat_v("Begin fitting the model\n", verbose, num_level=num_level+1)
    Fit.all <- banocc::mycapture(rstan::sampling(bayes_model, data=Data,
                                                 chains=chains, iter=iter,
                                                 warmup=warmup, thin=thin,
                                                 init=init, cores=cores))
    Fit <- Fit.all$output
    banocc::cat_v("End fitting the model\n", verbose, num_level=num_level+1)
    banocc::show_mcmc_time(Fit.all)
    cat("\n")

    if (!is.null(workflow.name)){
        banocc::cat_v("Outputting traceplots and act plots...", verbose,
                      num_level=num_level + 1)
        jpeg(paste0(workflow.name, "_traceplot.jpeg"))
        rstan::traceplot(Fit)
        dev.off()

        pdf(paste0(workflow.name, "_acf.pdf"))
        banocc::plot_stan_acf(Fit)
        dev.off()
        banocc::cat_v("Done.\n", verbose)
    }

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

    Estimates <-
        banocc::get_posterior_estimates(posterior_samples=post.samples.list,
                                        estimate_method="median",
                                        parameter.names="ln_Rho")
    dimnames(Estimates$ln_Rho) <- list(colnames(Data$C), colnames(Data$C))
    Estimates <- Estimates$ln_Rho

    
    workflow_return <- list(Data=Data, Fit=Fit, Fit.print=Fit.all$print.output,
                            CI.hpd=CI, Estimates.median=Estimates)

    if (!is.null(workflow.name)){
        workflow_return$prior.plots.normal    <- prior.plots.normal
        workflow_return$prior.plots.lognormal <- prior.plots.lognormal
    }
    banocc::cat_v("End stan_workflow\n", verbose, num_level=num_level)

    return(workflow_return)
}

#' @rdname stan_workflow
stan_workflow_sim <- function(bayes_model, C, nu, Lambda, alpha, beta,
                              eta, Rho.true,
                              cores = getOption("mc.cores", 1L),
                              chains = 4, iter = 2000, warmup = floor(iter/2),
                              thin = 1, init = 'random',
                              sd_mean=NULL, sd_var=NULL, n.prior=NULL,
                              workflow.name=NULL, verbose=FALSE, num_level=0){
    banocc::cat_v("Begin stan_workflow_sim\n", verbose, num_level=num_level)
    stan_workflow_return <-
        banocc::stan_workflow(bayes_model=bayes_model, n.prior=n.prior, C=C,
                              nu=nu, Lambda=Lambda, alpha=alpha, beta=beta,
                              eta=eta, workflow.name=workflow.name, cores=cores,
                              chains=chains, iter=iter, warmup=warmup, thin=thin,
                              init=init, sd_mean=sd_mean, sd_var=sd_var,
                              verbose=verbose, num_level=num_level + 1)

    if (!is.null(workflow.name)) {
        pdf(paste0(workflow.name, "_inference.pdf"))
        banocc::generate_inference_figure_MCMC(stan_workflow_return$Fit,
                                               Rho.true, ylim=c(-1, 1))
        dev.off()
    }

    banocc::cat_v("End stan_workflow_sim\n", verbose, num_level=num_level)
    return(stan_workflow_return)
}

check_nu <- function(nu, p, verbose=FALSE, num_level=0){
    banocc::cat_v("Begin check_nu\n", verbose, num_level=num_level)
    nu <- as.numeric(nu)
    nu <- banocc::check_vector("nu", nu, p, verbose, num_level=num_level+1)
    banocc::cat_v("End check_nu\n", verbose, num_level=num_level)
    return(nu)
}

check_alpha_beta <- function(alpha, beta, p, verbose=FALSE, num_level=0){
    banocc::cat_v("Begin check_alpha_beta\n", verbose, num_level=num_level)
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
    banocc::cat_v("End check_alpha_beta\n", verbose, num_level=num_level)
    return(list(alpha=alpha, beta=beta))
}

check_sd_mean_var <- function(sd_mean, sd_var, p, verbose=FALSE, num_level=0){
    banocc::cat_v("Begin check_sd_mean_var\n", verbose, num_level=num_level)
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
    banocc::cat_v("End check_sd_mean_var\n", verbose, num_level=num_level)
    return(list(sd_mean=sd_mean, sd_var=sd_var))
}

check_vector <- function(parm.name, parm, p, verbose=TRUE, num_level=0){
    banocc::cat_v("Begin check_vector...", verbose, num_level=num_level)
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
    banocc::cat_v("Done.\n", verbose)
    return(parm)
}

check_Lambda <- function(Lambda, p, verbose=FALSE, num_level=0){
    banocc::cat_v("Begin check_Lambda...", verbose, num_level=num_level)
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
            warn("'Lambda' is being recycled")
            num_rep <- ceiling(p / length(Lambda))
            Lambda <- diag(rep(Lambda, num_rep)[1:p])
        } else if (length(Lambda) == 1){
            Lambda <- diag(rep(Lambda, p))
        } else if (length(Lambda) > p){
            warning("'Lambda' has length > p; only first p elements will be ",
                    "used.")
            Lambda <- Lambda[1:p]
        }
    } else {
        stop("'Lambda' must be either a matrix or a vector")
    }
    banocc::cat_v("Done.\n", verbose)
    return(Lambda)
}
