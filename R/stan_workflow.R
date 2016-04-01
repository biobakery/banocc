#' Runs a stan workflow to fit the model and generate figures on SIMULATED DATA
#'
#' @inheritParams run_banocc
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
                          thin = 1, init = 'random', control=NULL,
                          sd_mean=NULL, sd_var=NULL, n.prior=NULL,
                          workflow.name=NULL, verbose=FALSE, num_level=0,
                          conf_alpha=0.05, get_min_width=FALSE){
    banocc::cat_v("Begin stan_workflow\n", verbose, num_level=num_level)
    workflow_return <- banocc::run_banocc(bayes_model=bayes_model, C=C, nu=nu,
                                          Lambda=Lambda, alpha=alpha, beta=beta,
                                          eta=eta, cores=cores, chains=chains,
                                          iter=iter, warmup=warmup, thin=thin,
                                          init=init, control=control,
                                          sd_mean=sd_mean, sd_var=sd_var,
                                          conf_alpha=conf_alpha,
                                          get_min_width=get_min_width,
                                          verbose=verbose,
                                          num_level=num_level+1)

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

    if (!is.null(workflow.name)){
        banocc::cat_v("Begin outputting trace and acf plots.\n", verbose,
                      num_level=num_level + 1)

        traceplot.name <- paste0(workflow.name, "_traceplot_%04d.jpeg")
        jpeg(traceplot.name)
        rstan::traceplot(workflow_return$Fit,
                         pars=attr(workflow_return$Fit, "model_pars"))
        dev.off()

        pdf(paste0(workflow.name, "_acf.pdf"))
        banocc::plot_stan_acf_set(workflow_return$Fit, NumChains=chains,
                                  p=Data$P, verbose=verbose,
                                  num_level=num_level+2)
        dev.off()
        banocc::cat_v("End outputting trace and acf plots.\n", verbose,
                      num_level=num_level+1)
    }
    
    if (!is.null(workflow.name) && !is.null(n.prior) && n.prior > 0){
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
                              thin = 1, init = 'random', control=NULL,
                              sd_mean=NULL, sd_var=NULL, n.prior=NULL,
                              workflow.name=NULL, verbose=FALSE, num_level=0){
    banocc::cat_v("Begin stan_workflow_sim\n", verbose, num_level=num_level)
    stan_workflow_return <-
        banocc::stan_workflow(bayes_model=bayes_model, n.prior=n.prior, C=C,
                              nu=nu, Lambda=Lambda, alpha=alpha, beta=beta,
                              eta=eta, workflow.name=workflow.name, cores=cores,
                              chains=chains, iter=iter, warmup=warmup, thin=thin,
                              init=init, control=control, sd_mean=sd_mean,
                              sd_var=sd_var, verbose=verbose,
                              num_level=num_level + 1)

    if (!is.null(workflow.name)) {
        pdf(paste0(workflow.name, "_inference.pdf"))
        banocc::generate_inference_figure_MCMC(stan_workflow_return$Fit,
                                               Rho.true, ylim=c(-1, 1),
                                               verbose=verbose,
                                               num_level=num_level+1)
        dev.off()
    }

    banocc::cat_v("End stan_workflow_sim\n", verbose, num_level=num_level)
    return(stan_workflow_return)
}
