#' Generate plots of the prior distributions for all parameters for a few
#'   special prior cases
#'
#' @inheritParams MCMCpack::riwish
#' @inheritParams plot_mu_prior
#' @inheritParams sample_mu_prior
#' @inheritParams get_vec_text
#' @inheritParams cat_v
#' @param n The number of samples to use in generating the plots
#' @inheritParams rlkj
#' @param gamma_mu Boolean: is mu a gamma distribution with parameters
#'   \code{nu} and \code{diag(Lambda)}?
#' @param transform one of \code{"none"}, \code{"normal:lognormal"}, or
#'   \code{"lognormal:normal"}, indicating a transformation of the parameters
#'   to be made. This allows priors corresponding to the parameter
#'   transformation to be plotted.
#' @param plot.Sigma Boolean: Should the prior distribution of
#'   \code{Sigma = diag(sigma) \%*\% R \%*\% diag(sigma)} be plotted?
#' @return A named list: \code{.plot} are the \code{ggplot} objects, which can
#'   be printed; \code{.data} are the data objects.    
#' @name get_prior_plot_sets
#'

#' @rdname get_prior_plot_sets
get_IW_prior_plots <-
function(n, v, S, nu, Lambda, gamma_mu = FALSE, lognorm_mu = FALSE,
         transform="none", add.stats=FALSE, digits=2, center.stat="mean",
         verbose=FALSE){
  if (!(transform %in% c('none', 'normal:lognormal', 'lognormal:normal'))){
    stop("transform must be one of: 'none', 'normal:lognormal', 'lognormal:normal'")
  }
  IW.sample <- banocc::sample_IW_prior(n, S, v, verbose)
  if (gamma_mu){
      mu.alpha <- nu^2/diag(Lambda)
      mu.beta  <- nu/diag(Lambda)
      mu.sample <- banocc::sample_sigma_prior(n, mu.alpha, mu.beta, verbose)
  } else {
      mu.sample <- banocc::sample_mu_prior(n, nu, Lambda, lognorm_mu, verbose)
  }
  
  if(transform=="normal:lognormal") {
    ln.prior.sample <- banocc::get_ln.prior.sample(IW.sample, mu.sample, verbose)
    Sigma.plot.stuff <- banocc::plot_Sigma_prior(ln.prior.sample$ln_Var,
                                                 add.stats=add.stats,
                                                 center.stat=center.stat,
                                                 digits=digits, verbose=verbose)
    mu.plot.stuff    <- banocc::plot_mu_prior(ln.prior.sample$ln_mean,
                                              add.stats=add.stats,
                                              center.stat=center.stat,
                                              digits=digits, verbose=verbose)
  } else if(transform=="lognormal:normal"){
    norm.prior.sample <- banocc::get_norm.prior.sample(IW.sample, mu.sample,
                                                       verbose)
    Sigma.plot.stuff  <- banocc::plot_Sigma_prior(norm.prior.sample$n_Var,
                                                  add.stats=add.stats,
                                                  center.stat=center.stat,
                                                  digits=digits, verbose=verbose)
    mu.plot.stuff     <- banocc::plot_mu_prior(norm.prior.sample$n_mean,
                                               add.stats=add.stats,
                                               center.stat=center.stat,
                                               digits=digits, verbose=verbose)
  } else {
    Sigma.plot.stuff <- banocc::plot_Sigma_prior(IW.sample, add.stats=add.stats,
                                                 center.stat=center.stat,
                                                 digits=digits, verbose=verbose)
    mu.plot.stuff    <- banocc::plot_mu_prior(mu.sample, add.stats=add.stats,
                                              center.stat=center.stat,
                                              digits=digits, verbose=verbose)
  }
  
  return(list(mu.plot=mu.plot.stuff$plot,
              mu.data=mu.plot.stuff$data,
              Sigma.plot=Sigma.plot.stuff$plot,
              Sigma.data=Sigma.plot.stuff$data))
}

#' @rdname get_prior_plot_sets
get_LKJ_prior_plots <- function(n, eta, alpha_vec,
                                beta_vec, nu, Lambda, gamma_mu = FALSE,
                                lognorm_mu = FALSE, transform="none",
                                add.stats=FALSE, digits=2, center.stat="mean",
                                plot.Sigma=FALSE, verbose=FALSE, num_level=FALSE){

  banocc::cat_v("Begin get_LKJ_prior_plots.\n", verbose, num_level=num_level)
  
  if (!(transform %in% c('none', 'normal:lognormal', 'lognormal:normal'))){
    stop("transform must be one of: 'none', 'normal:lognormal', 'lognormal:normal'")
  }
  
  p <- length(alpha_vec)
  
  LKJ.sample <- banocc::sample_lkj_prior(n=n, d=p, eta=eta, verbose=verbose,
                                         num_level=num_level+1)
  sigma.sample <- banocc::sample_sigma_prior(n=n, alpha_vec=alpha_vec,
                                             beta_vec=beta_vec, verbose=verbose,
                                             num_level=num_level+1)
  if (gamma_mu){
      mu.alpha <- nu^2/diag(Lambda)
      mu.beta  <- nu/diag(Lambda)
      mu.sample <- banocc::sample_sigma_prior(n=n, alpha_vec=mu.alpha,
                                              beta_vec=mu.beta, verbose=verbose,
                                              num_level=num_level + 1)
  } else {
      mu.sample <- banocc::sample_mu_prior(n=n, nu=nu, Lambda=Lambda,
                                           lognorm_mu=lognorm_mu,
                                           verbose=verbose,
                                           num_level=num_level+1)
  }
  warn_exponent(lapply(sigma.sample, "^", 2), verbose)

  banocc::cat_v("Getting Sigma sample...", verbose, num_level + 1)
  Sigma.sample <- mapply(function(R, sigma){
    diag(sigma) %*% R %*% diag(sigma)      
  }, LKJ.sample, sigma.sample, SIMPLIFY=FALSE)
  banocc::cat_v("Done.\n", verbose)
  
  if (transform=="normal:lognormal"){
    warn_exponent(Sigma.sample, verbose)
    warn_exponent(mu.sample, verbose)
    ln.prior.sample <- banocc::convert_prior_to_lognormal(Sigma.sample,
                                                          mu.sample,
                                                          verbose,
                                                          num_level=num_level+1)
    
    banocc::cat_v("Getting summary statistics...", verbose,
                  num_level=num_level+1)
    ln_mean.stats  <- banocc::get_prior_stats_to_print(ln.prior.sample$ln_mean)
    ln_corr.stats  <- banocc::get_prior_stats_to_print(ln.prior.sample$ln_corr)
    ln_sigma.stats <- banocc::get_prior_stats_to_print(ln.prior.sample$ln_sd)
    banocc::cat_v("Done.\n", verbose)
    
    to_return <- list(mu.stats    = ln_mean.stats,
                      R.stats     = ln_corr.stats,
                      sigma.stats = ln_sigma.stats)
    ## if (verbose){
    ##     cat("LN_MEAN median (IQR)\n")
    ##     cat(paste0(ln_mean.stats$median, "(", ln_mean.stats$IQR, ")\n"))

    ##     cat("LN_CORR median (IQR)\n")
    ##     cat(paste0(ln_corr.stats$median, "(", ln_corr.stats$IQR, ")\n"))

    ##     cat("LN_SIGMA median (IQR)\n")
    ##     cat(paste0(ln_sigma.stats$median, "(", ln_sigma.stats$IQR, ")\n"))
    ## }
    if (plot.Sigma){
      Sigma.plot.stuff <- banocc::plot_Sigma_prior(ln.prior.sample$ln_Var,
                                                   add.stats=add.stats,
                                                   center.stat=center.stat,
                                                   digits=digits,
                                                   verbose=verbose,
                                                   num_level=num_level+1)
    }
    if (verbose) cat("R.plot\n")
    R.plot.stuff <- banocc::plot_R_prior(ln.prior.sample$ln_corr,
                                         add.stats=add.stats,
                                         center.stat=center.stat, digits=digits,
                                         verbose=verbose,
                                         num_level=num_level+1)
    if (verbose) cat("sigma plot\n")
    sigma.plot.stuff <- banocc::plot_sigma_prior(ln.prior.sample$ln_sd,
                                                 add.stats=add.stats,
                                                 center.stat=center.stat,
                                                 digits=digits, verbose=verbose,
                                                 num_level=num_level+1)
    if (verbose) cat("mu plot\n")
    mu.plot.stuff    <- banocc::plot_mu_prior(ln.prior.sample$ln_mean,
                                              add.stats=add.stats,
                                              center.stat=center.stat,
                                              digits=digits, verbose=verbose,
                                              num_level=num_level+1)
  } else if(transform=="lognormal:normal"){
      to_return <- list()
      ## if (verbose){
      ##     print(mapply(function(Sigma, mu) Sigma * 1/(mu %*% t(mu)),
      ##                  Sigma.sample, mu.sample, SIMPLIFY = FALSE))
      ##     print(Sigma.sample)
      ##     print(lapply(mu.sample, function(mu) 1/(mu %*% t(mu))))
      ## }
      norm.prior.sample <-
          banocc::get_norm.prior.sample(Sigma.sample=Sigma.sample,
                                        mu.sample=mu.sample, verbose=verbose,
                                        num_level=num_level+1)

      if (plot.Sigma){
          Sigma.plot.stuff <- banocc::plot_Sigma_prior(norm.prior.sample$n_Var,
                                                       add.stats=add.stats,
                                                       center.stat=center.stat,
                                                       digits=digits,
                                                       verbose=verbose,
                                                       num_level=num_level+1)
      }
      R.plot.stuff <- banocc::plot_R_prior(norm.prior.sample$n_corr,
                                           add.stats=add.stats,
                                           center.stat=center.stat,
                                           digits=digits, verbose=verbose,
                                           num_level=num_level+1)
      sigma.plot.stuff <- banocc::plot_sigma_prior(norm.prior.sample$n_sd,
                                                   add.stats=add.stats,
                                                   center.stat=center.stat,
                                                   digits=digits,
                                                   verbose=verbose,
                                                   num_level=num_level+1)
      mu.plot.stuff    <- plot_mu_prior(norm.prior.sample$n_mean,
                                        add.stats=add.stats,
                                        center.stat=center.stat, digits=digits,
                                        verbose=verbose,
                                        num_level=num_level+1)
  } else {
    to_return <- list()

    if (plot.Sigma){
      Sigma.plot.stuff <- banocc::plot_Sigma_prior(Sigma.sample,
                                                   add.stats=add.stats,
                                                   center.stat=center.stat,
                                                   digits=digits,
                                                   verbose=verbose,
                                                   num_level=num_level+1)
    }
    R.plot.stuff <- banocc::plot_R_prior(LKJ.sample, add.stats=add.stats,
                                         center.stat=center.stat, digits=digits,
                                         verbose=verbose,
                                         num_level=num_level+1)
    sigma.plot.stuff <- banocc::plot_sigma_prior(sigma.sample,
                                                 add.stats=add.stats,
                                                 center.stat=center.stat,
                                                 digits=digits, verbose=verbose,
                                                 num_level=num_level+1)
    mu.plot.stuff    <- banocc::plot_mu_prior(mu.sample, add.stats=add.stats,
                                              digits=digits, verbose=verbose,
                                              center.stat=center.stat,
                                              num_level=num_level+1)
  }

  to_return$R.plot <- R.plot.stuff$plot
  to_return$R.data <- R.plot.stuff$data
  to_return$sigma.plot <- sigma.plot.stuff$plot
  to_return$sigma.data <- sigma.plot.stuff$data
  to_return$mu.plot <- mu.plot.stuff$plot
  to_return$mu.data <- mu.plot.stuff$data

  if (plot.Sigma){
      to_return$Sigma.plot <- Sigma.plot.stuff$plot
      to_return$Sigma.data <- Sigma.plot.stuff$data
  }
  banocc::cat_v("End get_LKJ_prior_plots.\n", verbose, num_level=num_level)
  return(to_return)
}
