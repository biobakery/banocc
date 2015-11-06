#' Convert mean and variance-covariance samples to lognormal parameters or
#'   normal parameters
#'
#' @param Sigma.sample A list of variance-covariance matrices
#' @param mu.sample A list of mean vectors
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @return A named list of converted parameters: \code{_mean} is the converted
#'   mean; \code{_Var} the converted variance; \code{_corr} the converted
#'   correlation matrix; and \code{_sd} the converted standard deviation
#'   vectors.
#' @examples
#' convert_prior_to_lognormal(diag(4), rep(0, 4))
#' convert_prior_to_normal(diag(4), rep(1, 4))
#' @name convert_prior_samples

#' @rdname convert_prior_samples
convert_prior_to_lognormal <-
function(Sigma.sample, mu.sample, verbose=FALSE, num_level=0){
  banocc::cat_v("Begin convert_prior_to_lognormal\n", verbose,
                num_level=num_level)
    
  banocc::cat_v("Generating lognormal parameters...", verbose, num_level+1)
  ln.prior.sample <- mapply(banocc::get_ln_params, mu.sample, Sigma.sample, 
                             SIMPLIFY=FALSE)
  param.names <- unlist(unique(lapply(ln.prior.sample, names)))
  names(param.names) <- param.names
  ln.prior.sample <- lapply(param.names,
                            function(name) lapply(ln.prior.sample, '[[', name))
  banocc::cat_v("Done.\n", verbose)
  
  banocc::cat_v("Getting lognormal correlation and standard deviation...",
                verbose, num_level=num_level+1)
  ln.prior.sample$ln_corr <- lapply(ln.prior.sample$ln_Var, cov2cor)
  ln.prior.sample$ln_sd   <- lapply(ln.prior.sample$ln_Var, 
                                    function(S) sqrt(diag(S)))
  banocc::cat_v("Done.\n", verbose)

  banocc::cat_v("End convert_prior_to_lognormal\n", verbose,
                num_level=num_level)
  return(ln.prior.sample)
}

#' @rdname convert_prior_samples
get_norm.prior.sample <-
function(Sigma.sample, mu.sample, verbose=FALSE, num_level=0){

  banocc::cat_v("Begin get_norm.prior.sample\n", verbose, num_level=num_level)
  
  banocc::cat_v("Generating normal parameters...", verbose,
                num_level=num_level+1)
  norm.prior.sample <- mapply(banocc::get_n_params, mu.sample, Sigma.sample,
                            SIMPLIFY=FALSE)
  param.names <- unlist(unique(lapply(norm.prior.sample, names)))
  names(param.names) <- param.names
  norm.prior.sample <- lapply(param.names, function(name)
                              lapply(norm.prior.sample, '[[', name))
  banocc::cat_v("Done.\n", verbose)
  
  banocc::cat_v("Getting normal correlation...", verbose, num_level=num_level+1)
  norm.prior.sample$n_corr <- lapply(norm.prior.sample$n_Var, cov2cor)
  norm.prior.sample$n_sd   <- lapply(norm.prior.sample$n_Var, 
                                     function(S) sqrt(diag(S)))
  banocc::cat_v("Done.\n", verbose)

  banocc::cat_v("End get_norm.prior.sample\n", verbose, num_level=num_level)
  
  return(norm.prior.sample)
}
