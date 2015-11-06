#' Sampling from various prior distributions
#'
#' @param n number of samples to draw from the distribution
#' @param verbose Boolean: print informative statements as
#' function runs?
#' @inheritParams MCMCpack::riwish
#' @inheritParams rlkj
#' @inheritParams cat_v
#' @param nu mean parameter for multivariate normal distribution
#' @param Lambda variance-covariance matrix for the multivariate normal
#'   distribution
#' @param lognorm_mu Boolean: exponentiate the mu samples after drawing?
#' @param alpha_vec,beta_vec vectors of alpha and beta parameters for a gamma
#'   distribution (see \code{\link{rgamma}}). The vectors must be paired
#'   (ie., the first element of \code{alpha_vec} and the first element of
#'   \code{beta_vec} means that the first sigma is distributed
#'   gamma(\code{alpha_vec[1]}, \code{beta_vec[1]}))
#' @name sample_prior
#' @importFrom MCMCpack riwish
#' @importFrom mvtnorm rmvnorm

#' @rdname sample_prior
sample_IW_prior <-
function(n, S, v, verbose=FALSE, num_level=0){
  banocc::cat_v("Begin sample_IW_prior...", verbose, num_level=num_level)
#  banocc::cat_v("Sampling IW prior...", verbose, num_level=num_level+1)
  Sigma.sample <- lapply(1:n, 
                         function(i){MCMCpack::riwish(S=S, v=v)}
                         )
  banocc::cat_v("Done\n", verbose)
  return(Sigma.sample)
}

#' @rdname sample_prior
sample_lkj_prior <-
function(n, d, eta, verbose=FALSE, num_level=0){
  banocc::cat_v("Begin sample_lkj_prior...", verbose, num_level=num_level)
#  banocc::cat_v("Sampling from LKJ prior...", verbose, num_level=num_level+1)
  lkj.sample <- lapply(1:n, function(i) banocc::rlkj(d, eta))
  banocc::cat_v("Done.\n", verbose)
  return(lkj.sample)
}

#' @rdname sample_prior
sample_mu_prior <-
function(n, nu, Lambda, lognorm_mu = FALSE, verbose=FALSE, num_level=0){
  banocc::cat_v("Begin sample_mu_prior...", verbose, num_level=num_level)
#  banocc::cat_v("Sampling from mu prior...", verbose, num_level)
  mu <- mvtnorm::rmvnorm(n,mean=nu,sigma=Lambda)
  if (lognorm_mu){
      mu <- exp(mu)
  }
  mu.sample <- as.list(as.data.frame(t(mu)))
  names(mu.sample) <- NULL
  banocc::cat_v("Done.\n", verbose)
  return(mu.sample)
}

#' @rdname sample_prior
sample_sigma_prior <-
function(n, alpha_vec, beta_vec, verbose=FALSE, num_level=0){
  if (length(alpha_vec) != length(beta_vec)){
    stop("Lengths of alpha_vec and beta_vec must equal.")
  }
  banocc::cat_v("Begin sample_sigma_prior...", verbose, num_level=num_level)
#  banocc::cat_v("Sampling variance...", verbose, num_level)
  sigma.sample <- lapply(1:n, function(i){
    mapply(function(a, b) rgamma(1, shape=a, rate=b), alpha_vec, beta_vec)
  })
  banocc::cat_v("Done.\n", verbose)
  return(sigma.sample)
}
