#' Convert normal and lognormal mean/variance parameters
#'
#' @param mean the numeric mean of the lognormal or
#'   normal distribution; for the lognormal distribution, must be
#'   greater than zero.
#' @param Var the numeric variance-covariance matrix of the lognormal or normal
#'   distribution.
#' @return A named list: \code{_mean} is the converted mean and \code{_Var}
#'   is the converted variance-covariance matrix.
#' @examples
#' get_n_params(rep(1, 4), diag(4))
#' get_ln_params(rep(1, 4), diag(4))
#' get_n_params(get_ln_params(rep(0, 4), diag(4)))
#' 
#' @name convert_params

#' @rdname convert_params
get_n_params <-
function(mean,Var){
 n_Var <- log(Var * 1/(mean %*% t(mean)) + 1)
 n_mean <- log(mean) - 0.5 * diag(n_Var)
 return(list(n_mean=n_mean, n_Var=n_Var))
}

#' @rdname convert_params
get_ln_params <-
function(mean,Var){
  ln_mean <- exp(mean + 0.5 * diag(Var))
  ln_Var <- ln_mean %*% t(ln_mean) * (exp(Var) - 1)
  return(list(ln_mean=ln_mean, ln_Var=ln_Var))
}
