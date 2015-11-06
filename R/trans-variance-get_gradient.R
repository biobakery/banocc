#' Estimate the gradient of a multivariate function at a particular domain
#'   value
#'
#' @param fn The function whose gradient will be estimated.  It must take two
#'   parameters: the first, a vector of parameter values; the second, a
#'   \code{stanfit} object.
#' @param param.u The parameter values at which the gradient will be estimated
#' @param Fit A \code{stanfit} object (to be used as an argument to \code{fn}).
#' @param epsilon The step size used in estimating the derivative
#'
#' @examples
#' myfn <-
#' function(param.u, Fit){
#'     param.c <- constrain_pars(Fit, param.u)
#'
#'     fn.param <- c(param.c$ln_mu,
#'                   param.c$ln_Rho[lower.tri(param.c$ln_Rho)],
#'                   sqrt(diag(param.c$ln_Sigma)))
#'     names(fn.param) <- c(rep("ln_mu", length(param.c$ln_mu)),
#'                          paste0("ln_Rho.", which(lower.tri(param.c$ln_Rho))),
#'                          rep("ln_sigma", length(diag(param.c$ln_Sigma))))
#'     return(fn.param)
#' }

get_gradient <-
function(fn, param.u, Fit, epsilon=1e-6){
    ### Columns are the original parameters, rows
    ### are the elements of fn; entries are derivatives of fn element i
    ### with respect to original parameter j

    I <- diag(length(param.u))
    grad.est <- apply(I, 1, function(i){
        (fn(param.u + i*epsilon, Fit) - fn(param.u, Fit))/epsilon
    }
                      )
    return(grad.est)
}
