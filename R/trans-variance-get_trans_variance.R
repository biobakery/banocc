#' Get the variance for a transformation of the parameters
#'
#' @param hessian.mat The Hessian matrix for the original parameters
#' @inheritParams get_gradient

get_trans_variance <-
function(hessian.mat, fn, param.u, Fit, epsilon=10e-6){
    ### param.u are the parameters at which the hessian.mat is calculated
    u.var <- solve(-hessian.mat)
    deriv <- banocc::get_gradient(fn, param.u, Fit, epsilon)

    trans.var <- deriv %*% u.var %*% t(deriv)
    return(trans.var)
}
