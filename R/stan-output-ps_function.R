# Apply a function to an array of posterior samples so that correct
#  dimension for matrix, vector and scalar parameters is returned
#
# @param func is the function to use
# @param name The parameter name
# @param posterior_samples The posterior samples as a list. \code{name}
#   must be one of the names of \code{posterior_samples}
#
ps_function <- function(name, posterior_samples, func, ...){
    is.mat <- length(dim(posterior_samples[[name]])) == 3
    is.vec <- length(dim(posterior_samples[[name]])) == 2
    if(is.mat){
        apply(posterior_samples[[name]], c(2, 3), FUN=func, ...)
    } else if (is.vec) {
        apply(posterior_samples[[name]], 2, FUN=func, ...)
    } else {
        func(posterior_samples[[name]], ...)
    }
}
