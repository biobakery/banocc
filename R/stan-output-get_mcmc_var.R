# Calculates the variance of the MCMC chains
# 
# @inheritParams get_posterior_quantities
# @inheritParams cat_v
# @return Returns a list of variance values for each parameter in
#   \code{parameter.names}
#

get_mcmc_var <- function(posterior_samples, parameter.names=c("m", "S"),
                         verbose=FALSE, num_level=0){
    cat_v("Begin get_mcmc_var\n", verbose, num_level=num_level)
    names(parameter.names) <- parameter.names
    mcmc_var_values <- lapply(parameter.names, function(name){
        is.mat <- length(dim(posterior_samples[[name]])) == 3
        is.vec <- length(dim(posterior_samples[[name]])) == 2
        if (is.mat){
            apply(posterior_samples[[name]], c(2, 3), var)
        } else if (is.vec){
            apply(posterior_samples[[name]], 2, var)
        } else {
            var(as.vector(posterior_samples[[name]]))
        }

    })
    cat_v("End get_mcmc_var\n", verbose, num_level=num_level)
    return(mcmc_var_values)
}
