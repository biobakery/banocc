# Get posterior estimates and quantiles
#
# @param posterior_samples A list of the posterior samples for each parameter.
#   Each element should be an array of the samples across all chains for a
#   parameter. The first dimension of the array should be the sample id; the
#   others, the parameter index.
# @param estimate_method Either \code{"mean"} or \code{"median"}, indicating
#   which estimate to use.
# @param parameter.names A list of parameters for which estimates or
#   quantiles will be returned
# @param list A boolean value indicating whether to return the quantiles as
#   an array or list.
# @return A named list where each element is the estimate or quantiles for a
#   parameter in \code{parameter.names} The quantiles for a parameter are
#   either an array with the first dimension being the quantile and the
#   others the parameter indices; or a list where each element is the
#   quantile of a parameter.
# @name get_posterior_quantities
# 

# @rdname get_posterior_quantities
get_posterior_estimates <-
function(posterior_samples,
         estimate_method="mean",
         parameter.names=c("m", "S")){
    names(parameter.names) <- parameter.names
    
    if(estimate_method=="mean"){
        posterior_estimates <- lapply(parameter.names, ps_function,
                                      posterior_samples=posterior_samples,
                                      func=mean)
    } else if (estimate_method=="median"){
        posterior_estimates <-
            get_posterior_quantiles(posterior_samples, probs=0.5,
                                    parameter.names=parameter.names)
    } else {
        stop('Invalid estimate_method; must be "mean" or "median"')
    }
    return(posterior_estimates)
}

# @rdname get_posterior_quantities
get_posterior_quantiles <-
function(posterior_samples, probs, list=FALSE,
         parameter.names=c("m", "S")
         ){
    names(parameter.names) <- parameter.names
    posterior_quantiles <- lapply(parameter.names, ps_function,
                                  posterior_samples=posterior_samples,
                                  func=quantile, probs=probs)
    posterior_quantiles <- lapply(posterior_quantiles, function(pq){
        if (length(probs) == 1 && is.vector(pq)){
            unname(pq)
        } else {
            pq
        }
    })

    if(!list){
        return(posterior_quantiles)
    } else if (length(probs) > 1){
        posterior_quantiles.list <- make_list(
            array_list=posterior_quantiles, parameter.names=parameter.names,
            posterior_samples=posterior_samples, elt_length=length(probs))
        return(posterior_quantiles.list)
    } else {
        posterior_quantiles.list <-
            lapply(parameter.names,
                   function(name) posterior_quantiles[[name]])
    }
}
