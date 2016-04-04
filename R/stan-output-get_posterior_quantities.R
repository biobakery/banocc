#' Get posterior estimates and quantiles
#'
#' @param posterior_samples A list of the posterior samples for each parameter.
#'   Each element should be an array of the samples across all chains for a
#'   parameter. The first dimension of the array should be the sample id; the
#'   others, the parameter index.
#' @param estimate_method Either \code{"mean"} or \code{"median"}, indicating
#'   which estimate to use.
#' @param parameter.names A list of parameters for which estimates or
#'   quantiles will be returned
#' @param list A boolean value indicating whether to return the quantiles as
#'   an array or list.
#' @return A named list where each element is the estimate or quantiles for a
#'   parameter in \code{parameter.names} The quantiles for a parameter are
#'   either an array with the first dimension being the quantile and the
#'   others the parameter indices; or a list where each element is the
#'   quantile of a parameter.
#' @name get_posterior_quantities
#' 

#' @rdname get_posterior_quantities
get_posterior_estimates <-
function(posterior_samples,
         estimate_method="mean",
         parameter.names=c("mu", "Sigma")){
    names(parameter.names) <- parameter.names
    
    if(estimate_method=="mean"){
        posterior_estimates <- lapply(parameter.names, function(name){
            is.mat <- length(dim(posterior_samples[[name]])) == 3
            if(is.mat){
                apply(posterior_samples[[name]], c(2, 3), mean)
            } else{
                apply(posterior_samples[[name]], 2, mean)
            }
        })
    } else if (estimate_method=="median"){
        posterior_estimates <-
            banocc::get_posterior_quantiles(posterior_samples, probs=0.5,
                                            parameter.names=parameter.names)
    } else {
        stop('Invalid estimate_method; must be "mean" or "median"')
    }
    return(posterior_estimates)
}

#' @rdname get_posterior_quantities
get_posterior_quantiles <-
function(posterior_samples, probs, list=FALSE,
         parameter.names=c("mu", "Sigma")
         ){
    names(parameter.names) <- parameter.names
    posterior_quantiles <- lapply(parameter.names, function(name){
        is.mat <- length(dim(posterior_samples[[name]])) == 3
        is.vec <- length(dim(posterior_samples[[name]])) == 2
        if(is.mat){
            apply(posterior_samples[[name]], c(2, 3), quantile, probs=probs)
        } else if (is.vec) {
            apply(posterior_samples[[name]], 2, quantile, probs=probs)
        } else {
            matrix(quantile(posterior_samples[[name]], probs=probs),
                   ncol=1)
        }
    })

    if(!list){
        return(posterior_quantiles)
    } else if (length(probs) > 1){
        posterior_quantiles.list <-
            lapply(parameter.names,
                   function(name) vector("list", length=length(probs))
                   )
        for (i in seq_along(probs)){
            for(name in parameter.names){
                is.mat <- length(dim(posterior_samples[[name]])) == 3
                if(is.mat){
                  posterior_quantiles.list[[name]][[i]] <-
                      posterior_quantiles[[name]][i, , ]
                } else {
                    posterior_quantiles.list[[name]][[i]] <-
                        posterior_quantiles[[name]][i, ]
                }
            }
        }
        return(posterior_quantiles.list)
    } else {
        posterior_quantiles.list <-
            lapply(parameter.names,
                   function(name) posterior_quantiles[[name]])
    }
}
