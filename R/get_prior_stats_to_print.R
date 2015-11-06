#' Get summary statistics for each parameter that can be printed.
#'
#' @param param.samples A list of parameter samples, which can be scalar or
#' vector or matrix.
#' 
#' @return A list of median and IQR for each parameter.
#' 
#' @examples
#' get_prior_stats_to_print(rnorm(10, 0, 1))
#' 

get_prior_stats_to_print <- function(param.samples){
    param.nomiss.id <- unlist(lapply(param.samples,
                                     function(p){
                                         if(any(is.na(p))){ FALSE } else {TRUE}
                                     }))
    param.nomiss <- param.samples[param.nomiss.id]
    n <- length(param.nomiss)

    ## param.avg <- Reduce('+', lapply(param.nomiss, '/', n))
    ## param.scaled.ss  <- Reduce('+', lapply(param.nomiss,
    ##                                        function(p){
    ##                                            (p / sqrt(n - 1))^2
    ##                                        }))
    ## param.var <- param.scaled.ss - (n / (n - 1)) * param.avg^2

    param.median         <- param.samples[[1]]
    param.upper.quartile <- param.samples[[1]]
    param.lower.quartile <- param.samples[[1]]
    for(k in seq_along(param.median)){
        param.vec.k <- unlist(lapply(param.nomiss, '[', k))
        param.median[k]         <- median(param.vec.k)
        param.upper.quartile[k] <- quantile(param.vec.k, probs=0.75)
        param.lower.quartile[k] <- quantile(param.vec.k, probs=0.25)
    }
    return(list(
        ## mean=param.avg, var = param.var,
        median=param.median, IQR = param.upper.quartile - param.lower.quartile
        ## upper.quartile=param.upper.quartile,
        ## lower.quartile=param.lower.quartile
        ))
}
