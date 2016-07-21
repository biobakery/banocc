# Get credible intervals for the output according to different paradigms
#
# @param type One of \code{"marginal.centered"} or \code{"marginal.hpd"},
#   indicating whether centeral or highest posterior density intervals are
#   desired
# @param conf The width of the credible interval (conf * 100%).
# @inheritParams get_posterior_quantiles
# @inheritParams cat_v
# @return Returns a list of the intervals for each parameter in
#   \code{parameter.names}.  If \code{list=TRUE}, then each parameter is a
#   list, with elements \code{lower} and \code{upper} being the lower and
#   upper bounds, respectively. If \code{list=FALSE}, then each parameter is
#   an array with the first dimension being the bounds and the remaining
#   dimensions being the parameter index or indices.
#
#' @importFrom coda as.mcmc
#' @importFrom coda HPDinterval
get_credible_intervals <- function(posterior_samples, list=FALSE,
                                   parameter.names=c("m", "S"),
                                   conf=0.95,
                                   type="marginal.hpd",
                                   verbose=FALSE, num_level=0){
    cat_v("Begin get_credible_intervals\n", verbose,
                  num_level=num_level)
    if (type == "marginal.centered"){
        alpha <- 1 - conf
        credible.intervals <-
            get_posterior_quantiles(posterior_samples, list=list,
                                            probs=c((alpha/2), 1 - (alpha/2)),
                                            parameter.names=parameter.names)
        if (!list){
            credible.intervals <- lapply(credible.intervals, function(param.ci){
                if (!is.null(dim(param.ci))){
                    dimnames(param.ci)[[1]] <- c("lower", "upper")
                } else {
                    names(param.ci) <- c("lower", "upper")
                }
                return(param.ci)
            })
        } else {
            credible.intervals <- lapply(credible.intervals, function(param.ci){
                names(param.ci) <- c("lower", "upper")
                return(param.ci)
            })
        }
    } else if (type == "marginal.hpd"){
        names(parameter.names) <- parameter.names
        credible.intervals <- lapply(parameter.names, function(name){
            is.mat <- length(dim(posterior_samples[[name]])) == 3
            is.vec <- length(dim(posterior_samples[[name]])) == 2
            if(is.mat){
                apply(posterior_samples[[name]], c(2, 3),
                      function(param){
                          coda::HPDinterval(coda::as.mcmc(param), prob=conf)
                      })
            } else if (is.vec){
                apply(posterior_samples[[name]], 2,
                      function(param){
                          coda::HPDinterval(coda::as.mcmc(param), prob=conf)
                      })
            } else {
                as.vector(coda::HPDinterval(
                    coda::as.mcmc(as.vector(posterior_samples[[name]])),
                    prob=conf))
            }
        })
        credible.intervals <- lapply(credible.intervals, function(param.ci){
            if (!is.null(dim(param.ci))){
                dimnames(param.ci)[[1]] <- c("lower", "upper")
            } else {
                names(param.ci) <- c("lower", "upper")
            }
            return(param.ci)
        })
        if (list){
            credible.intervals.list <-
                lapply(parameter.names,
                       function(name){
                           new.vec <-vector("list", length=2)
                           names(new.vec) <- c("lower", "upper")
                           return(new.vec)
                       })
            for (i in 1:2){
                for(name in parameter.names){
                    is.mat <- length(dim(posterior_samples[[name]])) == 3
                    is.vec <- length(dim(posterior_samples[[name]])) == 2
                    if(is.mat){
                        credible.intervals.list[[name]][[i]] <-
                            credible.intervals[[name]][i, , ]
                    } else if (is.vec){
                        credible.intervals.list[[name]][[i]] <-
                            credible.intervals[[name]][i, ]
                    } else {
                        credible.intervals.list[[name]][[i]] <-
                            credible.intervals[[name]][i]
                    }
                }
            }
            credible.intervals <- credible.intervals.list
        }
    } else {
        stop("'type' must be one of \"marginal.centered\" or \"marginal.hpd\"")
    }
    cat_v("End get_credible_intervals\n", verbose, num_level=num_level)
    return(credible.intervals)
}
