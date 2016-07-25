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
            credible.intervals <- lapply(credible.intervals, function(p.ci){
                if (!is.null(dim(p.ci))){
                    dimnames(p.ci)[[1]] <- c("lower", "upper")
                } else {
                    names(p.ci) <- c("lower", "upper")
                }
                return(p.ci)
            })
        } else {
            credible.intervals <- lapply(credible.intervals, function(p.ci){
                names(p.ci) <- c("lower", "upper")
                return(p.ci)
            })
        }
    } else if (type == "marginal.hpd"){
        names(parameter.names) <- parameter.names
        credible.intervals <- lapply(
            parameter.names, ps_function,
            posterior_samples=posterior_samples,
            func=function(param){
                coda::HPDinterval(coda::as.mcmc(as.vector(param)), prob=conf)
            })
        credible.intervals <- lapply(credible.intervals, function(p.ci){
            if (dim(p.ci)[[1]]==2){
                dimnames(p.ci)[[1]] <- c("lower", "upper")
            } else {
                p.ci <- p.ci[1, ]
                names(p.ci) <- c("lower", "upper")
            }
            return(p.ci)
        })
        if (list){
            credible.intervals.list <- make_list(
                array_list=credible.intervals,
                parameter.names=parameter.names,
                posterior_samples=posterior_samples, elt_length=2,
                elt_names=c("lower", "upper")) 
            credible.intervals <- credible.intervals.list
        }
    } else {
        stop("'type' must be one of \"marginal.centered\" or \"marginal.hpd\"")
    }
    cat_v("End get_credible_intervals\n", verbose, num_level=num_level)
    return(credible.intervals)
}
