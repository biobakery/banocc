# A function to get pseudo p-values based on posterior credible intervals
#
# @inheritParams get_posterior_quantiles
# @param null.value The null value to which the credible interval is being
#   compared
# @param step.size The step size by which to increment the interval to
#   calculate the pseudo p-value.
# @return A list of pseudo p-values for each parameter in
#   \code{parameter.names}

get_pseudo_p_values <- function(posterior_samples, parameter.names="ln_Rho",
                                null.value=0, step.size=0.01, verbose=FALSE){
    if (verbose) cat("Starting pseudo p-value function.\n")
    alpha <- 1 - step.size
    post.quant <-
        get_posterior_quantiles(posterior_samples,
                                        probs=c(alpha/2, 1-alpha/2), list=TRUE,
                                        parameter.names=parameter.names)
    p.values <- lapply(post.quant, function(elt){
        if (is.vector(elt[[1]])) vector(mode="numeric",
                                        length=length(elt[[1]]))
        if (is.matrix(elt[[1]])) matrix(0, ncol=ncol(elt[[1]]),
                                        nrow=nrow(elt[[1]]))
    })

    cover.null <- lapply(post.quant, function(elt){
        (elt[[1]] <= null.value) * (elt[[2]] >= null.value)
    })

    p.values <- mapply(function(p.vals, stop.cond){
        p.vals[intersect(which(stop.cond>0), which(!p.vals))] <- alpha
        return(p.vals)
    }, p.values, cover.null, SIMPLIFY=FALSE)

    ## Stop when all CIs fail to cover the null
    stop.condition <- all(unlist(lapply(cover.null, all)))

    if(verbose){
        cat(paste0("alpha = ", alpha, "\n"))
        cat(paste0("num_covered = ", toString(unlist(lapply(cover.null, sum))),
                   "\n"))
    }

    alpha <- alpha - step.size
    while(!stop.condition && alpha > step.size/2){
        if (verbose) cat(paste0("alpha = ", alpha, "\n"))
        post.quant <-
            get_posterior_quantiles(posterior_samples,
                                            probs=c(alpha/2, 1-alpha/2),
                                            list=TRUE,
                                            parameter.names=parameter.names)
        cover.null <- lapply(post.quant, function(elt){
            (elt[[1]] <= null.value) * (elt[[2]] >= null.value)
        })

        p.values <- mapply(function(p.vals, stop.cond){
            p.vals[intersect(which(stop.cond>0), which(!p.vals))] <- alpha
            return(p.vals)
        }, p.values, cover.null, SIMPLIFY=FALSE)

        stop.condition <- all(unlist(lapply(cover.null, all)))
        if (verbose) cat(paste0("num_covered = ",
                                toString(lapply(cover.null, sum)), "\n"))
        alpha <- alpha - step.size
    }
    if (verbose) cat("Done getting pseudo p-values.\n")
    return(p.values)
}

        
