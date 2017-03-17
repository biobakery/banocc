# Generate a list of samples when they have been extracted from a
#   \code{stanfit} object with \code{permuted=FALSE}.
#
# @param samples A set of samples extracted from a \code{stanfit} object with
#   \code{permuted=FALSE}
# @inheritParams rstan::sampling
# @param concatenate.chains Boolean: should all the chains be concatenated
#   together, or left separately?
# @inheritParams cat_v
# @return A named list of each parameter from the \code{stanfit} object. The
#   elements are either an array with the first dimension the iteration, the
#   second the chain, and the third the parameter index; or an array with the
#   first dimension the iteration (across all chains) and the remaining
#   dimensions the parameter index (as in \code{\link[rstan]{extract}} when
#   \code{permuted=TRUE}).

make_samples_list <- function(samples, thin=1, concatenate.chains=FALSE,
                              verbose=FALSE, num_level=0){
    cat_v("Begin make_samples_list.\n", verbose, num_level=num_level)
    cat_v("Getting parameter names...", verbose, num_level=num_level+1)
    if (thin >= dim(samples)[1]){
        warning("thin is greater than number of samples; only the ",
                "first sample from each chain will be used")
    }
    samples.to.keep <- seq(1, dim(samples)[1], thin)
    samples.thin        <- array(samples[samples.to.keep, , ],
                                 dim=c(length(samples.to.keep), dim(samples)[-1]))
    dimnames(samples.thin) <- dimnames(samples)
    samples.param.names <- dimnames(samples)$parameters
    
    param.names  <- unlist(lapply(strsplit(samples.param.names, '\\['), '[', 1))
    unique.names <- unique(param.names)

    samples.list        <- vector(mode="list", length=length(unique.names))
    names(samples.list) <- unique.names
    cat_v(paste0("Found ", length(unique.names), " unique names.\n"),
                  verbose)
    
    if (concatenate.chains){
        cat_v("Concatenating chains ...", verbose,
                      num_level=num_level+1)
        samples.thin.conc <- samples.thin[, 1, ]
        if (dim(samples.thin)[2] > 1){
            for (i in 2:dim(samples.thin)[2]){
                samples.thin.conc <- rbind(samples.thin.conc,
                                           samples.thin[, i, ])
            }
        }
        cat_v("Done.\n", verbose)
    }
    cat_v("Getting per parameter...", verbose, num_level=num_level+1)
    for (param in unique.names){
        cat_v(paste0(param, "..."), verbose)
        param.idx <- grep(paste0("^", param), samples.param.names)
        if (!concatenate.chains){
            if (length(dim(samples.thin)) == 3){
                samples.list[[param]] <-
                    array(samples.thin[, , param.idx],
                          dim=c(dim(samples.thin)[1:2], length(param.idx)))
            }
        } else {
            samples.list[[param]] <- samples.thin.conc[, param.idx]
            is.matrix <- all(grepl(",", samples.param.names[param.idx]))
            is.vec    <- all(grepl("\\[", samples.param.names[param.idx]))
            if (is.matrix){
                samples.list[[param]] <-
                    array(samples.list[[param]],
                          dim=c(nrow(samples.list[[param]]),
                              sqrt(length(param.idx)), sqrt(length(param.idx))))
            } else if (!is.vec) {
                samples.list[[param]] <- array(samples.list[[param]],
                                               dim=length(samples.list[[param]]))
            }
        }
    }
    cat_v("Done.\n", verbose)
    cat_v("End make_samples_list.\n", verbose, num_level=num_level)
    return(samples.list)
}
