# Get credible intervals for the output according to different paradigms
#
# @param precision The precision used to estimate the minimum width.
# @param null_value The null value to which to compare the credible interval.
#   Must be the same for all parameters.
# @inheritParams get_posterior_quantiles
# @inheritParams get_credible_intervals
# @inheritParams cat_v
# @return Returns a list of the minimum CI widths that include
#   \code{null_value} for each parameter in \code{parameter.names}.
#
get_min_width <- function(posterior_samples, parameter.names="W",
                          null_value=0, type="marginal.hpd", precision=0.01,
                          verbose=FALSE, num_level=0){
    cat_v("Begin get_min_width\n", verbose,
                  num_level=num_level)

    samples_dim <- lapply(posterior_samples, dim)
    min_width <- initialize_min_width(
        parameter.names=parameter.names, samples_dim=samples_dim,
        value=NA, verbose=verbose, num_level=num_level + 1)

    missing_min_width <- initialize_min_width(
        parameter.names=parameter.names, samples_dim=samples_dim,
        value=TRUE, verbose=verbose, num_level=num_level + 1)

    width <- precision
    any_still_missing <- Reduce(any, lapply(missing_min_width, any,
                                            na.rm=TRUE))


    while(any_still_missing && width < 1-precision){
        updated_width <- eval_width(
            min_width=min_width, missing_min_width=missing_min_width,
            posterior_samples=posterior_samples,
            parameter.names=parameter.names, width=width, type=type,
            null_value=null_value, verbose=verbose, num_level=num_level+1)
        min_width <- updated_width$min_width
        missing_min_width <- updated_width$missing_min_width

        width <- width + precision

        any_still_missing <- Reduce(any, lapply(missing_min_width, any,
                                                na.rm=TRUE))
    }

    if (any_still_missing){
        min_width <- update_min_width(
            min_width=min_width, which_to_update=missing_min_width,
            width=1, verbose=verbose, num_level=num_level+1)
    }

    cat_v("End get_min_width\n", verbose, num_level=num_level)
    return(min_width)
}

# Evaluate everything for a particular width and update min_width and
#   missing_min_width lists accordingly.
# 
# @param min_width A list; each elt is a parameter. The values in an element
#   are the minimum widths.
# @param missing_min_width The same structure as min_width, but the values
#   in an element are indicators of whether the min_width has been achieved.
# @param width The current width of credible interval to evaluate
# @param null_value The null value of the ``hypothesis test''
# @inheritParams get_credible_intervals
# @inheritParams cat_v
# 
eval_width <- function(min_width, missing_min_width, posterior_samples,
                       parameter.names, width, type, null_value,
                       verbose=FALSE, num_level=0
                       ){
    cat_v("Begin eval_width\n", verbose,
                  num_level=num_level)
    CI <- get_credible_intervals(posterior_samples=posterior_samples,
                                         list=FALSE,
                                         parameter.names=parameter.names,
                                         conf=width, type="marginal.hpd",#type,
                                         verbose=verbose,
                                         num_level=num_level+1)

    CI_include_null <- check_CI(CI=CI, null_value=null_value,
                                verbose=verbose, num_level=num_level+1)

    any_included <- Reduce(any, lapply(CI_include_null, any, na.rm=TRUE))

    if (any_included){
        which_to_update <- find_idx_to_update(
            CI_include_null=CI_include_null,
            missing_min_width=missing_min_width)
        need_to_update <- Reduce(any, lapply(which_to_update, any,
                                             na.rm=TRUE))
        if (need_to_update){
            updated_min_width <- update_min_width(
                min_width=min_width, which_to_update=which_to_update,
                width=width,
                verbose=verbose, num_level=num_level+1)
            updated_missing_min_width <- update_missing_min_width(
                missing_min_width=missing_min_width,
                which_to_update=which_to_update, width=width,
                verbose=verbose,
                num_level=num_level+1)
        } else {
            updated_min_width <- min_width
            updated_missing_min_width <- missing_min_width
        }
    } else {
        updated_min_width <- min_width
        updated_missing_min_width <- missing_min_width
    }

    cat_v("End eval_width\n", verbose,
                  num_level=num_level)
    return(list(min_width=updated_min_width,
                missing_min_width=updated_missing_min_width))
}

# Update the min_width list. Replaces every element to update with width.
update_min_width <- function(min_width, which_to_update, width,
                             verbose=FALSE, num_level=0){
    cat_v("Begin update_min_width.\n", verbose, num_level=num_level)
    updated_min_width <- lapply(names(min_width), function(name){
        new_width <- min_width[[name]]
        new_width[which_to_update[[name]]] <- width
        if (is.matrix(min_width[[name]])){
            matrix(new_width, ncol=ncol(min_width[[name]]),
                   nrow=nrow(min_width[[name]]))
        } else {
            new_width
        }
    })
    names(updated_min_width) <- names(min_width)
    cat_v("End update_min_width.\n", verbose, num_level=num_level)
    return(updated_min_width)
}

# Updated the missing_min_width list; elements are set to FALSE if (1)
#   they've just had a width added or (2) they have a missing width b/c
#   the CI width is always zero (eg., correlation diagonal elts.)
update_missing_min_width <- function(missing_min_width, which_to_update,
                                     width,
                                     verbose=FALSE, num_level=0){
    cat_v("Begin update_missing_min_width.\n", verbose,
                  num_level=num_level)
    updated_missing_min_width <- lapply(names(missing_min_width), function(name){
        new <- missing_min_width[[name]]
        new[which_to_update[[name]]] <- FALSE
        if (width >= 0.9 && any(is.na(which_to_update[[name]]))){
            new[is.na(which_to_update[[name]])] <- FALSE
        }
        if (is.matrix(missing_min_width[[name]])){
            matrix(new, ncol=ncol(missing_min_width[[name]]),
                   nrow=nrow(missing_min_width[[name]]))
        } else {
            new
        }
    })
    names(updated_missing_min_width) <- names(missing_min_width)
    cat_v("End update_missing_min_width.\n", verbose,
                  num_level=num_level)
    return(updated_missing_min_width)
}

# Find which elements need to have the min_width updated. These are those
#   whose CIs include the null (CI_include_null) AND haven't yet
#   been updated (missing_min_width)
find_idx_to_update <- function(CI_include_null, missing_min_width,
                               verbose=FALSE, num_level=0){
    cat_v("Begin find_idx_to_update.\n", verbose, num_level=num_level)
    which_to_update <- lapply(names(CI_include_null), function(name){
        update <- as.logical(CI_include_null[[name]] &
                             missing_min_width[[name]])
        if (is.vector(CI_include_null[[name]])){
            update
        } else {
            matrix(update, ncol=ncol(CI_include_null[[name]]),
                   nrow=nrow(CI_include_null[[name]]))
        }
    })
    names(which_to_update) <- names(CI_include_null)
    cat_v("End find_idx_to_update.\n", verbose, num_level=num_level)
    return(which_to_update)
}

# Checks a credible interval list (each elt is an array) against
#   a null value; returns a list. Each value in a list elt is
#   TRUE if the CI includes null_value, FALSE if it doesn't and NA
#   if the CI width is 0
check_CI <- function(CI, null_value,  verbose=FALSE, num_level=0){
    cat_v("Begin check_CI\n", verbose, num_level=num_level)
    CI_include_null <- lapply(names(CI), function(name){
        is.mat <- length(dim(CI[[name]])) == 3
        if (is.mat){
            apply_dim <- c(2, 3)
        } else {
            apply_dim <- 2
        }
        apply(CI[[name]], apply_dim, function(ci){
            if (isTRUE(all.equal(ci[1], ci[2],
                                 check.attributes=FALSE))){
                NA
            } else {
                (ci[1] < null_value) && (ci[2] > null_value)
            }
        })

    })
    names(CI_include_null) <- names(CI)

    cat_v("End check_CI\n", verbose, num_level=num_level)
    return(CI_include_null)
}

# Initializes a min_width list with names as parameter.names and
#   value in each entry of vector or matrix, as appropriate (determined by
#   samples_dim)
initialize_min_width <- function(parameter.names, samples_dim, value,
                                 verbose=verbose, num_level=0){
    
    cat_v("Begin initialize_min_width...", verbose,
                  num_level=num_level)
    min_width_init <- lapply(parameter.names, function(name){
        if (length(samples_dim[[name]]) == 3){
            matrix(value, ncol=samples_dim[[name]][2],
                   nrow=samples_dim[[name]][3])
        } else {
            rep(value, samples_dim[[name]][2])
        }
    })
    names(min_width_init) <- parameter.names
    cat_v("Done.\n", verbose)
    return(min_width_init)
}
