# Evaluate convergence of multiple chains using the 'Rhat' statistic
#
# @param b_stanfit The \code{stanfit} object from \code{rstan::sampling}.
# @inheritParams run_banocc
#

evaluate_convergence <- function(b_stanfit, verbose=FALSE, num_level=0){
    cat_v("Begin evaluating convergence\n", verbose,
          num_level=num_level+1)
    rhat_stat <- rstan::summary(b_stanfit)$summary[, "Rhat"]
    diag_elts <- grep("W.*\\[([0-9]*),[ ]?\\1\\]", names(rhat_stat))
    rhat_stat <- rhat_stat[-diag_elts]
    if (any(is.na(rhat_stat)) || max(rhat_stat) > 1.2){
        fit_converged <- FALSE
        warning(paste0("Fit has not converged as evaluated by the Rhat ",
                       "statistic. You might try a larger number of ",
                       "warmup iterations, different priors, or ",
                       "different initial values. See vignette for ",
                       "more on evaluating convergence."))
    } else {
        fit_converged <- TRUE
    }
    cat_v("End evaluating convergence\n", verbose, num_level=num_level+1)
    return(fit_converged)
}
