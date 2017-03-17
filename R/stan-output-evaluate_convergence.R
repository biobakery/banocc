# Evaluate convergence of multiple chains using the 'Rhat' statistic
#
# @param b_stanfit The \code{stanfit} object from \code{rstan::sampling}.
# @inheritParams run_banocc
#

evaluate_convergence <- function(b_stanfit, verbose=FALSE, num_level=0){
    cat_v("Begin evaluating convergence\n", verbose,
          num_level=num_level)
    rhat_stat <- rstan::summary(b_stanfit)$summary[, "Rhat"]
    if (any(is.na(rhat_stat)) || max(rhat_stat) > 1.2){
        fit_converged <- FALSE
        warning("Fit has not converged as evaluated by the Rhat ",
                "statistic. You might try a larger number of ",
                "warmup iterations, different priors, or ",
                "different initial values. See vignette for ",
                "more on evaluating convergence.")
    } else {
        fit_converged <- TRUE
    }
    cat_v("End evaluating convergence\n", verbose, num_level=num_level)
    return(fit_converged)
}
