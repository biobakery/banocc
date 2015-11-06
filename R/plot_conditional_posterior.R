#' Generate a conditional posterior plot based on a Stan fit
#'
#' @param Fit A \code{stanfit} object
#' @param estimates The posterior estimates (usually posterior median)
#' @param r.idx A vector of length two indicating the indices of the element in
#'   the correlation matrix whose plot will be generated
#' @param r The values of the correlation element over which to increment
#' @param ln Boolean: If TRUE, the parameter \code{ln_Rho} is plotted; if false,
#'   the parameter \code{Rho} is plotted.
#' @inheritParams sample_mu_prior
#'
#' @importFrom rstan log_prob
#' @importFrom rstan unconstrain_pars

plot_conditional_posterior <- function(Fit, estimates, r.idx, r=NULL, ln=FALSE,
                                       verbose=FALSE){

    if (verbose) cat("Getting initial R mat...")
    if (ln){
        R <- estimates$ln_Rho
    } else {
        R <- estimates$Rho
    }
    if (verbose) cat("Done.\n")

    if (is.null(r)){
        r <- seq(-1, 1, 0.01)
    }

    R.vals      <- list()
    idx.to.plot <- c()

    if (verbose) cat("Incrementing over r...")
    for(i in seq_along(r)){
        R.new.0 <- R
        R.new.0[r.idx[1], r.idx[2]] <- r[i]
        R.new.0[r.idx[2], r.idx[1]] <- r[i]
        if (ln){
            Sigma.ln <-
                diag(sqrt(diag(estimates$ln_Sigma))) %*% R.new.0 %*%
                    diag(sqrt(diag(estimates$ln_Sigma)))

            n.params   <- get_n_params(ln_mean = estimates$ln_mu,
                                       ln_Var = Sigma.ln)
            R.new <- cov2cor(n.params$n_Var)
        } else {
            R.new <- R.new.0
        }
        pos.def <- all(eigen(R.new, only.values=TRUE)$values > 0)

        if (pos.def){
            R.vals[[length(R.vals) + 1]] <- estimates

            L.new <- t(chol(R.new))
            R.vals[[length(R.vals)]]$L   <- L.new
            idx.to.plot <- c(idx.to.plot, i)
        }
    }
    if (verbose) cat("Done.\n")

    if (verbose) cat("Getting conditional log prob...")
    l.prob <- lapply(R.vals, function(pars){
        rstan::log_prob(Fit, rstan::unconstrain_pars(Fit, pars))
    })
    if (verbose) cat("Done.\n")

    if (verbose) cat("Plotting...")
    plot(r[idx.to.plot], unlist(l.prob), type="l",
         xlab=as.expression(substitute(paste(rho[id]),
                                       list(id=str_c(r.idx, collapse=",")))),
         ylab="conditional log prob")
    if (verbose) cat("Done.\n")
}
