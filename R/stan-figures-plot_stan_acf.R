#' Generate acf plots from a stanfit object
#' 
#' @param Fit A \code{stanfit} object
#' @param NumChains The number of chains used in the fitting
#' @param p The dimension of the correlation matrix
#' @param ln.params Boolean: generate plots of lognormal parameters?
#' @inheritParams rstan::sampling
#' @inheritParams cat_v
#' @inheritParams sample_mu_prior
#' @param samples.list An object as returned by \code{make_samples_list}.
#' @param is.matrix Boolean: is the parameter a matrix?
#' @param plot.lower Boolean: is the matrix symmetric?
#' @param plot.diag Boolean: include the diagnoal elements of a matrix?
#' @name plot_stan_acf
#' @return ACF plots are generated and printed on a per-parameter-element and
#'   per-chain basis.  Each parameter element is in a different plot; multiple
#'   chains are stacked side-by-side in a single plot.
#'
#' @importFrom rstan extract

#' @rdname plot_stan_acf
plot_stan_acf_set <- function(Fit, NumChains, p, ln.params=FALSE, thin=1,
                              verbose=FALSE, num_level=0){
    banocc::cat_v("Begin plot_stan_acf_set.\n", verbose, num_level=num_level)
    samples.list <-
        banocc::make_samples_list(rstan::extract(Fit, permuted=FALSE), thin=thin,
                                  verbose=verbose, num_level=num_level+1)
    if (ln.params){
        banocc::plot_stan_acf(samples.list, NumChains, p, "ln_L",
                             is.matrix = TRUE, plot.lower=TRUE, plot.diag=FALSE,
                             thin=thin, verbose=verbose, num_level=num_level+1)
        banocc::plot_stan_acf(samples.list, NumChains, p, "ln_mu",
                             is.matrix=FALSE, thin=thin, verbose=verbose,
                              num_level=num_level+1)
        banocc::plot_stan_acf(samples.list, NumChains, p, "ln_sigma",
                             is.matrix=FALSE, thin=thin, verbose=verbose,
                              num_level=num_level+1)
    } else {
        banocc::plot_stan_acf(samples.list, NumChains, p, "L", is.matrix = TRUE,
                             plot.lower=TRUE, plot.diag=FALSE, thin=thin,
                              verbose=verbose, num_level=num_level+1)
        banocc::plot_stan_acf(samples.list, NumChains, p, "mu", is.matrix=FALSE,
                             thin=thin, verbose=verbose, num_level=num_level+1)
        banocc::plot_stan_acf(samples.list, NumChains, p, "sigma",
                             is.matrix=FALSE, thin=thin, verbose=verbose,
                              num_level=num_level+1)
    }
    banocc::cat_v("End plot_stan_acf_set.\n", verbose, num_level=num_level)
}

#' @rdname plot_stan_acf
plot_stan_acf <- function(samples.list, NumChains, p, param.name, thin,
                          is.matrix, plot.lower=TRUE, plot.diag=FALSE,
                          verbose=FALSE, num_level=0){
    banocc::cat_v("Begin plot_stan_acf...", verbose, num_level=num_level)
    banocc::cat_v(paste0("for parameter ", param.name, "..."), verbose)
    par(mfrow=c(1, NumChains))
    if (is.matrix){
        if (plot.lower && !plot.diag){
            i.vals <- 1 + seq_len(p - 1)
        } else {
            i.vals <- seq_len(p)
        }
        for (i in i.vals){
            if (plot.lower && !plot.diag){
                j.vals <- 1:(i - 1)
            } else if (plot.lower){
                j.vals <- 1:i
            } else if (!plot.diag){
                j.vals <- setdiff(seq_len(p), i)
            } else {
                j.vals <- seq_len(p)
            }
            for (j in j.vals){
                for (k in seq_len(NumChains)){
                    acf(samples.list[[param.name]][, k, (j - 1) * p + i],
                        main="")
                    title(main=paste0(param.name, "[", i, ", ", j, "]\n",
                              "Chain = ", k, "; ",
                              "Thin = ", thin),
                          line=-1)
                }
            }
        }
    } else {
        for (i in seq_len(p)){
            for (k in seq_len(NumChains)){
                acf(samples.list[[param.name]][, k, i],
                    main="")
                title(main=paste0(param.name, "[", i, "]\n", "Chain = ", k, "; ",
                        "Thin = ", thin),
                      line=-1)
            }
        }
    }
    par(mfrow=c(1, 1))
    banocc::cat_v("Done.\n", verbose)
}
