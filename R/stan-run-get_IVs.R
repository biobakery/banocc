# Get initial values drawn from the prior
#
# @param chains The number of chains to run
# @param data A data object, which must be a list with elements:
#   \itemize{
#     \item \code{nu} The prior mean of mu
#     \item \code{Lambda} The prior variance/covariance matrix of mu
#     \item \code{eta} The parameter for the LKJ prior on the correlation
#       matrix
#     \item \code{alpha} The vector of shape parameters for the gamma prior on
#       sigma
#     \item \code{beta} The vector of scale parameters for the gamma prior on
#       sigma
#     \item \code{P} The number of features in the dataset
#   }
# @inheritParams cat_v
#

get_IVs <- function(chains, data, verbose=FALSE, num_level=0){
    cat_v("Begin get_IVs...", verbose, num_level=num_level)
    if (chains <= 0){
        stop("'chains' must be > 0")
    }
    IVs <- lapply(1:chains,
                  function(i) {
                      list(mu = as.vector(mvtnorm::rmvnorm(1, data$nu,
                               data$Lambda)),
                           L  = t(chol(banocc::rlkj(data$P, data$eta))),
                           sigma = mapply(function(a, b){
                               rgamma(1, shape=a, rate = b)
                               },
                               data$alpha, data$beta))
                  })
    cat_v("Done.\n", verbose)
    return(IVs)
}
