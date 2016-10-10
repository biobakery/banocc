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
#' @importFrom mvtnorm rmvnorm
#' @importFrom rmutil rlaplace

get_IVs <- function(chains, data, verbose=FALSE, num_level=0){
    cat_v("Begin get_IVs...", verbose, num_level=num_level)
    if (chains <= 0){
        stop("'chains' must be > 0")
    }
    IVs <- lapply(1:chains,
                  function(i) {
                      list(m = as.vector(mvtnorm::rmvnorm(1, data$n,
                               data$L)),
                           O = sample_O(data))
                  })
    cat_v("Done.\n", verbose)
    return(IVs)
}

sample_O <- function(data, num_tries=10){
    O <- matrix(0, ncol=data$P, nrow=data$P)
    diag(O) <- rexp(data$P, data$lambda/2)
    ntries <- 0
    while(any(eigen(O)$values <= 0) && ntries < num_tries){
      O[upper.tri(O)] <- rlaplace(choose(Data$P, 2), s=data$lambda)
      O[lower.tri(O)] <- t(O)[lower.tri(t(O))]        
    }

    if (all(eigen(O)$values > 0)) return(O)
    stop(paste0("After ", num_tries, "tries, could not sample PD O. ",
         "Try decreasing labda or specifying the initial values by hand."))
}
