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

get_IVs <- function(chains, data, verbose=FALSE, num_level=0){
    cat_v("Begin get_IVs...", verbose, num_level=num_level)
    if (chains <= 0){
        stop("'chains' must be > 0")
    }
    IVs <- lapply(1:chains,
                  function(i) {
                      list(m = as.vector(mvtnorm::rmvnorm(1, data$n,
                               data$L)),
                           O = sample_O(data),
                           lambda = rgamma(1, data$a, data$b))
                  })
    cat_v("Done.\n", verbose)
    return(IVs)
}

sample_O <- function(data, num_tries=10){
    lambda <- max(0.001, data$a / data$b - 2 * (data$a / data$b^2))
    O <- matrix(0, ncol=data$P, nrow=data$P)
    diag(O) <- rexp(data$P, lambda/2)
    O <- sample_O_tri(data, lambda, O)
    ntries <- 1
    while(any(eigen(O)$values <= 0) && ntries < num_tries){
      O <- sample_O_tri(data, lambda, O)
      ntries <- ntries + 1
    }

    if (all(eigen(O)$values > 0)) return(O)
    stop(paste0("After ", num_tries, " tries, could not sample PD O. ",
         "Try decreasing lambda or specifying the initial values by hand."))
}

sample_O_tri <- function(data, lambda, O){
    expsamp <- rexp(choose(data$P, 2), rate=1/lambda)
    O[upper.tri(O)] <- expsamp * sample(c(-1, 1), choose(data$P, 2),
                                        replace=TRUE)
    O[lower.tri(O)] <- t(O)[lower.tri(t(O))]        
    return(O)
}
