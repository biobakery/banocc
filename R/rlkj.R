# Generate a random sample from the LKJ distribution
#
# This function is based on code from
# \url{https://groups.google.com/forum/#!msg/stan-users/3gDvAs_qwN8/Xpgi2rPlx68J}.
# 
# @param d The dimension of the correlation matrix
# @param eta The scaling parameter of the LKJ distribution; must be > 1
#   (eta=1 means the distribution is uniform over d by d correlation matrices)
# @param cholesky Boolean: return the cholesky decomposition?

rlkj <- function(d, eta = 1, cholesky = FALSE) {
    if (d < 2){
        stop("Dimension of correlation matrix must be >= 2")
    }
    if (eta < 1){
        stop("The value of eta must be >= 1")
    }
    alpha <- eta + (d - 2) / 2
    L <- matrix(0, d, d)
    L[1,1] <- 1
    L[-1,1] <- partials <- rgbeta(d - 1, alpha)
    if(d == 2) {
      L[2,2] <- sqrt(1 - L[2,1]^2)
      if(cholesky) return(L)
      Sigma <- tcrossprod(L)
      return(Sigma)      
    }
    W <- log(1 - partials^2)
    for(i in 2:(d - 1)) {
      gap <- (i+1):d
      gap1 <- i:(d-1)
      alpha <- alpha - 0.5
      partials <- rgbeta(d - i, alpha)
      L[i,i] <- exp(0.5 * W[i-1])
      L[gap,i] <- partials * exp(0.5 * W[gap1])
      W[gap1] <- W[gap1] + log(1 - partials^2)
    }
    L[d,d] <- exp(0.5 * W[d-1])
    if(cholesky) return(L)
    Sigma <- tcrossprod(L)
    return(Sigma)      
  }


rgbeta <- function(d, shape) {
    if(shape == Inf)     rep(0, d)
    else if(shape > 0)  -1 + 2 * rbeta(d, shape, shape)
    else if(shape == 0) -1 + 2 * rbinom(d, 1, 0.5)
    else stop("shape must be non-negative")
  }
