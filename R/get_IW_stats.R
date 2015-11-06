#' Calculate summary statistics (mean, mode and variance) for an
#'   inverse-Wishart distribution
#'
#' @inheritParams MCMCpack::riwish
#' 
#' @return A named list: \code{Mean} is the mean, \code{Variance} is the
#'   variance, and \code{Mode} is the mode.

get_IW_stats <-
function(S, v){
  p <- nrow(S)
  IW.stats <- list()
  if(v > p + 1){
    IW.mean <- S/(v - p - 1)
    IW.stats$Mean <- IW.mean
  }
  IW.mode <- S/(v + p + 1)
  IW.stats$Mode <- IW.mode
  
  if(v > p + 3){
    denom  <- (v - p) * ((v - p - 1)^2) * (v - p - 3)
    num1   <- (v - p + 1) * S^2
    num2   <- (v - p - 1) * diag(S) %*% t(diag(S))
    IW.Var <- (num1 + num2) / denom
    IW.stats$Variance <- IW.Var
  }

  return(IW.stats)
}
