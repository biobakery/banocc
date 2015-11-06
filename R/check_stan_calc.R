#' Get as an R object any vectors or matrices printed by stan
#' 
#' @param stan_str A string of a single line of printed output from
#'   \code{\link[rstan]{sampling}} which contains a vector or matrix.
#'
#' @return Returns the vector or matrix as an R vector or matrix
#'
#' @examples
#' stan_str1 <- "mu = [-0.689748,-0.836567,0.815985]"
#' stan_str2 <- "Rho = [[1,-0.0601858,-0.708844],[-0.0601858,1,-0.385467],[-0.708844,-0.385467,1]]"
#' stan_str3 <- "mu_log_prob = -3.67753"
#'
#' get_R_from_stan(stan_str1)
#' get_R_from_stan(stan_str2)
#' get_R_from_stan(stan_str3)
#'
#' @importFrom stringr str_sub
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_split

get_R_from_stan <- function(stan_str){
  if (grepl("\\[", stan_str)){
    vec.part <- stringr::str_extract_all(stan_str, "[\\d\\.,\\[\\]\\-]+")[[1]]
    vec.no.brackets <- stringr::str_sub(tail(vec.part, 1), 2, -2)
    numeric.part    <- vec.no.brackets
  } else {
    split.equal <- stringr::str_split(stan_str, "=")[[1]]
    numeric.part <- tail(stringr::str_extract(split.equal, "[\\d\\.,\\-]+"), 1)
  }
  str.is.matrix <- grepl("\\[", numeric.part)
  if (str.is.matrix){
    numeric.part <- stringr::str_sub(numeric.part, 2, -2)
    split.by.row <- stringr::str_split(numeric.part, "\\],\\[")[[1]]
    numeric.vector <- as.numeric(unlist(stringr::str_split(split.by.row, ',')))
    numeric.matrix <- matrix(numeric.vector, byrow=TRUE,
                             nrow=length(split.by.row))
    return(numeric.matrix)
  } else {
    numeric.vector <- as.numeric(stringr::str_split(numeric.part, ",")[[1]])
    return(numeric.vector)
  }
}

#' Get the numeric value of the log-posterior from a stan string
#'
#' @param stan_str A string of a single line of printed output from
#'   \code{\link[rstan]{sampling}} which contains a log-posterior value
#'
#' @examples
#' stan_str <- "Adding L; lp__ = -5.87068"
#' get_stan_lp(stan_str)
#'
#' @importFrom stringr str_extract_all

get_stan_lp <- function(stan_str){
    numeric.str <-
        tail(stringr::str_extract_all(stan_str, "[0-9\\.\\-]+")[[1]], 1)
    return(as.numeric(numeric.str))
}

#' Calculate the components of the log-likelihood of C
#'
#' @param i The index of the subject whose log-likelihood is calculated.
#' @param C The data matrix: rows indicate samples, columns indicate features.
#' @param mu The current value of the normal mean.
#' @param Sigma The current value of the normal variance-covariance matrix.
#' 
#' @name ll_components

#' @rdname ll_components
norm.constant <- function(i, C){
    p <- ncol(C)
    return(-((p - 1) / 2) * log(2 * pi) - sum(log(C[i, ])))

}

#' @rdname ll_components
get_alpha <- function(i, mu, C){
    return(mu - log(C[i, ]))
}

#' @rdname ll_components
get_sigmasq.star <- function(Sigma, C){
    p <- ncol(C)
    j <- rep(1, p)
    return(1 / as.vector(t(j) %*% solve(Sigma) %*% j))
}

#' @rdname ll_components
get_mu.star <- function(i, mu, Sigma, C){
    p <- ncol(C)
    j <- rep(1, p)
    alpha <- get_alpha(i, mu, C)
    sigmasq.star <- get_sigmasq.star(Sigma, C)
    return(as.vector(t(j) %*% solve(Sigma) %*% alpha) * sigmasq.star)
}

#' @rdname ll_components
inc1 <- function(Sigma, C){
    sigmasq.star <- get_sigmasq.star(Sigma, C)
    return(0.5 * log(sigmasq.star))
}

#' @rdname ll_components
inc2 <- function(Sigma){
    return(-0.5 * log(det(Sigma)))
}

#' @rdname ll_components
inc3 <- function(i, mu, Sigma, C){
    alpha <- get_alpha(i, mu, C)
    return(as.vector(-0.5 * alpha %*% solve(Sigma) %*% alpha))
}

#' @rdname ll_components
inc4 <- function(i, mu, Sigma, C){
    mu.star <- get_mu.star(i, mu, Sigma, C)
    sigmasq.star <- get_sigmasq.star(Sigma, C)
    return(0.5 * mu.star^2 / sigmasq.star)
}

#' Check that the log-posterior as calculated by stan is approximately correct
#'
#' Prints out the change in log-posterior after each parameter subtracted from
#'   the value as calculated by an R function
#'
#' @param res The printed output from \code{\link[rstan]{sampling}} as a
#'   character vector; each element is a single line of output.
#' @param iter The iteration of stan from which to compare (since multiple
#'   iterations can be present in \code{res})
#'
#' @importFrom mvtnorm dmvnorm

check_stan_lik <- function(res, iter){
    one.iter.out <- res[grep("Starting lp__", res)[iter] + seq(0, 23)]
    cat(str_c(one.iter.out, collapse="\n"))
    mu    <- banocc::get_R_from_stan(one.iter.out[grep("^mu =", one.iter.out)])
    sigma <- banocc::get_R_from_stan(one.iter.out[grep("^sigma =",
                                                       one.iter.out)])
    L     <- banocc::get_R_from_stan(one.iter.out[grep("^L =", one.iter.out)])
    
    cat("\n\nChecking mu prior...\n")
    mu.prior <- mvtnorm::dmvnorm(mu, Data$nu, Data$Lambda, log=TRUE) +
        (Data$P / 2) * log(2 * pi)
    ## start.lp is non-zero b/c it is the Jacobian of transforming the
    ##  constrained parameters to unconstrained space
    start.lp <- banocc::get_stan_lp(one.iter.out[grep("^Starting",
                                                      one.iter.out)])
    mu.lp    <- banocc::get_stan_lp(one.iter.out[grep("^Adding mu",
                                                      one.iter.out)])
    cat("stan mu prior - round(mu prior, 5):\n")
    print((mu.lp - start.lp) - round(mu.prior, 5))
    
    cat("\nChecking L prior...\n")
    L.lp <- banocc::get_stan_lp(one.iter.out[grep("^Adding L", one.iter.out)])
    L.prior <- sum(sapply(2:Data$P, function(i){
        (Data$P - i + 2 * Data$eta - 2) * log(L[i, i])
    }))
    cat("stan L prior - round(L prior, 5):\n")
    print((L.lp - mu.lp) - round(L.prior, 5))
    
    cat("\nChecking sigma prior...\n")
    sigma.lp <- banocc::get_stan_lp(one.iter.out[grep("^Adding sigma",
                                                      one.iter.out)])
    sigma.prior <- sum(sapply(1:Data$P, function(i){
        dgamma(sigma[i], Data$alpha[i], Data$beta[i], log=TRUE) - 
            Data$alpha[i] * log(Data$beta[i]) + log(gamma(Data$alpha[i]))
    }))
    cat("stan sigma prior - round(sigma prior, 5):\n")
    print((sigma.lp - L.lp) - round(sigma.prior, 5))
    
    cat("\nChecking Sigma...\n")
    Sigma <- diag(sigma) %*% L %*% t(L) %*% diag(sigma)
    cat("stan Sigma - round(Sigma, 7):\n")
    print(banocc::get_R_from_stan(one.iter.out[grep("^Sigma = ",
                                                    one.iter.out)]) - 
        round(Sigma, 7))

    cat("\nSigma eigen values:\n")
    print(eigen(Sigma)$values)
    
    cat("\nChecking Sigma inverse...\n")
    Sigma.inv <- solve(Sigma)
    cat("stan Sigma inverse - round(Sigma inverse, 5):\n")
    print(banocc::get_R_from_stan(one.iter.out[grep("^invSigma = ",
                                                    one.iter.out)]) - 
        round(Sigma.inv, 5))
    
    cat("\nChecking sigmasq.star...\n")
    sigmasq.star <- banocc::get_sigmasq.star(Sigma, Data$C)
    cat("stan sigmasq.star - round(sigmasq.star, 8):\n")
    print(banocc::get_R_from_stan(one.iter.out[grep("^sigma_sq_star =",
                                                    one.iter.out)]) - 
        round(sigmasq.star, 8))

    cat("\nChecking mu.star...\n")
    mu.star <- sapply(1:Data$N, banocc::get_mu.star, mu=mu, Sigma=Sigma,
                      C=Data$C)
    cat("stan mu.star - round(mu.star, 8):\n")
    print(banocc::get_R_from_stan(one.iter.out[grep("^mu_star =",
                                                    one.iter.out)]) - 
        round(mu.star, 5))

    cat("\nChecking likelihood increments...\n")
    inc.1 <- sum(sapply(1:Data$N, function(i) banocc::inc1(Sigma, Data$C)))
    cat("stan inc.1 - round(inc.1, 4):\n")
    print(banocc::get_R_from_stan(one.iter.out[grep("^inc_1 =",
                                                    one.iter.out)]) - 
        round(inc.1, 4))
    
    inc.2 <- sum(sapply(1:Data$N, function(i) banocc::inc2(Sigma)))
    cat("stan inc.2 - round(inc.2, 4):\n")
    print(banocc::get_R_from_stan(one.iter.out[grep("^inc_2 =",
                                                    one.iter.out)]) - 
        round(inc.2, 4))
    
    cat("stan inc.1 + inc.2 likelihood - round(inc.1 + inc.2 likelihood, 5):\n")
    inc2.lp <- banocc::get_stan_lp(one.iter.out[grep("^Added inc_1",
                                                     one.iter.out)])
    print((inc2.lp - sigma.lp) - round(inc.1 + inc.2, 5))
    
    cat("stan inc.3 - round(inc.3, 2):\n")
    inc.3 <- sapply(1:Data$N, banocc::inc3, mu=mu, Sigma=Sigma, C=Data$C)
    print(banocc::get_R_from_stan(one.iter.out[grep("^inc_3i =",
                                                    one.iter.out)]) - 
        round(inc.3, 2))
    
    cat("stan inc.4 - round(inc.4, 2):\n")
    inc.4 <- sapply(1:Data$N, banocc::inc4, mu=mu, Sigma=Sigma, C=Data$C)
    print(banocc::get_R_from_stan(one.iter.out[grep("^inc_4i =",
                                                    one.iter.out)]) - 
        round(inc.4, 2))
    
    cat("stan inc.3 + inc.4 likelihood - (inc.3 + inc.4 likelihood):\n")
    inc4.lp <- banocc::get_stan_lp(one.iter.out[grep("^Added inc_3",
                                                     one.iter.out)])
    print((inc4.lp - inc2.lp) - (sum(inc.3) + sum(inc.4)))
}
