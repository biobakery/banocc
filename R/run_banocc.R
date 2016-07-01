#' Runs BAnOCC to fit the model and generate appropriate convergence metrics
#'   and inference.
#'
#' @param banocc_model The compiled stan model (as with
#'   \code{stan_model(model_code = banocc_model)}).
#' @param C The dataset (rows=samples, columns=features), which is nxp
#' @param n The prior mean for m; vectors of length less than p will be
#'   recycled.
#' @param L The prior variance-covariance for m (must be
#'   positive-definite with dimension pxp where p=number of features in C), or
#'   a vector of length p of variances for m. If a vector of length less than p
#'   is given, it will be recycled.
#' @param a,b Vectors of parameters for the gamma prior on the standard
#'   deviations. They must be of equal length. If a vector of length less than p
#'   is given, it will be recycled.
#' @param sd_mean,sd_var Vectors of means and variances for the gamma priors on
#'   the standard deviations. They must be of equal length. If a vector of length
#'   less than p is given, it will be recycled.
#' @param eta scale parameter of LKJ distribution; must be >= 1.
#' @param init The initial values as a list (see
#'   \code{\link[rstan]{sampling}}). Default value is NULL, which means that
#'   initial values are sampled from the priors.
#' @inheritParams rstan::sampling
#' @param get_min_width A boolean value: should the minimum CI width that
#'   includes zero be calculated?
#' @param conf_alpha The percentage of the posterior density outside the
#'   credible interval. That is, a \code{1-conf_alpha} * 100\% credible
#'   interval will be returned.
#' @param calc_snc Boolean: should the scaled neighborhood criterion be
#'   calculated?
#' @param verbose Print informative statements as the function executes?
#' @param num_level The number of the level (will determine the number of
#'   spaces to add to verbose output)
#'
#' @importFrom rstan extract
#' @importFrom rstan sampling
#'
#' @export

run_banocc <- function(banocc_model, C, n = rep(0, ncol(C)),
                       L = 10*diag(ncol(C)),
                       a = rep(1, ncol(C)), b = rep(0.5, ncol(C)),
                       eta = 1, cores = getOption("mc.cores", 1L),
                       chains = 4, iter = 50, warmup = floor(iter/2),
                       thin = 1, init = NULL, control=NULL,
                       sd_mean=NULL, sd_var=NULL, conf_alpha=0.05,
                       get_min_width=FALSE, calc_snc=FALSE, verbose=FALSE,
                       num_level=0){
    cat_v("Begin run_banocc\n", verbose, num_level=num_level)
    C <- check_C(C, verbose=verbose, num_level=num_level+1)
    Data <- list(C=C, N=nrow(C), P=ncol(C))
    
    Data$n <- check_n(n, Data$P, verbose, num_level=num_level+1)
    Data$L <- check_L(L, Data$P, verbose,
                      num_level=num_level+1)
    a_b <- get_a_b(a=a, b=b,
                   sd_mean=sd_mean, sd_var=sd_var,
                   p=Data$P,
                   verbose=verbose, num_level=num_level+1)
    Data$a   <- a_b$a
    Data$b   <- a_b$b
    Data$eta <- get_eta(eta)

    if (is.null(init)){
        init <- get_IVs(chains=chains, data=Data, verbose=verbose,
                                num_level=num_level + 1)
    }

    cat_v("Begin fitting the model\n", verbose, num_level=num_level+1)
    refresh <- ifelse(verbose, max(iter/10, 1), 0)
    Fit <- rstan::sampling(banocc_model, data=Data,
                           chains=chains, iter=iter,
                           warmup=warmup, thin=thin,
                           init=init, cores=cores,
                           control=control,
                           show_messages=FALSE,
                           refresh=refresh)
    cat_v("End fitting the model\n", verbose, num_level=num_level+1)

    post.samples.list <- rstan::extract(Fit)
    CI <- get_credible_intervals(posterior_samples=post.samples.list,
                                         list=TRUE,
                                         parameter.names=c("W"),
                                         conf=1-conf_alpha,
                                         type="marginal.hpd",
                                         verbose=verbose, num_level=num_level+1)

    dimnames(CI$W$lower) <- list(colnames(Data$C), colnames(Data$C))
    dimnames(CI$W$upper) <- list(colnames(Data$C), colnames(Data$C))
    CI <- CI$W

    if (get_min_width){
        min_width <- get_min_width(posterior_sample=post.samples.list,
                                           parameter.names=c("W"),
                                           null_value=0, type="marginal.hpd",
                                           precision=0.01, verbose=verbose,
                                           num_level=num_level + 1)
    } else {
        min_width <- list(Rho=NULL)
    }
    min_width <- min_width$W

    if (calc_snc){
        snc <- get_snc(posterior_samples=post.samples.list,
                               parameter.names=c("W"))
    } else {
        snc <- list(Rho=NULL)
    }
    snc <- snc$W

    Estimates <-
        get_posterior_estimates(posterior_samples=post.samples.list,
                                        estimate_method="median",
                                        parameter.names="W")
    dimnames(Estimates$W) <- list(colnames(Data$C), colnames(Data$C))
    Estimates <- Estimates$W

    
    return_object <- list(Data=Data, Fit=Fit, 
                          CI.hpd=CI, Estimates.median=Estimates,
                          Min.width=min_width, SNC=snc)

    cat_v("End run_banocc\n", verbose, num_level=num_level)

    return(return_object)
}

check_n <- function(n, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_n\n", verbose, num_level=num_level)
    n <- as.numeric(n)
    n <- check_vector("n", n, p, verbose, num_level=num_level+1)
    cat_v("End check_n\n", verbose, num_level=num_level)
    return(n)
}

check_a_b <- function(a, b, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_a_b\n", verbose, num_level=num_level)
    if (length(a) != length(b)){
        stop("'a' and 'b' must be of equal length")
    }
    a <- as.numeric(a)
    a <- check_vector("a", a, p, verbose,
                                  num_level=num_level + 1)
    if (any(a <= 0)){
        stop("'a' values must be positive")
    }
    b  <- as.numeric(b)
    b  <- check_vector("b", b, p, verbose,
                                  num_level=num_level + 1)
    if (any(b <= 0)){
        stop("'b' values must be positive")
    }
    cat_v("End check_a_b\n", verbose, num_level=num_level)
    return(list(a=a, b=b))
}

check_sd_mean_var <- function(sd_mean, sd_var, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_sd_mean_var\n", verbose, num_level=num_level)
    if (length(sd_mean) != length(sd_var)){
        stop("'sd_mean' and 'sd_var' must be of equal length")
    }
    sd_mean <- as.numeric(sd_mean)
    sd_mean <- check_vector("sd_mean", sd_mean, p, verbose,
                                    num_level=num_level + 1)
    if (any(sd_mean <= 0)){
        stop("'sd_mean' values must be positive")
    }
    sd_var  <- as.numeric(sd_var)
    sd_var  <- check_vector("sd_var", sd_var, p, verbose,
                                    num_level = num_level + 1)
    if (any(sd_var <= 0)){
        stop("'sd_var' values must be positive")
    }
    cat_v("End check_sd_mean_var\n", verbose, num_level=num_level)
    return(list(sd_mean=sd_mean, sd_var=sd_var))
}

check_vector <- function(parm.name, parm, p, verbose=TRUE, num_level=0){
    cat_v("Begin check_vector...", verbose, num_level=num_level)
    if (!is.vector(parm) || mode(parm)!="numeric"){
        stop(paste0("'", parm.name, "' must be a vector"))
    }
    if ((length(parm) < p) && (length(parm) > 1)){
        warning(paste0("recycling '", parm.name, "'"))
        parm <- rep(parm, ceiling(p / length(parm)))[1:p]
    } else if (length(parm) == 1){
        parm <- rep(parm, p)
    } else if (length(parm) > p){
        warning("length of '", parm.name,
                "' is > p; only using first p elements")
        parm <- parm[1:p]
    }
    cat_v("Done.\n", verbose)
    return(parm)
}

check_L <- function(L, p, verbose=FALSE, num_level=0){
    cat_v("Begin check_L...", verbose, num_level=num_level)
    if (!is.numeric(L)){
        stop("L must be numeric")
    }

    if (is.matrix(L)){
        if (any(dim(L) != p)){
            stop(paste0("L must be a square matrix with the same number of",
                        " columns as C"))
        }
        if (any(eigen(L)$values <= 0)){
            stop("'L' is not positive definite.")
        }
        if (!is.symmetric(L)){
            stop("'L' must be symmetric")
        }
    } else if (is.vector(L)){
        if ((length(L) < p) && (length(L) > 1)){
            warning("'L' is being recycled")
            num_rep <- ceiling(p / length(L))
            L <- diag(rep(L, num_rep)[1:p])
        } else if (length(L) == 1){
            L <- diag(rep(L, p))
        } else if (length(L) > p){
            warning("'L' has length > p; only first p elements will be ",
                    "used.")
            L <- diag(L[1:p])
        } else {
            L <- diag(L)
        }
    } else {
        stop("'L' must be either a matrix or a vector")
    }
    cat_v("Done.\n", verbose)
    return(L)
}

get_a_b <- function(a, b, p, sd_mean=NULL, sd_var=NULL,
                           verbose=FALSE, num_level=0){
    cat_v("Begin get_a_b\n", verbose=verbose,
                  num_level=num_level)
    if (!is.null(a) || !is.null(b)){
        if (is.null(a) || is.null(b)){
            stop(paste0("Must provide both 'a' and 'b' OR both 'sd_mean'",
                        " and 'sd_var'"))
        }
        a_b <- check_a_b(a, b, p, verbose,
                                               num_level=num_level+1)
    } else if (!is.null(sd_var) || !is.null(sd_mean)){
        if (is.null(sd_var) || is.null(sd_mean)){
            stop(paste0("Must provide both 'sd_mean' and 'sd_var' OR both ",
                        "'a' and 'b'"))
        }
        a_b <- list()
        sd_mean_var <- check_sd_mean_var(sd_mean, sd_var, p,
                                                 num_level=num_level+1)
        a_b$b  <- sd_mean_var$sd_mean / sd_mean_var$sd_var
        a_b$a <- sd_mean_var$sd_mean^2 / sd_mean_var$sd_var
    } else {
        stop(paste0("Must provide both 'a' and 'b' OR both 'sd_mean'",
                    " and 'sd_var'"))
    }
    cat_v("End get_a_b\n", verbose)
    return(a_b)
}

get_eta <- function(eta){
    if (eta < 1){
        stop("'eta' must be >= 1")
    } else {
        return(eta)
    }

}

check_C <- function(C, zero_adj=0.0001, verbose=FALSE, num_level=0){
    cat_v("Begin checking input data matrix\n", verbose=verbose,
          num_level=num_level)
    if (any(C < 0)){
        stop("Some values of C are < 0")
    }
    if (any(C > 1)){
        stop("Some values of C are > 1")
    }
    if (any(rowSums(C) - 1 > 1e-8)){
        stop(paste0("Some row sums of C are not 1; perhaps you transposed ",
                    "features and samples?"))
    }
    if (!is.data.frame(C) && !is.matrix(C)){
        stop("C must be a data frame or matrix")
    }
    C <- as.matrix(C)
    if (!is.numeric(C)){
        stop("C must be numeric")
    }
    if (any(C == 0)){
        warning(paste0("Some values of C are zero. ",
                       "Since zero-inflation is not yet implemented, these ",
                       "will be changed to ", zero_adj, "."))
        C <- adjust_zeros(C, zero_adj=zero_adj, verbose=verbose,
                          num_level=num_level+1)
    }
    cat_v("End checking input data matrix\n", verbose=verbose,
          num_level=num_level)
    return(C)
}
