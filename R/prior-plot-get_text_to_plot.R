#' get_text_to_plot
#'
#' @inheritParams get_Matrix_to_plot
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @param digits The number of digits to which to round summary statistics
#'   (if plotted)
#' @param center.stat One of \code{"mean"} or \code{"median"}; the summary
#'   statistic of central tendency to be plotted.
#' @return A list of summary statistics for each plot (determined by unique
#'   indicies of the parameter). Each element of the list is a data frame with
#'   columns \code{i} (and \code{k} in the case of \code{get_Matrix_text})
#'   indicating the index, \code{stat} indicating the value of the statistic, and
#'   \code{stat.type} indicating the name of the statistic (e.g., median or
#'   variance)
#' @name get_text_to_plot

#' @rdname get_text_to_plot
get_Matrix_text <-
function(Matrix.sample, digits=2, diag=TRUE, center.stat="mean", verbose=FALSE,
         num_level=0){
  if (!(center.stat %in% c("mean", "median"))){
      stop("center.stat must be one of: mean, median")
  }
  banocc::cat_v("Begin get_Matrix_text\n", verbose, num_level=num_level)
  p <- nrow(Matrix.sample[[1]])
  
  banocc::cat_v("Getting matrix parameter stats...", verbose,
                num_level=num_level+1)
  Matrix_text <- vector("list", length=(choose(p, 2) * 2 + p * diag))
  idx <- 0
  next.indices <- 1:p
  for (i in 1:p){
    Matrix.i <- Reduce(rbind, lapply(Matrix.sample, '[', i, ))
    if (!diag){
      next.indices <- setdiff(1:p, i)
    }
    for (k in next.indices){
      idx <- idx + 1
      x <- Matrix.i[,k]
      Matrix_text[[idx]] <-
              data.frame(i=rep(i, 2), k=rep(k, 2),
                         x=rep(0.3*max(unlist(Matrix.i)), 2))
      if(center.stat=="mean"){
          Matrix_text[[idx]]$stat <-
              c(signif(mean(x), digits), signif(var(x), digits))
          Matrix_text[[idx]]$stat.type <- c("mean", "variance")
      } else {
          Matrix_text[[idx]]$stat <-
              c(signif(median(x), digits), signif(var(x), digits))
          Matrix_text[[idx]]$stat.type <- c("median", "variance")
      }
    }
  }
  banocc::cat_v("Done.\n", verbose)
  banocc::cat_v("End get_Matrix_text\n", verbose, num_level=num_level)
  return(Matrix_text)
}

#' @rdname get_text_to_plot
get_vec_text <-
function(vec.sample, digits=2, center.stat="mean", verbose=FALSE,
         num_level=0){
  if (!(center.stat %in% c("mean", "median"))){
      stop("center.stat must be one of: mean, median")
  }
  banocc::cat_v("Begin get_vec_text\n", verbose, num_level=num_level)
  p <- length(vec.sample[[1]])
  
  banocc::cat_v("Getting vector parameter stats...", verbose,
                num_level=num_level+1)
  vec_text <- vector("list", length=p)
  for (i in 1:p){
    vec.i <- unlist(lapply(vec.sample, '[', i))
    vec_text[[i]] <- 
      data.frame(i=rep(i, 2),
                 x=rep(0.5*max(vec.i), 2)
                 )
    if (center.stat=="mean"){
        vec_text[[i]]$stat <-
            c(signif(mean(vec.i), digits), signif(var(vec.i), digits))
        vec_text[[i]]$stat.type <- c("mean", "variance")
    } else {
        vec_text[[i]]$stat <-
            c(signif(median(vec.i), digits), signif(var(vec.i), digits))
        vec_text[[i]]$stat.type <- c("median", "variance")
    }
  }
  banocc::cat_v("Done.\n", verbose)
  banocc::cat_v("End get_vec_text\n", verbose, num_level=num_level)
  return(vec_text)
}

#' @rdname get_text_to_plot
get_gamma_text <-
function(alpha_vec, beta_vec, digits=2, verbose=FALSE, num_level=0){

  banocc::cat_v("Begin get_gamma_text\n", verbose, num_level=num_level)
  p <- length(alpha_vec)
  
  banocc::cat_v("Getting gamma parameters...", verbose, num_level=num_level+1)
  gamma_text <- vector("list", length=p)
  gamma.means <- alpha_vec/beta_vec
  gamma.vars  <- alpha_vec/(beta_vec^2)
  for(i in 1:p){
    max.x <- gamma.means[i] + 4 * sqrt(gamma.vars[i])
    gamma_text[[i]] <- data.frame(i=i,
                                  stat=c(signif(gamma.means[i], digits),
                                         signif(gamma.vars[i], digits)),
                                  stat.type=c("mean", "variance"), x=0.8 * max.x)
  }
  banocc::cat_v("Done.\n", verbose)
  banocc::cat_v("End get_gamma_text\n", verbose, num_level=num_level)
  return(gamma_text)
}
