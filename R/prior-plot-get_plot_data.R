#' Get either a sample of vectors or matrices formatted for plotting
#'
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @param vec.sample A list of vectors
#' @param Matrix.sample A list of matrices
#' @param diag Boolean: include the diagonal elements of the matrix samples?
#' @return A data frame of indices and sample values (in the case of
#'   \code{get_gamma_to_plot}, pdf values for a range are returned)
#' @name get_plot_data
#'
#' 

#' @rdname get_plot_data
get_vec_to_plot <-
function(vec.sample, verbose=FALSE, num_level=0){

  banocc::cat_v("Begin get_vec_to_plot\n", verbose, num_level=num_level)
  p <- length(vec.sample[[1]])
  
  banocc::cat_v("Getting vector parameter df...", verbose,
                num_level=num_level+1)
  vec.to.plot <- vector("list", length=p)
  for (i in 1:p){
    vec.to.plot[[i]] <- data.frame(i=i, x=unlist(lapply(vec.sample, '[', i)))
  }
  
  banocc::cat_v("Done.\n", verbose)
  vec.to.plot <- Reduce(rbind, vec.to.plot)
  banocc::cat_v("End get_vec_to_plot\n", verbose, num_level=num_level)
  
  return(vec.to.plot)
}

#' @rdname get_plot_data
get_Matrix_to_plot <-
function(Matrix.sample, diag=TRUE, verbose=FALSE, num_level=0){

  banocc::cat_v("Begin get_Matrix_to_plot\n", verbose, num_level=num_level)
  p <- nrow(Matrix.sample[[1]])
  
  banocc::cat_v("Getting matrix parameter df...", verbose,
                num_level=num_level+1)
  Matrix.to.plot <- vector("list", length=(choose(p, 2) * 2) + p * diag)
  idx <- 0
  next.indices <- 1:p
  
  for (i in 1:p){
    Matrix.i <- Reduce(rbind, lapply(Matrix.sample, '[', i, ))
    if(!diag){
      next.indices <- setdiff(1:p, i)
    }
    for (k in next.indices){
      idx <- idx + 1
      val <- Matrix.i[, k]
      Matrix.to.plot[[idx]] <-data.frame(i=i, k=k, x=val)
    }
  }
  banocc::cat_v("Done.\n", verbose)
  Matrix.to.plot <- Reduce(rbind, Matrix.to.plot)

  banocc::cat_v("End get_Matrix_to_plot\n", verbose, num_level=num_level)
  return(Matrix.to.plot)
}

#' @rdname get_plot_data
get_gamma_to_plot <-
function(alpha_vec, beta_vec, verbose=FALSE, num_level=0){

  banocc::cat_v("Begin get_gamma_to_plot\n", verbose, num_level=num_level)
  
  banocc::cat_v("Getting gamma distribution df...", verbose,
                num_level=num_level+1)
  d <- length(alpha_vec)
  gamma.to.plot <- vector("list", length=d)
  for (i in 1:d){
    a.i <- alpha_vec[i]
    b.i <- beta_vec[i]
    avg.i <- a.i/b.i
    var.i <- a.i/(b.i^2)
    min.x <- max(avg.i - 4 * sqrt(var.i), 0)
    max.x <- avg.i + 4 * sqrt(var.i)
    x <- seq(min.x, max.x, (max.x - min.x)/1000)
    y <- dgamma(x, shape = a.i, rate = b.i)
    gamma.to.plot[[i]] <- data.frame(i=i, k=i, x=x, y=y)
  }
  banocc::cat_v("Done.\n", verbose)
  gamma.to.plot <- Reduce(rbind, gamma.to.plot)
  banocc::cat_v("End get_gamma_to_plot\n", verbose, num_level=num_level)
  return(gamma.to.plot)
}
