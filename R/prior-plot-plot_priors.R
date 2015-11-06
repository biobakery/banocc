#' Generate prior plots for several different possible types of priors
#'
#' @param Sigma.sample A list of variance-covariance matrices
#' @param mu.sample A list of vectors
#' @param sigma.sample A list of vectors; either \code{sigma_sample} or BOTH
#'   \code{alpha_vec} and \code{beta_vec} must be provided.
#' @param R.sample A list of correlation matrices
#' @param add.stats Boolean: add text of summary statistics to the plots?
#' @inheritParams get_vec_text
#' @inheritParams get_vec_plot
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @return A named list: \code{plot} is a \code{ggplot} object; \code{data} is
#'   the data frame which was used to generate the plot
#' @name plot_priors

#' @rdname plot_priors
plot_Sigma_prior <- function(Sigma.sample, add.stats=FALSE, digits=2,
                             center.stat="mean", verbose=FALSE, num_level=0){

      banocc::cat_v("Begin plot_Sigma_prior\n", verbose, num_level=num_level)
      Sigma.to_plot <- banocc::get_Matrix_to_plot(Matrix.sample=Sigma.sample,
                                                  diag=TRUE, verbose=verbose,
                                                  num_level=num_level+1)
      if (add.stats){
          Sigma_text    <- banocc::get_Matrix_text(Matrix.sample=Sigma.sample,
                                                   diag=TRUE, digits=digits,
                                                   center.stat=center.stat,
                                                   verbose=verbose,
                                                   num_level=num_level+1)
      } else {
          Sigma_text <- NULL
      }

      Sigma_plot <- banocc::get_Matrix_plot(Matrix_to_plot=Sigma.to_plot,
                                            Matrix_text=Sigma_text,
                                            verbose=verbose,
                                            num_level=num_level+1)
      
      banocc::cat_v("End plot_Sigma_prior\n", verbose, num_level=num_level)
      return(list(plot=Sigma_plot, data=Sigma.to_plot))
  }

#' @rdname plot_priors
plot_mu_prior <-
function(mu.sample, add.stats=FALSE, digits=2, center.stat="mean",
         verbose=FALSE, num_level=0){
  banocc::cat_v("Begin plot_mu_prior\n", verbose, num_level=num_level)

  mu.to_plot <- banocc::get_vec_to_plot(vec.sample=mu.sample, verbose=verbose,
                                        num_level=num_level+1)
  if (add.stats){
    mu_text <- banocc::get_vec_text(vec.sample=mu.sample, digits=digits,
                                    center.stat=center.stat, verbose=verbose,
                                    num_level=num_level+1)
  } else {
    mu_text <- NULL
  }

  mu_plot <- banocc::get_vec_plot(vec_to_plot=mu.to_plot, vec_text=mu_text,
                                  verbose=verbose, num_level=num_level+1)

  banocc::cat_v("End plot_mu_prior\n", verbose, num_level=num_level)
  return(list(plot=mu_plot, data=mu.to_plot))
}

#' @rdname plot_priors
plot_sigma_prior <-
function(sigma.sample=NULL, alpha_vec=NULL, beta_vec=NULL, mat_plot=TRUE,
         add.stats=FALSE, digits=2, center.stat="mean", verbose=FALSE,
         num_level=0){

  banocc::cat_v("Begin plot_sigma_prior\n", verbose, num_level=num_level)
  sigma_text <- NULL
  if(is.null(sigma.sample)){
    if(is.null(alpha_vec) || is.null(beta_vec)){
      stop("Must provide sigma.sample OR alpha_vec and beta_vec")
    }
    sampled <- FALSE
    sigma.to_plot <- banocc::get_gamma_to_plot(alpha_vec=alpha_vec,
                                               beta_vec=beta_vec,
                                               verbose=verbose,
                                               num_level=num_level+1)
    if (add.stats){
      sigma_text <- banocc::get_gamma_text(alpha_vec=alpha_vec,
                                           beta_vec=beta_vec, digits=digits,
                                           verbose=verbose,
                                           num_level=num_level+1)
    }
  } else {
    sampled <- TRUE
    sigma.to_plot <- banocc::get_vec_to_plot(vec.sample=sigma.sample,
                                             verbose=verbose, num_level=num_level)
    if (add.stats){
      sigma_text <- banocc::get_vec_text(vec.sample=sigma.sample, digits=digits,
                                         center.stat=center.stat, verbose=verbose,
                                         num_level=num_level+1)
    }
  }
  
  sigma_plot <- banocc::get_vec_plot(vec_to_plot=sigma.to_plot, vec_text=sigma_text,
                                     sampled=sampled, mat_plot=mat_plot,
                                     verbose=verbose, num_level=num_level+1)

  banocc::cat_v("End plot_sigma_prior\n", verbose, num_level=num_level)
  return(list(plot=sigma_plot, data=sigma.to_plot))
}

#' @rdname plot_priors
plot_R_prior <-
function(R.sample, add.stats=FALSE, digits=2, center.stat="mean", verbose=FALSE,
         num_level=0){
  banocc::cat_v("Begin plot_R_prior\n", verbose, num_level=num_level)
  R.to_plot <- banocc::get_Matrix_to_plot(Matrix.sample=R.sample, diag=FALSE,
                                          verbose=verbose, num_level=num_level+1)
  if (add.stats){
    R_text <- banocc::get_Matrix_text(Matrix.sample=R.sample, diag=FALSE,
                                      digits=digits, center.stat=center.stat,
                                      verbose=verbose, num_level=num_level+1)
  } else {
    R_text <- NULL
  }
  
  R_plot <- banocc::get_Matrix_plot(Matrix_to_plot=R.to_plot,
                                    Matrix_text=R_text, xlim=c(-1, 1),
                                    verbose=verbose, num_level=num_level+1)
  banocc::cat_v("End plot_R_prior\n", verbose, num_level=num_level)
  return(list(plot=R_plot, data=R.to_plot))  
}
