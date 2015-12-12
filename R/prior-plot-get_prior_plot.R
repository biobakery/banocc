#' Get density plots for either matrix or vector parameters
#'
#' @param Matrix_to_plot A data frame with columns \code{i} and \code{k}
#'   indicating index and column \code{x} indicating a value of the parameter for
#'   that index.
#' @param vec_to_plot A data frame with columns \code{i} indicating index and
#'   column \code{x} indicating a value of the parameter for that index. May
#'   optionally contain column \code{y}, in which case \code{y} is a pdf
#'   evaluated at \code{x} (if \code{sampled=FALSE})
#' @param xlim Optional limits on the domain values of plot
#' @inheritParams add_vec_text
#' @inheritParams cat_v
#' @param mat_plot Boolean: Make a matrix of plots with \code{facet_grid} or
#'   wrap the plots with \code{facet_wrap}?
#' @param sampled Boolean: is this an estimated density based on sample, or an
#'   actual calculated density? (if \code{sampled=FALSE}, then \code{vec_to_plot}
#'   must contain columns \code{x} and \code{y})
#' @inheritParams sample_mu_prior
#' @return Returns a \code{ggplot} object which is the final plot of the samples.
#' @name get_prior_plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap

#' @rdname get_prior_plot
get_Matrix_plot <-
function(Matrix_to_plot, xlim = NULL, Matrix_text=NULL, verbose=FALSE,
         num_level=0){

  banocc::cat_v("Begin get_Matrix_plot\n", verbose, num_level=num_level)
  add_stats <- !is.null(Matrix_text)
  p <- max(Matrix_to_plot$i)

#  if (verbose) print(head(Matrix_to_plot))
  
  banocc::cat_v("Generating plot...", verbose, num_level=num_level+1)
  Matrix_plot <- ggplot2::ggplot(aes(x=x), data=Matrix_to_plot) + 
    ggplot2::geom_density(alpha=0.2, colour="black", fill="red") +
        ggplot2::facet_grid(k~i,scales="free") ## +
  ##ggplot2::geom_histogram(aes(y=..density..), binwidth=hist.binwidth)
  if (!is.null(xlim)) Matrix_plot <- Matrix_plot + ggplot2::xlim(xlim[1], xlim[2])
  banocc::cat_v("Done.\n", verbose)
  
  if (add_stats){
    Matrix_plot <- banocc::add_matrix_text(Matrix.plot=Matrix_plot,
                                           Matrix_text=Matrix_text, p=p,
                                           verbose=verbose,
                                           num_level=num_level+1)
  }
  banocc::cat_v("End get_Matrix_plot\n", verbose, num_level=num_level)
  return(Matrix_plot)
}

#' @rdname get_prior_plot
get_vec_plot <-
function(vec_to_plot, vec_text=NULL, sampled=TRUE, mat_plot=FALSE,
         verbose=FALSE, num_level=0){

  banocc::cat_v("Begin get_vec_plot\n", verbose, num_level=num_level)
  add_stats <- !is.null(vec_text)
  
  banocc::cat_v("Generating plot...", verbose, num_level=num_level+1)
  if (mat_plot) vec_to_plot$k <- vec_to_plot$i
  
  if(sampled){
    vec_plot <- ggplot2::ggplot(aes(x=x), data=vec_to_plot) + 
      ggplot2::geom_density(alpha=0.2, colour="black", fill="red")
  } else {
    vec_plot <- ggplot2::ggplot(aes(x=x, y=y), data=vec_to_plot) + 
      ggplot::geom_line(fill="red")
  }
  
  if (mat_plot){
    vec_plot <- vec_plot + ggplot2::facet_grid(k ~ i, scales="free")
    if (add_stats){
      vec_text <- lapply(vec_text, function(t){ t$k <- t$i; return(t)})
    }
  } else {
    vec_plot <- vec_plot + ggplot2::facet_wrap(~ i, scales="free")
  }
  banocc::cat_v("Done.\n", verbose)
  
  if (add_stats){
    if (mat_plot){
      vec_plot <- banocc::add_matrix_text(vec_plot, vec_text,
                                          p=max(vec_to_plot$i),
                                          verbose=verbose,
                                          num_level=num_level+1)
    } else{
      vec_plot <- banocc::add_vec_text(vec_plot, vec_text, verbose=verbose,
                                       num_level=num_level+1)
    }
  }
  banocc::cat_v("End get_vec_plot\n", verbose, num_level=num_level)
  return(vec_plot)
}
