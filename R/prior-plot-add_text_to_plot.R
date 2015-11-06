#' Add summary statistics to a plot
#'
#' @param Matrix.plot A \code{ggplot} plot faceted with \code{facet_grid(i ~k)}
#' @param vec.plot A \code{ggplot} plot faceted with \code{facet_wrap(~ i)}
#' @param p The dimension of the matrix
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @param Matrix_text A list of summary statistics for each parameter index. This
#'   should be in the same format as that returned by
#'   \code{\link{get_Matrix_text}} 
#' @param vec_text A list of summary statistics for each
#'   parameter index. This should be in the same format as that returned by
#'   \code{\link{get_vec_text}}.
#' @return Returns a \code{ggplot} plot with the text added
#' @name add_text_to_plot
#'
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 geom_text

#' @rdname add_text_to_plot
add_matrix_text <-
function(Matrix.plot, Matrix_text, p, verbose=FALSE, num_level=0){
  banocc::cat_v("Begin add_matrix_text\n", verbose, num_level=num_level)
  banocc::cat_v("Building plot...", verbose, num_level=num_level+1)
  ggb <- ggplot2::ggplot_build(Matrix.plot)
  banocc::cat_v("simplifying output...", verbose)
  ggb.data <- ggb$data[[1]]

  banocc::cat_v("Getting panels...", verbose)
  panels <- matrix(seq(1, p^2), ncol=p)
  used.panels <- as.numeric(unique(ggb.data$PANEL))
  banocc::cat_v("Getting text per panel...", verbose)
  Matrix_text <- lapply(seq_along(used.panels), function(i){
    panel <- used.panels[i]
    row.idx <- panel %% p
    if(!row.idx){
      row.idx <- p
    }
    plot.subset <- subset(ggb.data, 
                          ggb.data$group %in% panels[row.idx,]
    )
    
    Matrix_text[[i]]$y <- max(plot.subset$y) * c(0.9, 0.8)
    return(Matrix_text[[i]])
  })
  banocc::cat_v("Concatenating text...", verbose)
  Matrix_text <- Reduce(rbind, Matrix_text)
  banocc::cat_v("Plotting text...", verbose)
  Matrix.plot <- Matrix.plot +
    ggplot2::geom_text(aes(label=stat, x=x, y=y, colour=stat.type),
                       data=Matrix_text, hjust=0, vjust=1)
  banocc::cat_v("End add_matrix_text", verbose, num_level=num_level)
  return(Matrix.plot)
}

#' @rdname add_text_to_plot
add_vec_text <-
function(vec.plot, vec_text, verbose=FALSE, num_level=0){
  banocc::cat_v("Begin add_vec_text\n", verbose, num_level=num_level)
  ggb <- ggplot2::ggplot_build(vec.plot)
  ggb.data <- ggb$data[[1]]
  
  panels <- unique(ggb.data$group)
  vec_text <- lapply(panels, function(i){
    plot.subset <- subset(ggb.data, ggb.data$group==i)
    vec_text[[i]]$y <- max(plot.subset$y) * c(0.9, 0.8)
    return(vec_text[[i]])
  })
  vec_text    <- Reduce(rbind, vec_text)
  vec.plot <- vec.plot + 
    ggplot2::geom_text(aes(label=stat, x=x, y=y, colour=stat.type),
                       data=vec_text, hjust=0, vjust=1)
  banocc::cat_v("End add_vec_text\n", verbose, num_level=num_level)
  return(vec.plot)
}
