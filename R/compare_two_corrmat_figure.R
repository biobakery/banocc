#' Generate a heatmap figure comparing two correlation matrices
#'
#' @inheritParams melt_two_corr_mat
#' @inheritParams get_two_corrmat_df
#' @inheritParams make_two_corrmat_figure
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' 
#' @examples
#' Rho.lower <- matrix(c(1, 0, 0, 0.5, 0,
#'                       0, 1, 0, 0, 0,
#'                       0, 0, 1, 0, 0,
#'                       0.5, 0, 0, 1, -0.5,
#'                       0, 0, 0, -0.5, 1), ncol=5)
#' X <- rmvnorm(5, rep(0, 5), Rho.lower)
#' Rho.upper <- cov2cor(X %*% t(X))
#' truestatus_mat <- as.logical(Rho.lower != 0 - diag(5))
#' title <- "Estimating based on 5 samples"
#' categories <- c("Truth", "Estimate (5 samples)")
#' compare_two_corrmat_figure(Rho.lower, Rho.upper, title, truestatus_mat,
#'                            categories=categories)
#' compare_two_corrmat_figure(Rho.lower, Rho.upper, title, truestatus_mat,
#'                            categories=categories, legend.position="vertical")
#'

compare_two_corrmat_figure <-
  function(Rho.lower, Rho.upper, title=NA, truestatus_mat=NA, verbose=FALSE, 
           categories=c("Truth", "Estimate"),
           corr_label=expression(rho["jj\'"]),
           include_text=TRUE, text.size=20, corr.text.size=6,
           legend.position="horizontal", num_level=0)
{
    banocc::cat_v("Begin compare_two_corrmat_figure.\n", verbose, num_level)

    Rho.df <- banocc::get_two_corrmat_df(Rho.lower=Rho.lower,
                                         Rho.upper=Rho.upper,
                                         truestatus_mat=truestatus_mat,
                                         verbose=verbose,
                                         num_level=num_level + 1)
    Rho.plot <- banocc::make_two_corrmat_figure(Rho.df, title=title,
                                                verbose=verbose,
                                                categories=categories,
                                                corr_label=corr_label,
                                                include_text=include_text,
                                                text.size=text.size,
                                                corr.text.size=corr.text.size,
                                                legend.position=legend.position,
                                                num_level=num_level + 1)
    
    banocc::cat_v("Printing plot...", verbose, num_level + 1)
    print(Rho.plot)
    banocc::cat_v("Done.\n", verbose)
    banocc::cat_v("End compare_two_corrmat_figure.\n", verbose, num_level)
}

#' Make a heatmap comparing two matrices when they are in data frame format.
#'
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @param Rho.df The data frame; each row is an (i,j) element of one of the two
#'   correlation matrices.  It must have the following columns:
#'   \itemize{
#'     \item \code{Var1, y}: the (i,j) indices as factors
#'     \item \code{cor.value}: the value of the correlation of the (i,j)th
#'       element
#'   }
#'   It can also have the optional column \code{truestatus}, which is TRUE if
#'   i and j are supposed to be associated and NA if not.
#' @param title The title of the plot
#' @param categories A character vector of length two. The first element
#'   indicates the label for the lower triangle; the second indicates the label
#'   for the upper triangle.
#' @param corr_label The label for the legend (can be an expression or a
#'   string).
#' @param include_text Boolean: should the values of the correlation matrix be
#'   added as text to the plot?
#' @param text.size The base size of the text for \code{ggplot} (in pt size).
#' @param legend.position One of either \code{"horizontal"} or
#'   \code{"vertical"}; indicates whether the legends will be positioned at the
#'   bottom (\code{"horizontal"}) or left (\code{"vertical"}) of the plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 guide_legend
#' @importFrom stringr str_c
#' 
make_two_corrmat_figure <-
  function(Rho.df, title=NA, verbose=FALSE, 
           categories=c("Truth", "Estimate"),
           corr_label=expression(rho["jj\'"]),
           include_text=TRUE, text.size=20,
           corr.text.size=6,
           legend.position="horizontal", num_level=0)
{
    banocc::cat_v("Begin make_two_corrmat_figure.\n", verbose, num_level)
    if (length(categories) !=2){
      stop("Must specify EXACTLY two categories.")
    }
    if (!(legend.position %in% c("horizontal", "vertical"))){
      stop("legend.position must be one of \"horizontal\" or \"vertical\".")
    }

    if("truestatus" %in% names(Rho.df)){      
      banocc::cat_v("Generating basic plot...", verbose, num_level + 1)
      p_Rho <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = Rho.df, 
                           ggplot2::aes(x=Var1, y=y, fill=cor.value,
                                        colour=truestatus, size=truestatus)) +
        ggplot2::scale_size_manual(
            values=c("TRUE"=2), guide=ggplot2::guide_legend(title="True Assn",
                                    title.position="top")) +
        ggplot2::scale_colour_manual(
            values=c("TRUE" = rgb(1, 1, 0)),
            guide=ggplot2::guide_legend(title="True Assn",
                                        title.position="top"))
      banocc::cat_v("Done.\n", verbose)
    } else {
      banocc::cat_v("Generating basic plot...", verbose, num_level + 1)
      p_Rho <- ggplot2::ggplot() + 
        ggplot2::geom_tile(data = Rho.df, ggplot2::aes(x=Var1, y=y,
                                                       fill=cor.value))
      banocc::cat_v("Done.\n", verbose)
    }
    
    banocc::cat_v("Adding legend...", verbose, num_level + 1)
    triangle.legend <- paste0(c("Upper ", "Lower "), "Triangle = ",
                              rev(categories))
    p_Rho <- p_Rho + ggplot2::coord_fixed() +
      ggplot2::scale_x_discrete(name=stringr::str_c(triangle.legend,
                                                    collapse="\n"),
                                breaks=NULL) +
      ggplot2::scale_y_discrete(name="", breaks=NULL)
      ## ggplot2::geom_rect(ggplot2::aes(xmin=0.4, xmax=p+0.6, ymin=-0.6,
      ##                                 ymax=0.4),
      ##                    fill="white") +
      ## ggplot2::annotate("text", label=str_c(triangle.legend, collapse="\n"),
      ##                   x=p/2 + 0.5, y=-0.1, colour="black", 
      ##                   size=(5/14) * (text.size - 4))
    if (legend.position=="horizontal"){
        colourbar <- ggplot2::guide_colourbar(
            title=corr_label, title.position="top",
            label.theme=ggplot2::element_text(angle=270, size=text.size - 4))
        p_Rho <- p_Rho  +
            ggplot2::theme(panel.background=ggplot2::element_rect(fill="black",
                                                         colour="white"),
                           panel.border=ggplot2::element_blank(),
                           plot.background=ggplot2::element_blank(),
                           panel.grid = ggplot2::element_blank(),
                           text=ggplot2::element_text(size=text.size),
                           legend.title=ggplot2::element_text(hjust=-1),
                           legend.title.align=0.5,
                           legend.position="bottom",
                           legend.box="horizontal",
                           axis.title.y=ggplot2::element_blank(),
                           title=ggplot2::element_text(vjust=1.2),
                           strip.background=ggplot2::element_rect(fill="white",
                               colour="black"))
    } else if (legend.position=="vertical"){
        colourbar <- ggplot2::guide_colourbar(title=corr_label,
                                              title.position="top")
        p_Rho <- p_Rho +
            ggplot2::theme(panel.background=ggplot2::element_rect(fill="black",
                                                            colour="white"),
                           panel.border=ggplot2::element_blank(),
                           plot.background=ggplot2::element_blank(),
                           panel.grid = ggplot2::element_blank(),
                           text=ggplot2::element_text(size=text.size),
                           legend.title=ggplot2::element_text(hjust=-1),
                           legend.title.align=0.5,
                           axis.title.y=ggplot2::element_blank(),
                           title=ggplot2::element_text(vjust=1.2),
                           strip.background=ggplot2::element_rect(fill="white",
                               colour="black"))
    }
    p_Rho <- p_Rho +
        ggplot2::scale_fill_gradient2(midpoint=0,
                                      mid = rgb(0.8,0.8,0.8),
                                      low="blue", high="red",
                                      limits=c(-1.01, 1.01), na.value="black",
                                      guide=colourbar, space="Lab")
      

    if(length(title)==1 && !is.na(title)){
      p_Rho <- p_Rho + ggplot2::ggtitle(title)
    }
    banocc::cat_v("Done.\n", verbose)
    
    if (include_text){
        banocc::cat_v("Adding text...", verbose, num_level + 1)
        p_Rho <- p_Rho + 
          ggplot2::geom_text(colour="black", size=corr.text.size, data=Rho.df,
                             ggplot2::aes(x=Var1, y=y,
                                          label=round(cor.value, 3)))
        banocc::cat_v("Done.\n", verbose)
    }
    
    return(p_Rho)
    banocc::cat_v("End make_two_corrmat_figure.\n", verbose, num_level)
}

#' Get a data frame for use with \code{make_two_corrmat_figure} from two
#'   matrices
#' @inheritParams melt_two_corr_mat
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @param truestatus_mat An optional matrix indicating by Boolean values which
#'   elements of Rho.lower and Rho.upper should be non-zero.  The dimensions
#'   must match those of Rho.lower and Rho.upper.
get_two_corrmat_df <-
    function(Rho.lower, Rho.upper, truestatus_mat=NA, verbose=FALSE,
             num_level=0)
{
    banocc::cat_v("Begin get_two_corrmat_df.\n", verbose, num_level)
    if (any(dim(Rho.lower)!=dim(Rho.upper))){
      stop("Dimensions of Rho.lower and Rho.upper must match.")
    }

    p <- ncol(Rho.lower)
    Rho.melt <- banocc::melt_two_corr_mat(Rho.lower, Rho.upper, verbose=verbose, num_level=num_level + 1)
    rownames(Rho.lower) <- paste0("f", seq_len(nrow(Rho.lower)))
    idx.mat <- sapply(rownames(Rho.lower),
                      function(a) paste(a, rownames(Rho.lower)))
    
    if(is.matrix(truestatus_mat) || is.data.frame(truestatus_mat)){
      banocc::cat_v("Adding true status...", verbose, num_level + 1)
      truestatus_melt <- banocc::get_melt_dataset(truestatus_mat, cor_mat=TRUE)
      
      diag.elt        <- which(paste(truestatus_melt$Var1, 
                                     truestatus_melt$Var2) %in% diag(idx.mat))
      truestatus_melt <- truestatus_melt[-diag.elt,]
      
      truestatus_melt$value[!truestatus_melt$value] <- NA
      names(truestatus_melt) <- c("Var1", "Var2", "truestatus")
      Rho.melt <- banocc::mymerge(Rho.melt, truestatus_melt)
      banocc::cat_v("Done.\n", verbose)
    }
    banocc::cat_v("End get_two_corrmat_df.\n", verbose, num_level)
    return(Rho.melt)
}
