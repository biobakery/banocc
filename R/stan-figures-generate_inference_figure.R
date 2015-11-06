#' Generate a figure of the estimates and credible intervals of a
#'   correlation matrix
#'
#' Only \code{generate_inference_figure_MCMC} can generate credible intervals
#'   of arbitrary confidence; \code{generate_inference_figure_optim} only
#'   produces 95\% credible intervals.
#'
#' @inheritParams get_inference_MCMC_df
#' @inheritParams rstan::extract
#' @inheritParams get_credible_intervals
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @inheritParams make_samples_list
#' @inheritParams corr_inference_figure
#' @inheritParams generate_estimate_figure_optim
#' @inheritParams rstan::unconstrain_pars
#' @inheritParams get_gradient

#' @name generate_inference_figure
#'
#' @importFrom rstan unconstrain_pars
#' @importFrom rstan extract

#' @rdname generate_inference_figure
generate_inference_figure_MCMC <-
function(stanfit, Rho.true, type="marginal.hpd", conf=0.95, ylimits=NULL,
         ln.Rho=TRUE, thin=1, Rho.comp=NULL, segment_size=2, text.size=18,
         axis.ticks.x=TRUE, axis.text.size.x=NULL,
         verbose=FALSE, num_level=0){
    banocc::cat_v("Begin generate_inference_figure_MCMC\n", verbose, num_level)
    CI.df <- banocc::get_inference_MCMC_df(stanfit, Rho.true, type=type,
                                           conf=conf,
                                   Sigma.comp=Rho.comp,
                                   ln.Rho=ln.Rho, thin=thin,
                                   verbose=verbose, num_level=num_level + 1)
    CI.plot <- banocc::corr_inference_figure(CI.df, type=type, conf=conf,
                                             ylimits=ylimits,
                                             segment_size=segment_size,
                                             text.size=text.size,
                                             axis.ticks.x=axis.ticks.x,
                                             axis.text.size.x=axis.text.size.x,
                                             verbose=verbose,
                                             num_level=num_level + 1)

    banocc::cat_v("Printing plot...", verbose, num_level + 1)
    print(CI.plot)
    banocc::cat_v("Done.\n", verbose)
    
    banocc::cat_v("End generate_inference_figure_MCMC\n", verbose,)
}

#' @rdname generate_inference_figure
generate_inference_figure_optim <-
function(Fit, stanfit, fn, Rho.true,
         epsilon=1e-6, ylimits=NULL, text.size=20, verbose=FALSE,
         num_level=0){

    banocc::cat_v("Begin generate_inference_figure_optim\n", verbose, num_level)
    banocc::cat_v("Getting estimates...", verbose, num_level + 1)
    R.est.elt <- Fit$par[grep("ln_Rho", names(Fit$par))]
    p         <- sqrt(length(R.est.elt))
    R.est     <- matrix(R.est.elt, ncol=p)
    banocc::cat_v("Done.\n", verbose)

    banocc::cat_v("Getting unconstrained parameters...", verbose, num_level + 1)
    pars <- list(
        mu = Fit$par[grep("^mu", names(Fit$par))],
        L  = matrix(Fit$par[grep("^L", names(Fit$par))], ncol=p),
        sigma = Fit$par[grep("^sigma", names(Fit$par))]
        )
    banocc::cat_v("Done.\n", verbose)

    banocc::cat_v("Getting transformed variance...", verbose, num_level + 1)
    param.u <- rstan::unconstrain_pars(stanfit, pars)
    trans.var <- banocc::get_trans_variance(Fit$hessian, fn, param.u=param.u,
                                            stanfit, epsilon)
    R.var.elt <- trans.var[grep("ln_Rho", colnames(trans.var)),
                           grep("ln_Rho", rownames(trans.var))]
    R.var <- 0 * diag(p)
    R.var[lower.tri(R.var)] <- diag(R.var.elt)
    R.var[upper.tri(R.var)] <- t(R.var)[upper.tri(R.var)]
    banocc::cat_v("Done.\n", verbose)

    banocc::cat_v("Getting approximate 95% CI...", verbose, num_level + 1)
    R.ci <- list(
        lower = R.est - qnorm(0.95) * R.var,
        upper = R.est + qnorm(0.95) * R.var
        )
    R.ci$lower[which(R.ci$lower < -1)] <- -1
    R.ci$lower[which(R.ci$lower > 1)] <- 1
    R.ci$upper[which(R.ci$upper < -1)] <- -1
    R.ci$upper[which(R.ci$upper > 1)] <- 1
    banocc::cat_v("Done.\n", verbose)

    inference.df <- banocc::get_inference_df(Rho.true, R.est, R.ci,
                                             verbose=verbose,
                                             num_level=num_level + 1)
    R.plot <- banocc::corr_inference_figure(inference.df, ylimits=ylimits,
                                            text.size=text.size,
                                            verbose=verbose,
                                            num_level=num_level + 1)
    banocc::cat_v("Printing plot...", verbose, num_level + 1)
    print(R.plot)
    banocc::cat_v("Done.\n", verbose)
    
    banocc::cat_v("End generate_inference_figure_optim\n", verbose, num_level)
}


#' Get the data frames used for plotting MCMC inference figure
#' @inheritParams rstan::extract
#' @inheritParams get_credible_intervals
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @inheritParams make_samples_list
#' @param Rho.true The true correlation matrix
#' @param Rho.comp An optional correlation matrix in the composition
#' @param stanfit An object of class \code{stanfit}.


get_inference_MCMC_df <- function(stanfit, Rho.true, type="marginal.hpd",
                                  conf=0.95, Rho.comp=NULL,
                                  ln.Rho=TRUE, thin=1, verbose=FALSE,
                                  num_level=0){
    banocc::cat_v("Begin get_inference_MCMC_df\n", verbose, num_level)
    banocc::cat_v("Begin getting posterior samples\n", verbose, num_level + 1)
    posterior.samples <- rstan::extract(stanfit,
                                        permuted= FALSE)
    post.samples.list <- banocc::make_samples_list(posterior.samples, thin=thin,
                                                   concatenate.chains=TRUE,
                                                   verbose)
    banocc::cat_v("End getting posterior samples.\n", verbose, num_level + 1)

    banocc::cat_v("Getting median estimate...", verbose, num_level + 1)
    if (ln.Rho){
        estimates.median <-
            banocc::get_posterior_estimates(post.samples.list,
                                            estimate_method="median",
                                            parameter.names="ln_Rho")

        R.est.median <- estimates.median$ln_Rho
    } else {
        estimates.median <-
            banocc::get_posterior_estimates(post.samples.list,
                                            estimate_method="median",
                                            parameter.names="Rho")

        R.est.median <- estimates.median$Rho

    }
    banocc::cat_v("Done.\n", verbose)
    
    if (ln.Rho){
        credible.intervals <-
            banocc::get_credible_intervals(post.samples.list, type=type,
                                           parameter.names="ln_Rho", conf=conf,
                                           list=TRUE)

        R.ci <- credible.intervals$ln_Rho
    } else {
        credible.intervals <-
            banocc::get_credible_intervals(post.samples.list, type=type,
                                           parameter.names="Rho", conf=conf,
                                           list=TRUE)

        R.ci <- credible.intervals$Rho
    }
    
    names(R.ci) <- c("lower", "upper")
    R.CI.df <- banocc::get_inference_df(Rho.true, R.est.median, R.ci,
                                        Rho.comp, corr=TRUE, verbose,
                                        num_level=num_level + 1)
    banocc::cat_v("End get_inference_MCMC_df\n", verbose, num_level)
    return(R.CI.df)
}

#' Get the data frame for the inference plot from the correlation or covariance
#'   matrices
#' @param Sigma.true The true value of the correlation or variance/covariance matrix
#' @param Sigma.median The estimated value of the correlation or variance/covariance matrix (assumed to be the posterior median)
#' @param Sigma.ci A list with two values: \code{lower} is the lower pointwise
#'   limits of the 95\% credible interavls and \code{upper} is the upper
#'   pointwise limits of the 95\% credible intervals.  Both elements are
#'   matrices.
#' @param Sigma.comp An optional matrix of observed correlation values in the
#'   composition
#' @param corr Boolean: Are these correlation matrices?
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#'
#' @importFrom reshape2 melt
#' 
get_inference_df <- function(Sigma.true, Sigma.median, Sigma.ci, Sigma.comp=NULL,
                             corr=TRUE, verbose=FALSE, num_level=0){
  banocc::cat_v("Begin get_inference_df\n", verbose, num_level)
      
  Sigma.ci$median <- Sigma.median
  Sigma.ci$true   <- Sigma.true
  if(!is.null(Sigma.comp)) Sigma.ci$comp <- Sigma.comp
  
  ## Remove diagonal if corr=TRUE; keep if corr=FALSE
  get.diag <- corr

  banocc::cat_v("Getting appropriate elements...", verbose, num_level + 1)
  Sigma.ci.tri <- lapply(Sigma.ci, function(mat){
      mat[lower.tri(mat, diag=get.diag)] <- NA
      return(mat)
  })
  banocc::cat_v("Done.\n", verbose)

  banocc::cat_v("Getting data frame to plot...", verbose, num_level + 1)
  Sigma.ci.melt <- mapply(function(S, val){
      na.omit(reshape2::melt(S, value.name=val))
      },
                          Sigma.ci.tri, names(Sigma.ci.tri), SIMPLIFY=FALSE)
  
  Sigma.ci.df          <- Reduce(banocc::mymerge, Sigma.ci.melt)
  Sigma.ci.df$x        <- paste0("(", Sigma.ci.df$Var1, ", ", 
                                 Sigma.ci.df$Var2, ")")
  Sigma.ci.df$var.type <- apply(Sigma.ci.df, 1, function(row){
      ifelse(row[1]==row[2], "variance", "covariance")
  })
  banocc::cat_v("Done.\n", verbose)


  banocc::cat_v("Get rectangle boundaries...", verbose, num_level + 1)
  ## The minimum and maximum x values for the emphasis rectangles
  Sigma.ci.df$xmin <- seq_len(nrow(Sigma.ci.df)) - 0.5
  Sigma.ci.df$xmax <- Sigma.ci.df$xmin + 1
  Sigma.ci.df$xmin.t <- Sigma.ci.df$xmin + 0.25
  Sigma.ci.df$xmax.t <- Sigma.ci.df$xmax - 0.25
  banocc::cat_v("Done.\n", verbose)
  
    banocc::cat_v("End get_inference_df\n", verbose, num_level)
  return(Sigma.ci.df)
}

#' Generates an interval plot for correlation/covariance matrix elements
#'
#' @param title A title for the plot
#' @param ylimits A vector of length two giving the lower and upper limits of
#'   the y-axis
#' @param segment_size The size of the segments drawn in the plot
#' @param inference.df The data frame to use in plotting. Must have the
#'   following columns:
#'   \itemize{
#'     \item \code{median}: The median posterior value for each element
#'     \item \code{true}: The true value for each element
#'     \item \code{lower},\code{upper}: The upper and lower credible interval
#'       limits for each element
#'     \item \code{x}: The x-axis label for each element
#'     \item \code{xmin}, \code{xmax}: The minimum and maximum x-axis values
#'       for the bounding box of each element
#'     \item \code{xmin.t},\code{xmax.t}: The minimum and maximum x-axis values
#'       for the true value and estimate tick marks
#'   }
#'   Can optionally have a column called \code{comp}, which is the
#'   compositional correlation for each element. If \code{corr} is \code{TRUE},
#'   must also have the columns:
#'   \itemize{
#'     \item \code{var.type}: "variance" or "covariance" depending on whether
#'       the element is on the diagonal or not
#'   }
#' @param text.size The size of the plot text to use
#' @param axis.ticks.x A boolean indicating whether to put the coordinate
#'   values on the x-axis.
#' @param axis.text.size.x If specified when \code{x.axis.ticks} is \code{TRUE},
#'   then gives the text size of the x-axis tick labels.
#' @inheritParams get_inference_df
#' @inheritParams sample_mu_prior
#' @inheritParams SPIn::SPIn
#' @inheritParams get_credible_intervals
#' @inheritParams cat_v
#' @return Returns a \code{ggplot} plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw


corr_inference_figure <-
function(inference.df,
         conf=0.95, type="marginal.hpd",
         corr=TRUE, title=NULL, axis.ticks.x=TRUE, axis.text.size.x=NULL,
         ylimits=NULL, segment_size=2, verbose=FALSE, num_level=0,
         text.size=20){
  ## Takes a matrix of true values, a median matrix and a list of
  ##  95% Credible Interval matrices (named "lower" and "upper") and
  ##  makes a CI plot

  banocc::cat_v("Begin corr_inference_figure\n", verbose, num_level)

  ci.name <- paste0(100*conf, "% ",
                    ifelse(type=="marginal.centered", "Central", "HPD"),
                    " CI")
  point.cols <- c("blue", "green", "black")
  point.labels <- c("True Value", "Post. Median", ci.name)
  names(point.labels) <- c("True Value", "Estimate", "95% CI")
  names(point.cols) <- c("True Value", "Estimate", "95% CI")

  if("comp" %in% names(inference.df)){
      point.cols <- c(point.cols, "Observed" = "red")
      point.labels <- c(point.labels, "Observed" = "Observed")
  }
  fill.cols  <- c("variance" = "grey", "covariance"="white")

  banocc::cat_v("Making initial plot...", verbose, num_level + 1)
  gg <- ggplot2::ggplot(ggplot2::aes(y=lower, x=x, xmin=xmin, xmax=xmax, 
                   ymin=lower - 0.1 * abs(lower), 
                   ymax=upper + 0.1 * abs(upper)), 
               data=inference.df) + ggplot2::xlab("(i, j)")
  banocc::cat_v("Done.\n", verbose)

  banocc::cat_v("Adding labels...", verbose, num_level + 1)
  if(!corr){
    gg <- gg + ggplot2::geom_rect(ggplot2::aes(fill=var.type), data=inference.df) + 
      ggplot2::scale_fill_manual("", values=fill.cols) +
      ggplot2::ylab(expression(sigma["(i, j)"]))
  } else {
    gg <- gg + ggplot2::ylab(expression(rho["X,(i, j)"]))
  }
  banocc::cat_v("Done.\n", verbose)

  banocc::cat_v("Adding CI segments and estimates...", verbose, num_level + 1)
  gg <- gg + ggplot2::geom_segment(ggplot2::aes(xend=x, yend=upper, colour="95% CI"),
                                   cex=segment_size) + 
    ggplot2::geom_segment(ggplot2::aes(x=xmin.t, xend=xmax.t, y=true, yend=true, 
                     colour="True Value"), cex=segment_size) + 
    ggplot2::geom_segment(ggplot2::aes(y=median, yend=median, xend=xmax.t, x=xmin.t, 
                     colour="Estimate"), cex=segment_size) +
    ggplot2::scale_colour_manual("", values=point.cols, labels=point.labels) +
    ggplot2::theme_bw()
  if (axis.ticks.x){
      if (is.null(axis.text.size.x)){
          gg <- gg +
              ggplot2::theme(text=ggplot2::element_text(size=text.size),
                             strip.background=ggplot2::element_rect(fill="white",
                                 colour="black"))          
      } else {
          gg <- gg +
              ggplot2::theme(text=ggplot2::element_text(size=text.size),
                             strip.background=ggplot2::element_rect(fill="white",
                                 colour="black"),
                             axis.text.x=ggplot2::element_text(size=axis.text.size.x))          
      }
  } else {
      gg <- gg +
          ggplot2::theme(text=ggplot2::element_text(size=text.size),
                         strip.background=ggplot2::element_rect(fill="white",
                             colour="black"),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank())
  }
  if ("comp" %in% names(inference.df)){
      gg <- gg + ggplot2::geom_segment(ggplot2::aes(y=comp, yend=comp, xend=xmax.t,
                                           x=xmin.t, colour="Observed"),
                                       cex=2)
  }
  banocc::cat_v("Done.\n", verbose)

  if (!is.null(title)){
    banocc::cat_v("Adding title...", verbose, num_level + 1)
    gg <- gg + ggplot2::ggtitle(title)
    banocc::cat_v("Done.\n", verbose)
  }
  if (!is.null(ylimits)){
      banocc::cat_v("Correcting y limits...", verbose, num_level + 1)
      gg <- gg + ggplot2::ylim(ylimits)
      banocc::cat_v("Done.\n", verbose)
  }

  banocc::cat_v("End corr_inference_figure", verbose, num_level)
  return(gg)
}
