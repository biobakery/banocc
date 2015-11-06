#' Generate a melted dataset by combining two correlation matrices
#'
#' @param Rho.lower The lower correlation matrix
#' @param Rho.upper The upper correlation matrix; must match \code{Rho.lower} in dimension
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @return A data frame with the following columns: \code{Var1} and \code{Var2}
#'   indicating the column/row indices of the correlation matrix;
#'   \code{cor.value} indicating the value of the correlation at those indices;
#'   \code{triangle} with values \code{"lower"} for the lower triangle elements
#'   and \code{"upper"} for the upper triangle elements.
#' @examples
#' Rho.lower <- matrix(c(1, 0, 0, 0.5, 0,
#'                       0, 1, 0, 0, 0,
#'                       0, 0, 1, 0, 0,
#'                       0.5, 0, 0, 1, -0.5,
#'                       0, 0, 0, -0.5, 1), ncol=5)
#' X <- rmvnorm(5, rep(0, 5), Rho.lower)
#' Rho.upper <- cov2cor(X %*% t(X))
#' melt_two_corr_mat(Rho.lower, Rho.upper)

melt_two_corr_mat <- 
  function(Rho.lower, Rho.upper, verbose=FALSE, num_level=0){
    if (any(dim(Rho.lower)!=dim(Rho.upper))){
      stop("Dimensions of Rho.lower and Rho.upper must match.")
    }
    banocc::cat_v("Begin melt_two_corr_mat\n", verbose, num_level)
    banocc::cat_v("Melting both datasets...", verbose, num_level + 1)
    p <- nrow(Rho.lower)
    
    rownames(Rho.lower) <- paste0("f", seq(1, p))
    
    Rho.lower.melt  <- banocc::get_melt_dataset(Rho.lower, cor_mat=TRUE)
    Rho.upper.melt  <- banocc::get_melt_dataset(Rho.upper, cor_mat=TRUE)
    
    idx.mat <- sapply(rownames(Rho.lower), 
                      function(a) paste(a, rownames(Rho.lower)))
    
    Rho.melt <- banocc::mymerge(Rho.lower.melt, Rho.upper.melt)
    names(Rho.melt) <- c("Var1", "Var2", "lower", "upper")
    Rho.melt$y <- factor(Rho.melt$Var2,
                         levels=rev(levels(Rho.melt$Var2)), ordered=TRUE)
    Rho.melt$triangle <- ""
    banocc::cat_v("Done.\n", verbose)
    
    banocc::cat_v("Combining correlations...", verbose, num_level + 1)
    Rho.melt$cor.value <- NA
    Rho.melt.idx       <- paste(Rho.melt$Var1, Rho.melt$Var2)
    Rho.melt.l.idx     <- which(Rho.melt.idx %in% idx.mat[lower.tri(idx.mat)])
    Rho.melt.u.idx     <- which(Rho.melt.idx %in% idx.mat[upper.tri(idx.mat)])
    
    Rho.melt$cor.value[Rho.melt.l.idx] <- Rho.melt$lower[Rho.melt.l.idx]
    Rho.melt$cor.value[Rho.melt.u.idx] <- Rho.melt$upper[Rho.melt.u.idx]
    
    Rho.melt$triangle[Rho.melt.l.idx] <- "lower"
    Rho.melt$triangle[Rho.melt.u.idx] <- "upper"

    Rho.melt$lower <- NULL
    Rho.melt$upper <- NULL
    
    Rho.melt <- na.omit(Rho.melt)
    banocc::cat_v("Done.\n", verbose)
    banocc::cat_v("End melt_two_corr_mat\n", verbose, num_level)

    return(Rho.melt)
  }


#' Generate a melted correlation dataset from a matrix using either a dataset
#'   or a correlation matrix
#'
#' @param mat Either a correlation matrix, or a dataset from which the
#'   correlation will be calculated. If a dataset, samples must be rows.
#' @param cor_mat Boolean: is \code{mat} a correlation matrix?
#' @return A melted data frame (as from \code{\link{reshape2::melt}})
#'
#' @importFrom reshape2 melt

get_melt_dataset <- function(mat, cor_mat=FALSE){
  colnames(mat) <- paste0("f", seq(1, ncol(mat)))
  if (cor_mat){
    rownames(mat) <- colnames(mat)
    mat.melt <- reshape2::melt(mat)
  } else {
    mat.melt <- reshape2::melt(cor(mat))
  }
  return(mat.melt)
}

#' A function to merge two melted datasets
#'
#' @param x,y melted datasets, each with columns \code{Var1} and \code{Var2}
#' @return A combined data frame combined by combinations of \code{Var1} and
#'   \code{Var2}
#' @examples
#' x <- reshape2::melt(diag(5))
#' y <- reshape2::melt(diag(5) + 0.01 - diag(rep(5, 0.01)))
#' mymerge(x, y)
mymerge <-
function(x, y){
  merge(x, y, by=c("Var1", "Var2"))
}
