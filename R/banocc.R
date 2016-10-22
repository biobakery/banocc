#' banocc: A package for Bayesian ANalysis of Compositional Correlation
#'
#' BAnOCC is a package for inferring correlations between features in
#' compositional data, where each sample sums to one. It provides one
#' object, \code{banocc_model} and one function, \code{run_banocc}
#' 
#' @section banocc objects:
#'
#' banocc_model has the \code{stan} model code to be compiled using
#' \code{rstan::stan}.
#' 
#' @section banocc functions:
#'
#' \code{run_banocc} takes a compiled model, and returns the `stanfit` object
#' resulting from a call to \code{rstan::sampling}
#'
#' \code{get_banocc_output} takes a `stanfit` object or the output of
#'   \code{run_banocc} and returns a list with the posterior median and
#'   credible interval estimates
#'
#' @docType package
#' @name banocc
#'
NULL
