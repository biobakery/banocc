#' Cat output only if verbose
#'
#' @param verbose Print informative statements as the function executes?xo
#' @param num_level The number of the level (will determine the number of spaces to add)
#'
#' @importFrom stringr str_c
#' @importFrom stringr str_dup

cat_v <- function(s_string, verbose, num_level=0){
   if (num_level > 0) s_string <- stringr::str_c(stringr::str_dup("  ", num_level), s_string)
   if(verbose) cat(s_string)
}

#' A function that captures any printed output from an evaluated R expression
#'   and returns it as an element of a list
#'
#' @param ... The expression to be evaluated
#' @param file An optional string or connection to which to write the output
#' @return A list; \code{output} is the actual output of \code{...}, and
#'   \code{print.output} is the printed output from evaluating \code{...}.

mycapture <-
function(..., file=NULL){
      closeit <- TRUE
        if (is.null(file))
                file <- textConnection("rval", "w", local = TRUE)
        else if (is.character(file))
                file <- file(file, if (append)
                                   "a"
                                              else "w")
        else if (inherits(file, "connection")) {
                if (!isOpen(file))
                          open(file, if (append)
                                       "a"
                                          else "w")
                    else closeit <- FALSE
            }
        else stop("'file' must be NULL, a character string or a connection")
        sink(file)
        output <- list(output=eval(...))
        sink()
        if (closeit){
                close(file)
            }
        output$print.output = rval
        return(output)
  }


#' Check if a square matrix is symmetric
#'
#' @param mat The square matrix to be checked
#'
is.symmetric <- function(mat){
    if (nrow(mat) != ncol(mat)){
        stop("'mat' must be square.")
    }

    return(all(t(mat) == mat))
}


#' Check if all the elements of a vector are the same within a certain
#'   tolerance
#' @param vec A numeric vector
#' @param tol The tolerance
#' 

vec.equal <- function(vec, tol=1e-6){
    vec.mat <- matrix(rep(vec, length(vec)), ncol=length(vec))
    diff.mat <- vec.mat - t(vec.mat)
    return(all(abs(diff.mat) < tol))
}

#' Get elements of a vector corresponding to a particular pattern
#'
#' @param print.output A vector of strings
#' @inheritParams base::grep
#' @return A vector of strings, each one of which contains \code{pattern}

grep_from_output <-
function(print.output, pattern){
    return(print.output[grep(pattern, print.output)])
}
