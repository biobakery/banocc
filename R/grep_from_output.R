#' Get elements of a vector corresponding to a particular pattern
#'
#' @param print.output A vector of strings
#' @inheritParams base::grep
#' @return A vector of strings, each one of which contains \code{pattern}

grep_from_output <-
function(print.output, pattern){
    return(print.output[grep(pattern, print.output)])
}
