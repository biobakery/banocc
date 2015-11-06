#' Warn if any values are too large/small for exponentiation
#'
#' @param sample a list of samples from a particular distribution
#' @param verbose Boolean: print the indices and values of the samples that are
#'   too large or small
#' @examples
#' warn_exponent(1)
#' warn_exponent(710)
#' warn_exponent(list(-800, 0, 800, 12, -1e-16))
#' warn_exponent(list(-800, 0, 800, 12, -1e-16), verbose=TRUE)

warn_exponent <-
    function(sample, verbose=FALSE){
        banocc::warn_large_sample(sample, verbose)
        banocc::warn_small_sample(sample, verbose)
    }

warn_large_sample <-
function(sample, verbose=FALSE){
    warning.text <- paste0(" samples have > 709 or < -745 ",
                           "exponentiating them might give Inf or zero.")
    large.sample <- unlist(lapply(sample,
                                  function(s) any((s > 709) + (s < -745))))
    if (any(large.sample)){
        warning(paste0(sum(large.sample), warning.text))
        if (verbose){
            print(which(large.sample))
            print(sample[large.sample])
        }
    }
}

warn_small_sample <-
function(sample, verbose=FALSE){
    warning.text <- paste0(" samples have values between -5e-8 and 5e-7; ",
                           "exponentiating them might give 1.00.")
    small.sample <- unlist(lapply(sample,
                                  function(s) any((s < 5e-7) * (s >= 0) +
                                                  (s > -5e-8) * (s <= 0))))
    if(any(small.sample)){
        warning(paste0(sum(small.sample), warning.text))
        if (verbose){
            print(which(small.sample))
            print(sample[small.sample])
        }
    }
}
