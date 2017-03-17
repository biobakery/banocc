# Cat output only if verbose
#
# @param s_string The string to output if \code{verbose} is \code{TRUE}
# @param verbose Print informative statements as the function executes?xo
# @param num_level The number of the level (will determine the number of spaces to add)
#

cat_v <- function(s_string, verbose, num_level=0){
   if (num_level > 0) s_string <- paste0(c(rep("  ", num_level), s_string),
                                         collapse="")
   if(verbose) message(s_string)
}

# Check if a square matrix is symmetric
#
# @param mat The square matrix to be checked
#
is.symmetric <- function(mat){
    if (nrow(mat) != ncol(mat)){
        stop("'mat' must be square.")
    }

    return(all(t(mat) == mat))
}


#  Check if all the elements of a vector are the same within a certain
#   tolerance
# @param vec A numeric vector
# @param tol The tolerance
# 

vec.equal <- function(vec, tol=1e-6){
    vec.mat <- matrix(rep(vec, length(vec)), ncol=length(vec))
    diff.mat <- vec.mat - t(vec.mat)
    return(all(abs(diff.mat) < tol))
}

# Get elements of a vector corresponding to a particular pattern
#
# @param print.output A vector of strings
# @inheritParams base::grep
# @return A vector of strings, each one of which contains \code{pattern}

grep_from_output <-
function(print.output, pattern){
    return(print.output[grep(pattern, print.output)])
}
