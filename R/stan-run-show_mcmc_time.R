#' Display the total time it took to run each chain of a \code{stan} run
#'
#' @param Fit.all The captured output (using \code{\link{mycapture}}) of the
#'   generation of a \code{stanfit} object
#' @return prints the lines of printed output corresponding to the time
#' @importFrom stringr str_c

show_mcmc_time <-
function(Fit.all){
    time.idx <- grep("(Total)", Fit.all$print.output)
    cat(stringr::str_c(Fit.all$print.output[time.idx], collapse="\n"))
}
