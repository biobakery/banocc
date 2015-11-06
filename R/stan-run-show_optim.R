#' Display relevant parts of the printed outputs from running
#'   \code{stan::optimizing}
#'
#' @inheritParams show_mcmc_time
#' @param exit Boolean: print the exit status of the run?
#' @param startEnd Boolean: print the starting and ending portions of the run?
#' @name show_optim
#'
#' @importFrom stringr str_c

#' @rdname show_optim
show_optim_end <-
function(Fit.all){
    header.idx <- grep("Iter", Fit.all$print.output)
    end.idx    <- grep("terminated", Fit.all$print.output) - 1
    tables <- paste0(Fit.all$print.output[header.idx], "\n",
                     Fit.all$print.output[end.idx])
    cat(stringr::str_c(tables, collapse="\n"))
}

#' @rdname show_optim
show_optim_exit <-
function(Fit.all){
    exit.idx <- grep("terminated", Fit.all$print.output)
    cat(stringr::str_c(Fit.all$print.output[exit.idx + c(0, 1)], collapse="\n"))
}

#' @rdname show_optim
show_optim_start <-
function(Fit.all){
    start.idx <- grep("^initial", Fit.all$print.output)
    cat(stringr::str_c(Fit.all$print.output[start.idx], collapse="\n"))
}

#' @rdname show_optim
show_selected_optim_output <-
function(Fit.all, exit=TRUE, startEnd=TRUE){
    if(startEnd){
        banocc::show_optim_start(Fit.all)
        cat("\n")
        banocc::show_optim_end(Fit.all)
        cat("\n")
    }
    if(exit){
        banocc::show_optim_exit(Fit.all)
    }
}
