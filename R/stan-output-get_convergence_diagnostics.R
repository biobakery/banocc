#' Get several convergence diagnostics for each parameter/chain from the fit.
#'
#' @param samples A set of samples extracted from a \code{stanfit} object with
#'   \code{permuted=FALSE}
#' @param calculate_raftery Boolean: should the raftery diagnostic be calculated?
#' @param stan.output An list that contains the element "Fit", which is a
#'   stanfit object
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#' @return A named list where each element is a matrix of parameters (rows) by
#'   convergence diagnostics (columns). The name of the element represents
#'   the type of convergence diagnostic
#'
#' @importFrom coda geweke.diag
#' @importFrom coda raftery.diag
#' @importFrom coda heidel.diag
#' @importFrom rstan summary

get_convergence_diagnostics <- function(samples, stan.output,
                                        calc_raftery=TRUE,
                                        calc_heidelberg=TRUE,
                                        geweke_fracs=c(0.1, 0.5),
                                        verbose=FALSE, num_level=0){
    banocc::cat_v("Begin get_convergence_diagnostics\n",
                 verbose, num_level=num_level)

    num_chains <- dim(samples)[2]
    # Order rownames so they match for all statistics
    sample_names_order <- dimnames(samples)[[3]][order(dimnames(samples)[[3]])]
    samples <- samples[, , sample_names_order]

    banocc::cat_v("Getting non-zero variance parameters...",
                  verbose, num_level=num_level + 1)
    non_zero_var <- apply(samples, c(2, 3), var) > 1e-16
    banocc::cat_v("Done.\n", verbose)
    

    banocc::cat_v("Begin calculating geweke statistic...",
                 verbose, num_level=num_level + 1)
    geweke.stat <- array(NA, dim=c(dim(samples)[3], num_chains))
    dimnames(geweke.stat) <- list(dimnames(samples)[[3]],
                                  paste0("Geweke_z,chain",
                                         seq_len(num_chains)))
    for (i in seq_len(num_chains)){
        banocc::cat_v(paste0("chain ", i, "..."), verbose)
        non_zero_var_i <- non_zero_var[i, ]
        if (sum(non_zero_var_i) > 0){
            ## for (k in which(non_zero_var_i)){
            ##     print(k)
            ##     print(samples[, i, k])
            ##     print(all(samples[, i, k] > 1e8))
            ##     print(as.logical(non_zero_var_first[i, k] + non_zero_var_second[i, k]))
            ##     print(dimnames(samples)[[3]][k])
            ##     coda::geweke.diag(samples[, i, k],
            ##                       frac1=geweke_fracs[1],
            ##                       frac2=geweke_fracs[2])
            ## }
            geweke.stat[non_zero_var_i, i] <- tryCatch(
                coda::geweke.diag(samples[, i, non_zero_var_i],
                                  frac1=geweke_fracs[1],
                                  frac2=geweke_fracs[2])$z,
                error=function(cond){
                    NA
                }
                )
                
        }
    }
    banocc::cat_v("Done.\n", verbose)


    if (calc_raftery){
        banocc::cat_v("Begin calculating raftery statistic...",
                      verbose, num_level=num_level + 1)
        raftery.stat <- array(NA, dim=c(dim(samples)[3], num_chains))
        dimnames(raftery.stat) <- list(dimnames(samples)[[3]],
                                       paste0("Raftery_I,chain",
                                              seq_len(num_chains)))
        for (i in seq_len(num_chains)){
            banocc::cat_v(paste0("chain ", i, "..."), verbose)
            non_zero_var_i <- non_zero_var[i, ]
            if (sum(non_zero_var_i) > 0){
                raftery.stat[non_zero_var_i, i] <-
                    coda::raftery.diag(samples[, i, non_zero_var_i])$resmatrix[, 4]
            }
        }
        banocc::cat_v("Done.\n", verbose)
    } else {
        raftery.stat <- NULL
    }


    if (calc_heidelberg){
        banocc::cat_v("Begin calculating Heidelberg and Welch statistic...",
                      verbose, num_level=num_level + 1)
        heidelberg.stat <- array(NA, dim=c(dim(samples)[3], 6 * num_chains))
        dimnames(heidelberg.stat)[[1]] <- dimnames(samples)[[3]]
        for (i in seq_len(num_chains)){
            banocc::cat_v(paste0("chain ", i, "..."), verbose)
            non_zero_var_i <- non_zero_var[i, ]
            if (sum(non_zero_var_i)){
                h <- coda::heidel.diag(samples[, i, non_zero_var_i])
                heidelberg.stat[non_zero_var_i, (6 * (i-1) + 1:6)] <- h
                h.names <- colnames(h)
            }
        }
        dimnames(heidelberg.stat)[[2]] <- paste0("),chain", seq_len(num_chains)) %>%
            sapply(., function(s) paste0("Heidelberg_(", h.names, s)) %>%
                as.vector
        banocc::cat_v("Done.\n", verbose)
    } else {
        heidelberg.stat <- NULL
    }

    
    banocc::cat_v("Begin obtaining Rhat statistic...",
                 verbose, num_level=num_level + 1)
    rhat.stat <- array(NA, dim=c(dim(samples)[3], 1))
    summary.fit <- rstan::summary(stan.output$Fit)$summary
    idx_order <- order(rownames(summary.fit))
    rhat.stat[, 1] <- summary.fit[idx_order, "Rhat"]
    dimnames(rhat.stat) <- list(rownames(summary.fit)[idx_order],
                                "Rhat")
    # Order rownames so they match with other stats
    if (any(dimnames(rhat.stat)[[1]] != dimnames(samples)[[3]])){
        stop("Parameter names do not match for Rhat and other statistics")
    }
    banocc::cat_v("Done.\n", verbose)


    banocc::cat_v("End get_convergence_diagnostics.\n",
                 verbose, num_level=num_level)

    return(list(geweke = geweke.stat,
                raftery = raftery.stat,
                heidelberg = heidelberg.stat,
                rhat = rhat.stat))
}

#' Get the convergence diagnostics as a table
#' 
#' @param convergence_diagnostics A list of convergence diagnostics. Each element
#'   should be a matrix of parameters (rows) by diagnostic statistic (columns)
#' @param output_file The name of the file to which to write all diagnostics
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
get_convergence_diagnostics_table <- function(convergence_diagnostics,
                                              params_as_rows=TRUE,
                                              verbose=FALSE, num_level=0){
    banocc::cat_v("Begin get_convergence_diagnostics_table...", verbose,
                  num_level=num_level)
    convergence_matrix <- Reduce(cbind, convergence_diagnostics)

    if (!params_as_rows){
        convergence_matrix <- t(convergence_matrix)
        convergence_matrix_to_file <- as.data.frame(convergence_matrix)
    }
    convergence_table <- as.data.frame(convergence_matrix)
    convergence_table <- cbind(rownames(convergence_matrix),
                               convergence_table)
    if (params_as_rows){
        colnames(convergence_table)[1] <- "Param"
    } else {
        colnames(convergence_table)[1] <- "Statistic"
    }
    rownames(convergence_table) <- NULL
    
    banocc::cat_v("Done.\n", verbose)
    return(convergence_table)
}

#' Write the convergence diagnostics to a file
#'
#' @param params_as_rows Boolean; if TRUE, the rows of the output file will
#'   be the parameters and the columns the diagnostic statistics, while if FALSE
#'   the transpose will be output.
#' @inheritParams get_convergence_diagnostics_table
#' @inheritParams sample_mu_prior
#' @inheritParams cat_v
#'
#' 
output_convergence_diagnostics <- function(convergence_diagnostics, output_file,
                                           params_as_rows = TRUE,
                                           verbose=FALSE, num_level=0){
    banocc::cat_v("Begin output_convergence_diagnostics\n", verbose,
                  num_level=num_level)
    convergence_table <- banocc::get_convergence_diagnostics_table(
        convergence_diagnostics=convergence_diagnostics,
        params_as_rows=params_as_rows, verbose=verbose,
        num_level=num_level+1)
    write.table(convergence_table, file=output_file, quote=FALSE,
                col.names=TRUE, row.names=FALSE, sep="\t")
    banocc::cat_v("End output_convergence_diagnostics\n", verbose,
                  num_level=num_level)
}

#' Read in convergence diagnostics from a file
#'
#' @param input_file The name of the file with the convergence diagnostics.
#'   These are assumed to have parameter names as the first column and
#'   the convergence diagnostics on the remaining columns
#' @inheritParams readr::read_delim
#' @importFrom magrittr %>%
read_in_convergence_diagnostics <- function(input_file, delim="\t",
                                            verbose=FALSE,
                                            num_level=0){
    banocc::cat_v("Begin read_in_convergence_diagnostics\n", verbose=verbose,
                  num_level=num_level)
    diagnostic_stats <- readr::read_delim(input_file, delim=delim)
    if (colnames(diagnostic_stats) %>% tail(1) %>% is.na){
        banocc::cat_v("Adding Rhat colname...", verbose=verbose,
                      num_level=num_level + 1)
        colnames(diagnostic_stats)[ncol(diagnostic_stats)] <- "Rhat,ChainAll"
        banocc::cat_v("Done.\n", verbose)
    }
    num_chains <- diagnostic_stats %>% colnames %>%
        stringr::str_extract("chain\\d") %>% na.omit %>%
        stringr::str_split("chain") %>% unlist %>%
        stringr::str_extract("\\d") %>% na.omit %>%
        as.numeric %>% max

    banocc::cat_v("Are Heidelberg colnames in correct order?...",
                  verbose=verbose, num_level=num_level + 1)
    if (num_chains > 1){
        stest_cols <- diagnostic_stats %>% colnames %>%
            stringr::str_detect("(stest)")
        chain1_cols <- diagnostic_stats %>% colnames %>%
            stringr::str_detect("chain1")
        chain2_cols <- diagnostic_stats %>% colnames %>%
            stringr::str_detect("chain2")
        heidelberg_ordering <- (stest_cols * (chain1_cols + chain2_cols)) %>%
            rle
        heidelberg_incorrect_order <- (heidelberg_ordering$lengths > 1 &
                                       heidelberg_ordering$values == 1) %>% any
    } else {
        heidelberg_incorrect_order <- FALSE
    }
    banocc::cat_v(ifelse(heidelberg_incorrect_order, "No...", "Yes.\n"),
                  verbose=verbose)

    if (heidelberg_incorrect_order){
        banocc::cat_v("Correcting...", verbose=verbose)
        heidelberg_idx <- grep("eidelberg", colnames(diagnostic_stats))
        heidelberg_names <- colnames(diagnostic_stats)[heidelberg_idx]
        h_names <- stringr::str_extract(heidelberg_names, "\\(.*\\)") %>%
            stringr::str_sub(2L,-2L) %>%
            unique
        heidelberg_names <- paste0("),chain", seq_len(num_chains)) %>%
            sapply(., function(s) paste0("Heidelberg_(", h_names, s)) %>%
            as.vector
        colnames(diagnostic_stats)[heidelberg_idx] <- heidelberg_names
        banocc::cat_v("Done.\n", verbose=verbose)
    }
    banocc::cat_v("End read_in_convergence_diagnostics\n", verbose=verbose,
                  num_level=num_level)
    return(diagnostic_stats)
}
