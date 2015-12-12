#' Tidy the convergecne diagnostics table
#'
#' @inheritParams manipulate_convergence_diagnostics
#' @param important_stats A vector of the names of statistics that will be
#'   plotted; used to filter stat names
#' @param important_params A vector of the names of the key parameters to
#'   look at; used to make an indicator of important variables.
#' @inheritParams cat_v
#' @inheritParams sample_mu_prior
#' 
plot_convergence_diagnostics <- function(convergence_diagnostics_table,
                                         filter=TRUE,
                                         important_stats=c("Rhat", "Geweke_z"),
                                         important_params=c("L", "sigma","mu",
                                             "ln_Rho", "lp__"),
                                         verbose=FALSE, num_level=0){
    banocc::cat_v("Begin plot_convergence_diagnostics.\n", verbose=verbose,
                  num_level=num_level)
    convergence_diagnostics_df <- banocc::get_convergence_diagnostics_df(
        convergence_diagnostics_table = convergence_diagnostics_table,
        filter=filter,
        verbose=verbose, num_level=num_level + 1
        )
    less_important_params <- setdiff(unique(convergence_diagnostics_df$param_name),
                                    important_params)
    convergence_diagnostics_df <- convergence_diagnostics_df %>%
        dplyr::mutate(important_param=ifelse(param_name %in% important_params,
                          "Important", "Less Important"),
                      param_name_factor = factor(param_name,
                          levels=c(important_params,
                              less_important_params))) %>%
        dplyr::filter(stat_name %in% important_stats)

    Rhat_lines <- data.frame(stat_name_factor=factor("Rhat"),
                             yintercept=c(1.1, 1.2),
                             xintercept=c(1.1, 1.2),
                             severity=c("High", "Extreme"))
    probs <- c(0.025, 0.005)
    Geweke_lines <- data.frame(stat_name_factor=factor("Geweke_z"),
                               yintercept=as.vector(sapply(c(-1, 1), "*",
                                   qnorm(probs))),
                               xintercept=as.vector(sapply(c(-1, 1), "*",
                                   qnorm(probs))),
                               severity=rep(c("High", "Extreme"), 2))

    diagnostic_plot <- ggplot2::ggplot(data=convergence_diagnostics_df,
                                       ggplot2::aes(x=param_name_factor,
                                                    y=stat_value)) +
        ggplot2::facet_grid(stat_name_factor ~ important_param, scales="free") + 
        ## ggplot2::facet_wrap(~stat_name_factor, scales="free_y") +
        ## ggplot2::geom_point(ggplot2::aes(colour=chain_factor)) +
        ggplot2::geom_boxplot(ggplot2::aes(colour=chain_factor)) +
        ggplot2::geom_hline(data=rbind(Rhat_lines, Geweke_lines),
                            show_guide=TRUE,
                            ggplot2::aes(yintercept=yintercept,
                                         linetype=severity))

    banocc::cat_v("End plot_convergence_diagnostics.\n", verbose=verbose,
                  num_level=num_level)
    return(diagnostic_plot)
}


#' Manipulate the convergence diagnostics to get useful tables
#'
#' @param convergence_diagnostics_table A table of convergence diagnostics;
#'   columns are the diagnostics, rows are the parameters
#' @param convergence_diagnostics_table_tidy A tidy version of the
#'   convergence diagnostics table as output by
#'   \code{tidy_convergence_diagnostics}
#' @param filter Boolean: perform the filtering?
#' @inheritParams cat_v
#' @inheritParams sample_mu_prior
#' @name manipulate_convergence_diagnostics

#' @rdname manipulate_convergence_diagnostics
tidy_convergence_diagnostics <- function(convergence_diagnostics_table,
                                         verbose=FALSE, num_level=0){

    banocc::cat_v("Begin tidy_convergence_diagnostics...",
                  verbose=verbose, num_level=num_level + 1)
    convergence_diagnostics_table_tidy <- convergence_diagnostics_table %>%
        dplyr::mutate(param_name=stringr::str_extract(Param, "[a-zA-Z_]+"),
                      param_idx=stringr::str_extract(Param, "[\\d,]+")) %>%
        tidyr::gather(stat_ID, stat_value,
                      -c(Param, param_name, param_idx)) %>%
        dplyr::mutate(stat_name=stringr::str_extract(stat_ID,
                          "[a-zA-Z_\\(\\)]+"),
                      chain=stringr::str_extract(stat_ID, "chain[\\dAll]+"),
                      param_idx1=stringr::str_split_fixed(param_idx, ",",
                          n=2)[, 1],
                      param_idx2=stringr::str_split_fixed(param_idx, ',',
                          n=2)[, 2])
    banocc::cat_v("Done.\n", verbose=verbose)
    return(convergence_diagnostics_table_tidy)
}

#' @rdname manipulate_convergence_diagnostics
filter_convergence_diagnostics <- function(convergence_diagnostics_table_tidy,
                                           verbose=FALSE, num_level=0){
    banocc::cat_v("Begin filter_convergence_diagnostics...",
                  verbose=verbose, num_level=num_level+1)
    convergence_diagnostics_table_filter <-
        convergence_diagnostics_table_tidy %>%
        dplyr::filter(param_idx1>=param_idx2) %>%
        dplyr::filter(param_idx1 != param_idx2 | !is.na(stat_value))
    banocc::cat_v("Done.\n", verbose=verbose)
    return(convergence_diagnostics_table_filter)
}

#' @rdname manipulate_convergence_diagnostics
get_convergence_diagnostics_df <- function(convergence_diagnostics_table,
                                           filter=TRUE,
                                    verbose=FALSE, num_level=0){
    banocc::cat_v("Begin get_convergence_diagnostics_df.\n", verbose,
                  num_level=num_level)

    convergence_diagnostics_table_tidy <-
        banocc::tidy_convergence_diagnostics(
            convergence_diagnostics_table=convergence_diagnostics_table,
            verbose=verbose, num_level=num_level + 1)
    if (filter){
        convergence_diagnostics_table_filter <-
            banocc::filter_convergence_diagnostics(
                convergence_diagnostics_table_tidy=convergence_diagnostics_table_tidy,
                verbose=verbose, num_level=num_level + 1
                )
    } else {
        convergence_diagnostics_table_filter <-
            convergence_diagnostics_table_tidy
    }
    convergence_diagnostics_df <- convergence_diagnostics_table_filter %>%
            dplyr::mutate(stat_name_factor = factor(stat_name,
                          levels=c("Rhat", "Geweke_z",
                              "Heidelberg_(stest)", "Heidelberg_(start)",
                              "Heidelberg_(pvalue)",
                              "Heidelberg_(htest)", "Heidelberg_(mean)",
                              "Heidelberg_(halfwidth)")),
                       chain_factor = factor(chain,
                           levels=sort(unique(chain))))


    banocc::cat_v("End get_convergence_diagnostics_df.\n", verbose=verbose,
                  num_level=num_level)

    return(convergence_diagnostics_df)
}
