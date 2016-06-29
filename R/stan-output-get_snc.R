# Calculates the scaled neighborhood criterion for the output
#
# @inheritParams get_posterior_quantities
# @inheritParams cat_v
# @return Returns a list of the SNC values for each parameter in
#   \code{parameter.names}.
#

get_snc <- function(posterior_samples, parameter.names=c("m", "S"),
                    verbose=FALSE, num_level=0){
    cat_v("Begin get_snc\n", verbose, num_level=num_level)
    names(parameter.names) <- parameter.names
    snc_values <- lapply(parameter.names, function(name){
        is.mat <- length(dim(posterior_samples[[name]])) == 3
        is.vec <- length(dim(posterior_samples[[name]])) == 2
        if (is.mat){
            apply(posterior_samples[[name]], c(2, 3), calc_snc)
        } else if (is.vec){
            apply(posterior_samples[[name]], 2, calc_snc)
        } else {
            calc_snc(as.vector(posterior_samples[[name]]))
        }
    })
    cat_v("End get_snc\n", verbose, num_level=num_level)
    return(snc_values)
}

# Calculate the scaled neighborhood criterion for a vector of
#   posterior samples
calc_snc <- function(sample_vec){
    if (!is.vector(sample_vec) || !(mode(sample_vec)=="numeric")){
        stop("'sample_vec' must be a numeric vector")
    }
    if (sd(sample_vec) == 0){
        return(NA)
    } else {
        return(sum(abs(sample_vec) <= sd(sample_vec)) / length(sample_vec))
    }
}
