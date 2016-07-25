# Take a list of arrays of posterior quantities and make the elements
#   lists instead of arrays
#
# @param array_list The list of arrays of posterior quantities
# @inheritParams get_posterior_quantiles
#

make_list <- function(array_list, parameter.names, posterior_samples,
                      elt_length=2, elt_names=NULL){
    list_list <- lapply(parameter.names, function(name){
        new.vec <- vector("list", length=elt_length)
        names(new.vec) <- elt_names
        return(new.vec)
    })
    for (i in seq_len(elt_length)){
        for (name in parameter.names){
            is.mat <- length(dim(posterior_samples[[name]])) == 3
            is.vec <- length(dim(posterior_samples[[name]])) == 2
            if(is.mat){
                list_list[[name]][[i]] <- array_list[[name]][i, , ]
            } else if (is.vec){
                list_list[[name]][[i]] <- array_list[[name]][i, ]
            } else {
                list_list[[name]][[i]] <- unname(array_list[[name]][i])
            }
        }
    }
    return(list_list)
}
