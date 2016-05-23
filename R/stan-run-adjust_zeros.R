# Replace all zeros in a data frame with a small constant; possibly
#   renormalize by samples
#
# @param data The data frame to be adjusted; assumed to have rows for samples
# @param renorm Whether to renormalize the data after adjusting the zeros
# @param zero_adj The amount to use instead of zero.
# @inheritParams cat_v

adjust_zeros <- function(data, renorm=TRUE, zero_adj=0.0001, verbose=FALSE, num_level=0){
   cat_v("Begin adjust_zeros...", verbose, num_level=num_level)
   data.adj <- as.matrix(data)
   data.adj[which(data.adj==0)] <- zero_adj
   if (renorm) {
       data.adj <- data.adj / rowSums(data.adj)
   }
   if (is.data.frame(data)){
       data.adj <- as.data.frame(data.adj)
   }
   cat_v("Done.\n", verbose)
   return(data.adj)
}
