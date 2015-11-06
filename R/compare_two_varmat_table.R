#' Make a latex table comparing one symmetric matrix with another
#'
#' @param Mat.lower One symmetric matrix
#' @param Mat.upper Second symmetric matrix
#' @return Returns a string that, when used with \code{cat}, forms a LaTeX table
#'   comparing \code{Mat.lower} (in the lower triangle) with
#'   \code{Mat.upper} (in the upper triangle).
#' @importFrom stringr str_c

compare_two_varmat_table <-
function(Mat.lower,Mat.upper){
  p.true <- ncol(Mat.lower)
  comparison_portion <- stringr::str_c(
    c(
      stringr::str_c(c(paste0("\\multirow{", p.true + 1, "}{*}{$\\bhSigma$}"),
              "", signif(Mat.lower[1, ], 3)), collapse=" & "),
        sapply( seq_len(p.true - 1) + 1, function(j){
            stringr::str_c(c("", round(Mat.upper[j - 1, seq_len(j - 1)], 3),
                    "",
                    signif(Mat.lower[j, -seq_len(j - 1)], 3)), collapse=" & ")
        } ),
      stringr::str_c(c("", round(Mat.upper[p.true, ], 2), ""),
                     collapse=" & ")
    ),
    collapse=paste0("\\\\\n\\cline{2-", p.true + 2, "}\n")
  )
  tail_portion <- "\\\\\n\\hline\n\\end{tabular}"
  head_portion <- paste0("\\begin{tabular}{|",
                         stringr::str_c(rep("c", p.true + 2),
                                        collapse="|"), "|}\n",
                         "\\hline\n",
                         "& \\multicolumn{", p.true + 1,
                         "}{c|}{$\\bSigma$}\\\\\n",
                         "\\hline\n")
  return(stringr::str_c(head_portion, comparison_portion, tail_portion))
}
