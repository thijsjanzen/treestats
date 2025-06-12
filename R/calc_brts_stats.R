#' Apply all tree statistics related to branching times to a single tree.
#' @param phylo phylo object
#' @return list with statistics
#' @export
#' @description this function applies all tree statistics based on
#' branching times to a single tree (more or less ignoring topology), being:
#' \itemize{
#'   \item{gamma}
#'   \item{pigot's rho}
#'   \item{mean branch length}
#'   \item{nLTT with empty tree}
#'   \item{treeness}
#'   \item{var branch length}
#'   \item{mean internal branch length}
#'   \item{mean external branch length}
#'   \item{var internal branch length}
#'   \item{var external branch length}
#' }
#'
calc_brts_stats <- function(phylo) {

  stats <- list()

  stats$gamma                  <- try_stat(phylo, treestats::gamma_statistic)
  stats$pigot_rho              <- try_stat(phylo, treestats::pigot_rho)
  stats$nltt_base              <- try_stat(phylo, treestats::nLTT_base)

  stats$treeness               <- try_stat(phylo, treestats::treeness)

  stats$mean_branch_length     <- try_stat(phylo, treestats::mean_branch_length)
  stats$var_branch_length      <- try_stat(phylo, treestats::var_branch_length)

  stats$mean_branch_length_int <- try_stat(phylo,
                                           treestats::mean_branch_length_int)
  stats$mean_branch_length_ext <- try_stat(phylo,
                                           treestats::mean_branch_length_ext)
  stats$var_branch_length_int  <- try_stat(phylo,
                                           treestats::var_branch_length_int)
  stats$var_branch_length_ext  <- try_stat(phylo,
                                           treestats::var_branch_length_ext)
  stats <- unlist(stats)
  stats <- stats[order(names(stats))]

  return(stats)
}
