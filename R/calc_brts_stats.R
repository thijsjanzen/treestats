#' function to apply all tree statistics related to branching times
#' to a single tree.
#' @param phylo phylo object
#' @return list with statistics
#' @export
#' @description this function applies all tree statistics based on
#' branching times to a single tree, being:
#' \itemize{
#'   \item{gamma}
#'   \item{pigot's rho}
#'   \item{mean branch length}
#'   \item{nLTT with empty tree}
#' }
#'
calc_brts_stats <- function(phylo, normalize = FALSE) {

  stats <- list()

  stats$gamma              <- treestats::gamma_statistic(phylo)
  stats$pigot_rho          <- treestats::pigot_rho(phylo)
  stats$mean_branch_length <- treestats::mean_branch_length(phylo)
  stats$nltt_base          <- treestats::nLTT_base(phylo)

  return(stats)
}
