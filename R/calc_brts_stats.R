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
  stats$pigot_rho              <- treestats::pigot_rho(phylo)
  stats$nltt_base              <- treestats::nLTT_base(phylo)

  stats$treeness               <- treestats::treeness(phylo)

  stats$mean_branch_length     <- treestats::mean_branch_length(phylo)
  stats$var_branch_length      <- treestats::var_branch_length(phylo)

  stats$mean_branch_length_int <- treestats::mean_branch_length_int(phylo)
  stats$mean_branch_length_ext <- treestats::mean_branch_length_ext(phylo)
  stats$var_branch_length_int  <- treestats::var_branch_length_int(phylo)
  stats$var_branch_length_ext  <- treestats::var_branch_length_ext(phylo)
  stats <- unlist(stats)
  stats <- stats[order(names(stats))]

  return(stats)
}
