#' function to apply all available tree statistics to a single tree
#' @param phylo phylo object
#' @return list with statistics
#' @export
#' @description this function applies all available tree statistics to a
#' single tree, being:
#' \itemize{
#'   \item{gamma}
#'   \item{sackin}
#'   \item{colless}
#'   \item{beta}
#'   \item{blum}
#'   \item{crown age}
#'   \item{tree height}
#'   \item{pigot's rho}
#'   \item{mean branch length}
#'   \item{number of lineages}
#'   \item{laplacian spectrum}
#'   \item{nLTT with empty tree}
#'   \item{phylogenetic diversity}
#'   \item{avgLadder index}
#'   \item{cherries}
#'   \item{ILnumber}
#'   \item{pitchforks}
#'   \item{stairs}
#' }
calc_all_stats <- function(phylo) {

  stats <- list()

  stats$gamma              <- treestats::gamma_statistic(phylo)
  stats$sackin             <- treestats::sackin(phylo)
  stats$colless            <- treestats::colless(phylo)
  stats$beta               <- treestats::beta_statistic(phylo)
  stats$blum               <- treestats::blum(phylo)
  stats$crown_age          <- treestats::crown_age(phylo)
  stats$tree_height        <- treestats::tree_height(phylo)
  stats$pigot_rho          <- treestats::pigot_rho(phylo)
  stats$mean_branch_length <- treestats::mean_branch_length(phylo)
  stats$number_of_lineages <- treestats::number_of_lineages(phylo)
  stats$nltt_base          <- treestats::nLTT_base(phylo)
  stats$phylogenetic_div   <- treestats::phylogenetic_diversity(phylo)
  stats$avgLadder          <- treestats::avgLadder(phylo)
  stats$cherries           <- treestats::cherries(phylo)
  stats$ILnumber           <- treestats::ILnumber(phylo)
  stats$pitchforks         <- treestats::pitchforks(phylo)
  stats$stairs             <- treestats::stairs(phylo)

  temp_stats <- treestats::calc_lapl_spectrum(phylo)

  stats$laplac_spectrum_a  <- temp_stats$asymmetry
  stats$laplac_spectrum_p  <- temp_stats$peakedness

  return(stats)
}
