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
  stats$avgLadder          <- treestats::avgLadder(phylo) #nolint
  stats$cherries           <- treestats::cherries(phylo)
  stats$ILnumber           <- treestats::ILnumber(phylo) #nolint
  stats$pitchforks         <- treestats::pitchforks(phylo)
  stats$stairs             <- treestats::stairs(phylo)

  temp_stats <- treestats::calc_lapl_spectrum(phylo)

  stats$laplac_spectrum_a  <- temp_stats$asymmetry
  stats$laplac_spectrum_p  <- temp_stats$peakedness
  stats$laplac_sepctrum_e  <- log(temp_stats$principal_eigenvalue)

  stats$B1           <- treestats::b1(phylo)
  stats$B2           <- treestats::b2(phylo)
  stats$aPP          <- treestats::area_per_pair(phylo)
  stats$aLD          <- treestats::average_leaf_depth(phylo)
  stats$Ibased       <- treestats::mean_i(phylo)
  stats$ewColless    <- treestats::ew_colless(phylo)
  stats$maxDelW      <- treestats::max_del_width(phylo)
  stats$maxDepth     <- treestats::max_depth(phylo)
  stats$maxWidth     <- treestats::max_width(phylo)
  stats$rogers       <- treestats::rogers(phylo)
  stats$stairs2      <- treestats::stairs2(phylo)
  stats$totCoph      <- treestats::tot_coph(phylo)
  stats$varLeafDepth <- treestats::var_leaf_depth(phylo)
  stats$symNodes     <- treestats::sym_nodes(phylo)

  stats$mpd          <- treestats::mean_pair_dist(phylo)
  stats$psv          <- treestats::psv(phylo)
  stats$vpd          <- treestats::var_pair_dist(phylo)
  stats$mntd         <- treestats::mntd(phylo)
  stats$J            <- treestats::entropy_j(phylo)

  return(stats)
}
