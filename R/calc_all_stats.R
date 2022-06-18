#' function to apply all available tree statistics to a single tree
#' @param phylo phylo object
#' @return list with statistics
#' @export
#' @description this function applies all tree statistics available in
#' this package to a
#' single tree, being:
#' \itemize{
#'   \item{gamma}
#'   \item{sackin (yule corrected)}
#'   \item{colless (yule corrected)}
#'   \item{Aldous' beta statistic}
#'   \item{Blum}
#'   \item{crown age}
#'   \item{tree height}
#'   \item{pigot's rho}
#'   \item{mean branch length}
#'   \item{number of lineages}
#'   \item{nLTT with empty tree}
#'   \item{phylogenetic diversity}
#'   \item{avgLadder index}
#'   \item{cherries}
#'   \item{ILnumber}
#'   \item{pitchforks}
#'   \item{stairs}
#'   \item{stairs2}
#'   \item{laplacian spectrum}
#'   \item{B1}
#'   \item{B2}
#'   \item{area per pair (aPP) }
#'   \item{average leaf depth (aLD)}
#'   \item{I statistic}
#'   \item{ewColless}
#'   \item{max Delta Width (maxDelW)}
#'   \item{maximum of Depth}
#'   \item{variance of Depth}
#'   \item{maximum Width}
#'   \item{Rogers}
#'   \item{total Cophenetic distance}
#'   \item{symmetry Nodes}
#'   \item{mean of pairwise distance (mpd)}
#'   \item{variance of pairwise distance (vpd)}
#'   \item{Phylogenetic Species Variability (psv)}
#'   \item{mean nearest taxon distance (mntd)}
#'   \item{J statistic of entropy}
#'   \item{rquartet index}
#'   \item{Wiener index}
#'   \item{max betweenness}
#'   \item{max closeness}
#'   \item{diameter, without branch lenghts}
#'   \item{maximum eigen vector value}
#' }
#'
#' For the Laplacian spectrum properties, four properties of the eigenvalue
#' distribution are returned: 1) asymmetry, 2) peakedness, 3) log(principal
#' eigenvalue) and 4) eigengap.
#' Please notice that for some very small or very large trees, some of the
#' statistics can not be calculated. The function will report an NA for this
#' statistic, but will to break, to facilitate batch analysis of large numbers
#' of trees.
calc_all_stats <- function(phylo) {

  stats <- list()

  stats$gamma              <- treestats::gamma_statistic(phylo)
  stats$sackin             <- treestats::sackin(phylo, normalization = "yule")
  stats$colless            <- treestats::colless(phylo, normalization = "yule")
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
  stats$il_number           <- treestats::ILnumber(phylo) #nolint
  stats$pitchforks         <- treestats::pitchforks(phylo)
  stats$stairs             <- treestats::stairs(phylo)

  temp_stats <- tryCatch(expr = {treestats::calc_lapl_spectrum(phylo) },
                         error = function(e) {return(NA) })

  if (length(temp_stats) == 5) {
    stats$laplac_spectrum_a  <- temp_stats$asymmetry
    stats$laplac_spectrum_p  <- temp_stats$peakedness
    stats$laplac_spectrum_e  <- log(temp_stats$principal_eigenvalue)
    stats$laplac_spectrum_g  <- temp_stats$eigengap[[1]]
  } else {
    stats$laplac_spectrum_a  <- NA
    stats$laplac_spectrum_p  <- NA
    stats$laplac_spectrum_e  <- NA
    stats$laplac_spectrum_g  <- NA
  }


  stats$b1           <- treestats::b1(phylo)
  stats$b2           <- treestats::b2(phylo)
  stats$area_per_pair          <- treestats::area_per_pair(phylo)
  stats$average_leaf_depth          <- treestats::average_leaf_depth(phylo)
  stats$i_stat       <- treestats::mean_i(phylo)
  stats$ew_colless    <- treestats::ew_colless(phylo)
  stats$max_del_width      <- treestats::max_del_width(phylo)
  stats$max_depth     <- treestats::max_depth(phylo)
  stats$max_width     <- treestats::max_width(phylo)
  stats$rogers       <- treestats::rogers(phylo)
  stats$stairs2      <- treestats::stairs2(phylo)
  stats$tot_coph      <- tryCatch(expr = {treestats::tot_coph(phylo)},
                                 error = function(e) {return(NA)})

  stats$var_depth <- treestats::var_leaf_depth(phylo)
  stats$symmetry_nodes     <- treestats::sym_nodes(phylo)

  stats$mpd          <- tryCatch(expr = {treestats::mean_pair_dist(phylo)},
                                 error = function(e) {return(NA)})
  stats$psv          <- tryCatch(expr = {treestats::psv(phylo)},
                                 error = function(e) {return(NA)})
  stats$vpd          <- tryCatch(expr = {treestats::var_pair_dist(phylo)},
                                 error = function(e) {return(NA)})
  stats$mntd         <- tryCatch(expr = {treestats::mntd(phylo)},
                                 error = function(e) {return(NA)})
  stats$j_stat       <- tryCatch(expr = {treestats::entropy_j(phylo)},
                                 error = function(e) {return(NA)})

  stats$rquartet     <- treestats::rquartet(phylo)


  stats$wiener          <- treestats::wiener(phylo)
  stats$max_betweenness <- treestats::max_betweenness(phylo)
  stats$max_closeness   <- treestats::max_closeness(phylo, weight = TRUE)
  stats$diameter        <- treestats::diameter(phylo, weight = FALSE)
  stats$eigenvector     <- max(treestats::eigen_vector(phylo)$eigenvector)

  return(stats)
  }
