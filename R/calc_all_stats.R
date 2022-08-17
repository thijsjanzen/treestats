#' function to apply all available tree statistics to a single tree
#' @param phylo phylo object
#' @param normalize if set to TRUE, results are normalized (if possible) under
#' either the  Yule expectation (if available), or the number of tips
#' @return list with statistics
#' @export
#' @description this function applies all tree statistics available in
#' this package to a
#' single tree, being:
#' \itemize{
#'   \item{gamma}
#'   \item{Sackin}
#'   \item{Colless}
#'   \item{Aldous' beta statistic}
#'   \item{Blum}
#'   \item{crown age}
#'   \item{tree height}
#'   \item{Pigot's rho}
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
calc_all_stats <- function(phylo, normalize = FALSE) {

  stats <- list()

  stats$gamma              <- treestats::gamma_statistic(phylo)
  stats$sackin             <- treestats::sackin(phylo,
                                                normalization =
                                                  ifelse(normalize,
                                                         "yule", "none"))
  stats$colless            <- treestats::colless(phylo,
                                                 normalization =
                                                   ifelse(normalize,
                                                          "yule", "none"))

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
  stats$cherries           <- treestats::cherries(phylo,
                                                  normalization =
                                                    ifelse(normalize,
                                                           "yule", "none"))

  stats$il_number           <- treestats::ILnumber(phylo,
                                                   normalization =
                                                     ifelse(normalize,
                                                            "tips", "none"))

  stats$pitchforks         <- treestats::pitchforks(phylo,
                                                    normalization =
                                                      ifelse(normalize,
                                                             "tips", "none"))
  stats$stairs             <- treestats::stairs(phylo)

  temp_stats <- tryCatch(expr = {treestats::calc_lapl_spectrum(phylo) }, #nolint
                         error = function(e) {return(NA) }) #nolint

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


  stats$b1           <- treestats::b1(phylo,
                                      normalization =
                                        ifelse(normalize,
                                               "tips", "none"))

  stats$b2           <- treestats::b2(phylo,
                                      normalization =
                                        ifelse(normalize,
                                               "yule", "none"))

  stats$area_per_pair <- treestats::area_per_pair(phylo,
                                                  normalization =
                                                    ifelse(normalize,
                                                           "yule", "none"))

  stats$average_leaf_depth  <- treestats::average_leaf_depth(phylo,
                                      normalization = ifelse(normalize,
                                                             "yule", "none"))

  stats$i_stat       <- treestats::mean_i(phylo)
  stats$ew_colless    <- treestats::ew_colless(phylo)
  stats$max_del_width <- treestats::max_del_width(phylo,
                                                  normalization =
                                                   ifelse(normalize,
                                                          "tips", "none"))


  stats$max_depth     <- treestats::max_depth(phylo,
                                              normalization =
                                                ifelse(normalize,
                                                       "tips", "none"))
  stats$max_width     <- treestats::max_width(phylo,
                                              normalization =
                                                ifelse(normalize,
                                                       "tips", "none"))
  stats$rogers       <- treestats::rogers(phylo,
                                          normalization =
                                            ifelse(normalize,
                                                   "tips", "none"))
  stats$stairs2      <- treestats::stairs2(phylo)
  stats$tot_coph      <- tryCatch(expr = {
                          treestats::tot_coph(phylo,
                                              normalization =
                                                ifelse(normalize,
                                                       "yule", "none"))},
                                  error = function(e) {return(NA)}) #nolint

  stats$var_depth <- treestats::var_leaf_depth(phylo,
                                               normalization =
                                                 ifelse(normalize,
                                                        "yule", "none"))
  stats$symmetry_nodes  <- treestats::sym_nodes(phylo,
                                                normalization =
                                                  ifelse(normalize,
                                                    "tips", "none"))

  stats$mpd          <- tryCatch(expr = {
                          treestats::mean_pair_dist(phylo,
                                                    normalization =
                                                     ifelse(normalize,
                                                            "tips", "none"))},
                                 error = function(e) {return(NA)}) #nolint
  stats$psv          <- tryCatch(expr = {treestats::psv(phylo, #nolint
                                                        normalization =
                                                        ifelse(normalize,
                                                              "tips", "none"))},
                                 error = function(e) {return(NA)}) #nolint
  stats$vpd          <- tryCatch(expr = {treestats::var_pair_dist(phylo)},  #nolint
                                 error = function(e) {return(NA)})     #nolint
  stats$mntd         <- tryCatch(expr = {treestats::mntd(phylo)},      #nolint
                                 error = function(e) {return(NA)})     #nolint
  stats$j_stat       <- tryCatch(expr = {treestats::entropy_j(phylo)}, #nolint
                                 error = function(e) {return(NA)})     #nolint

  stats$rquartet     <- treestats::rquartet(phylo,
                                            normalization =
                                              ifelse(normalize,
                                                     "yule", "none"))


  stats$wiener          <- treestats::wiener(phylo, normalize = normalize)
  stats$max_betweenness <- treestats::max_betweenness(phylo,
                                                      normalization =
                                                        ifelse(normalize,
                                                               "tips", "none"))
  stats$max_closeness   <- treestats::max_closeness(phylo,
                                                    normalization =
                                                      ifelse(normalize,
                                                             "tips", "none"))

  stats$diameter        <- treestats::diameter(phylo,
                                               normalization =
                                                 ifelse(normalize,
                                                        "minmax", "none"))
  stats$eigenvector     <- max(treestats::eigen_vector(phylo)$eigenvector)

  return(stats)
}
