#' function to apply all available tree statistics to a single tree
#' @param phylo phylo object
#' @param normalize if set to TRUE, results are normalized (if possible) under
#' either the  Yule expectation (if available), or the number of tips
#' @return List with statistics
#' @export
#' @description this function applies all tree statistics available in
#' this package to a
#' single tree, being:
#' \itemize{
#'   \item gamma
#'   \item Sackin
#'   \item Colless
#'   \item Aldous' beta statistic
#'   \item Blum
#'   \item crown age
#'   \item tree height
#'   \item Pigot's rho
#'   \item number of lineages
#'   \item nLTT with empty tree
#'   \item phylogenetic diversity
#'   \item avgLadder index
#'   \item cherries
#'   \item ILnumber
#'   \item pitchforks
#'   \item stairs
#'   \item stairs2
#'   \item laplacian spectrum
#'   \item B1
#'   \item B2
#'   \item area per pair (aPP)
#'   \item average leaf depth (aLD)
#'   \item I statistic
#'   \item ewColless
#'   \item max Delta Width (maxDelW)
#'   \item maximum of Depth
#'   \item variance of Depth
#'   \item maximum Width
#'   \item Rogers
#'   \item total Cophenetic distance
#'   \item symmetry Nodes
#'   \item mean of pairwise distance (mpd)
#'   \item variance of pairwise distance (vpd)
#'   \item Phylogenetic Species Variability (psv)
#'   \item mean nearest taxon distance (mntd)
#'   \item J statistic of entropy
#'   \item rquartet index
#'   \item Wiener index
#'   \item max betweenness
#'   \item max closeness
#'   \item diameter, without branch lenghts
#'   \item maximum eigen vector value
#'   \item mean branch length
#'   \item variance of branch length
#'   \item mean external branch length
#'   \item variance of external branch length
#'   \item mean internal branch length
#'   \item variance of internal branch length
#'   \item number of imbalancing steps
#'   \item j_one statistic
#' }
#'
#' For the Laplacian spectrum properties, four properties of the eigenvalue
#' distribution are returned: 1) asymmetry, 2) peakedness, 3) log(principal
#' eigenvalue) and 4) eigengap.
#' Please notice that for some very small or very large trees, some of the
#' statistics can not be calculated. The function will report an NA for this
#' statistic, but will not break, to facilitate batch analysis of large numbers
#' of trees.
#' @rawNamespace import(Rcpp)
#' @rawNamespace import(nloptr)
#' @rawNamespace useDynLib(treestats)
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

  stats$beta          <- tryCatch(expr = {treestats::beta_statistic(phylo)}, #nolint
                                  error = function(e) {return(NA) }) #nolint

  stats$blum               <- treestats::blum(phylo, normalization = normalize)
  stats$crown_age          <- treestats::crown_age(phylo)
  stats$tree_height        <- treestats::tree_height(phylo)
  stats$pigot_rho          <- treestats::pigot_rho(phylo)

  stats$number_of_lineages <- treestats::number_of_lineages(phylo)
  stats$nltt_base          <- treestats::nLTT_base(phylo)
  stats$phylogenetic_div   <- treestats::phylogenetic_diversity(phylo)
  stats$avg_ladder         <- treestats::avg_ladder(phylo) #nolint
  stats$max_ladder         <- treestats::max_ladder(phylo)
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

  temp_stats <- tryCatch(expr = {treestats::laplacian_spectrum(phylo) }, #nolint
                         error = function(e) {return(NA) }) #nolint

  if (length(temp_stats) == 5) {
    stats$laplace_spectrum_a  <- temp_stats$asymmetry
    stats$laplace_spectrum_p  <- temp_stats$peakedness
    stats$laplace_spectrum_e  <- log(temp_stats$principal_eigenvalue)
    stats$laplace_spectrum_g  <- temp_stats$eigengap[[1]]
  } else {
    stats$laplace_spectrum_a  <- NA
    stats$laplace_spectrum_p  <- NA
    stats$laplace_spectrum_e  <- NA
    stats$laplace_spectrum_g  <- NA
  }

  stats$imbalance_steps  <-
    treestats::imbalance_steps(phylo,
                               normalization = normalize)

  stats$j_one        <- treestats::j_one(phylo)

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
  stats$tot_coph      <-  treestats::tot_coph(phylo,
                                              normalization =
                                                ifelse(normalize,
                                                       "yule", "none"))

  stats$var_depth <- treestats::var_leaf_depth(phylo,
                                               normalization =
                                                 ifelse(normalize,
                                                        "yule", "none"))
  stats$symmetry_nodes  <- treestats::sym_nodes(phylo,
                                                normalization =
                                                  ifelse(normalize,
                                                    "tips", "none"))

  stats$mpd          <- treestats::mean_pair_dist(phylo,
                                                    normalization =
                                                     ifelse(normalize,
                                                            "tips", "none"))

  stats$psv          <- treestats::psv(phylo, #nolint
                                                        normalization =
                                                        ifelse(normalize,
                                                              "tips", "none"))
  stats$vpd          <- tryCatch(expr = {treestats::var_pair_dist(phylo)},  #nolint
                                 error = function(e) {return(NA)})     #nolint
  stats$mntd         <- treestats::mntd(phylo)

  stats$j_stat       <- treestats::entropy_j(phylo)

  stats$rquartet     <- treestats::rquartet(phylo,
                                            normalization =
                                              ifelse(normalize,
                                                     "yule", "none"))


  stats$wiener          <- treestats::wiener(phylo, normalization = normalize)
  stats$max_betweenness <- treestats::max_betweenness(phylo,
                                                      normalization =
                                                        ifelse(normalize,
                                                               "tips", "none"))
  stats$max_closeness   <- treestats::max_closeness(phylo,
                                                    normalization =
                                                      ifelse(normalize,
                                                             "tips", "none"))

  stats$diameter        <- treestats::diameter(phylo)
  stats$eigenvector     <- max(treestats::eigen_vector(phylo)$eigenvector)


  stats$mean_branch_length <- treestats::mean_branch_length(phylo)
  stats$var_branch_length  <- treestats::var_branch_length(phylo)

  stats$mean_branch_length_int <- treestats::mean_branch_length_int(phylo)
  stats$mean_branch_length_ext <- treestats::mean_branch_length_ext(phylo)
  stats$var_branch_length_int <- treestats::var_branch_length_int(phylo)
  stats$var_branch_length_ext <- treestats::var_branch_length_ext(phylo)

  return(stats)
}
