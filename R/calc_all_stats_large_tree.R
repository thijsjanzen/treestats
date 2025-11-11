#' Apply all available tree statistics to a single tree, but exclude those
#' statistics that require large amounts of memory.
#' @param phylo phylo object
#' @param normalize if set to TRUE, results are normalized (if possible) under
#' either the  Yule expectation (if available), or the number of tips
#' @return List with statistics
#' @export
#' @description this function applies all tree statistics available in
#' this package to a single tree, excluding statistics that require large
#' amounts of memory, because they are based on a distance matrix. The
#' statistics EXCLUDED are:
#' \itemize{
#'   \item laplacian spectrum
#'   \item variance of pairwise distance (vpd)
#'   \item maximum eigen vector value
#'   \item minimum eigenvalue of the Laplacian matrix
#'   \item maximum eigenvalue of the Laplacian matrix
#'   \item minimum eigenvalue of the adjacency matrix
#'   \item maximum eigenvalue of the adjacency matrix
#' }
#'
calc_all_stats_large_tree <- function(phylo, normalize = FALSE) {

  stats <- list()

  stats$gamma              <- try_stat(phylo, treestats::gamma_statistic)

  stats$sackin             <- try_stat(phylo, treestats::sackin,
                                       normalize, c("yule", "none"))

  stats$colless            <- try_stat(phylo, treestats::colless,
                                       normalize, c("yule", "none"))
  stats$colless_corr       <- try_stat(phylo, treestats::colless_corr,
                                       normalize, c("yule", "none"))
  stats$colless_quad       <- try_stat(phylo, treestats::colless_quad,
                                       normalize, c("yule", "none"))

  stats$beta               <- try_stat(phylo, treestats::beta_statistic)

  stats$blum               <- try_stat(phylo, treestats::blum,
                                       normalize, c(TRUE, FALSE))

  stats$crown_age          <- try_stat(phylo, treestats::crown_age)
  stats$tree_height        <- try_stat(phylo, treestats::tree_height)
  stats$pigot_rho          <- try_stat(phylo, treestats::pigot_rho)

  stats$number_of_lineages <- try_stat(phylo, treestats::number_of_lineages)
  stats$nltt_base          <- try_stat(phylo, treestats::nLTT_base)
  stats$phylogenetic_div   <- try_stat(phylo, treestats::phylogenetic_diversity)

  stats$avg_ladder         <- try_stat(phylo, treestats::avg_ladder)
  stats$max_ladder         <- try_stat(phylo, treestats::max_ladder)

  stats$cherries           <- try_stat(phylo, treestats::cherries,
                                       normalize, c("yule", "none"))

  stats$double_cherries    <- try_stat(phylo, treestats::double_cherries)
  stats$four_prong         <- try_stat(phylo, treestats::four_prong)

  stats$il_number          <- try_stat(phylo, treestats::ILnumber,
                                       normalize, c("tips", "none"))

  stats$pitchforks         <- try_stat(phylo, treestats::pitchforks,
                                       normalize, c("tips", "none"))

  stats$stairs              <- try_stat(phylo, treestats::stairs)

  stats$imbalance_steps  <- try_stat(phylo, treestats::imbalance_steps,
                                     normalize)

  stats$j_one              <- try_stat(phylo, treestats::j_one)

  stats$b1                 <- try_stat(phylo, treestats::b1,
                                       normalize, c("tips", "none"))

  stats$b2                 <- try_stat(phylo, treestats::b2,
                                       normalize, c("yule", "none"))

  stats$area_per_pair      <- try_stat(phylo, treestats::area_per_pair,
                                       normalize, c("yule", "none"))

  stats$average_leaf_depth <- try_stat(phylo, treestats::average_leaf_depth,
                                       normalize, c("yule", "none"))

  stats$i_stat             <- try_stat(phylo, treestats::mean_i)

  stats$ew_colless         <- try_stat(phylo, treestats::ew_colless)

  stats$max_del_width      <- try_stat(phylo, treestats::max_del_width,
                                       normalize, c("tips", "none"))

  stats$max_depth          <- try_stat(phylo, treestats::max_depth,
                                       normalize, c("tips", "none"))

  stats$avg_vert_depth     <- try_stat(phylo, treestats::avg_vert_depth)

  stats$max_width          <- try_stat(phylo, treestats::max_width,
                                       normalize, c("tips", "none"))

  stats$mw_over_md         <- try_stat(phylo, treestats::mw_over_md)

  stats$tot_path           <- try_stat(phylo, treestats::tot_path_length)

  stats$tot_internal_path  <- try_stat(phylo, treestats::tot_internal_path)

  stats$rogers             <- try_stat(phylo, treestats::rogers,
                                       normalize, c("tips", "none"))


  stats$stairs2            <- try_stat(phylo, treestats::stairs2)

  stats$tot_coph           <- try_stat(phylo, treestats::tot_coph,
                                       normalize, c("yule", "none"))

  stats$var_depth          <- try_stat(phylo, treestats::var_leaf_depth,
                                       normalize, c("yule", "none"))

  stats$symmetry_nodes     <- try_stat(phylo, treestats::sym_nodes,
                                       normalize, c("tips", "none"))

  stats$mpd                <- try_stat(phylo, treestats::mean_pair_dist,
                                       normalize, c("tips", "none"))

  stats$psv                <- try_stat(phylo, treestats::psv,
                                       normalize, c("tips", "none"))

  stats$vpd                <- try_stat(phylo, treestats::var_pair_dist)

  stats$mntd               <- try_stat(phylo, treestats::mntd)

  stats$j_stat             <- try_stat(phylo, treestats::entropy_j)

  stats$rquartet           <- try_stat(phylo, treestats::rquartet,
                                       normalize, c("yule", "none"))

  stats$wiener             <- try_stat(phylo, treestats::wiener, normalize)

  stats$max_betweenness    <- try_stat(phylo, treestats::max_betweenness,
                                       normalize, c("tips", "none"))

  local_closeness <- function(tree, w, n) {
    return(treestats::max_closeness(tree, weight = w,
                                    normalization =
                                      ifelse(n == TRUE, "tips", "none")))
  }

  stats$max_closeness      <- try_stat(phylo, function(x) {
    local_closeness(x, FALSE, normalize)})

  stats$max_closenessW     <- try_stat(phylo, function(x) {
    local_closeness(x, TRUE, normalize)})

  stats$diameter           <- try_stat(phylo, treestats::diameter)

  if (1 == 2) {
  stats$eigen_centrality   <- try_stat(phylo,
                                       function(x) {
                                         return(max(treestats::eigen_centrality(x,
                                                                                weight = FALSE)$eigenvector))}) #nolint

  stats$eigen_centralityW  <- try_stat(phylo,
                                       function(x) {
                                         return(max(treestats::eigen_centrality(x,
                                                                                weight = TRUE)$eigenvector))}) #nolint
}

  stats$mean_branch_length <- try_stat(phylo, treestats::mean_branch_length)
  stats$var_branch_length  <- try_stat(phylo, treestats::var_branch_length)

  stats$mean_branch_length_int <- try_stat(phylo,
                                           treestats::mean_branch_length_int)
  stats$mean_branch_length_ext <- try_stat(phylo,
                                           treestats::mean_branch_length_ext)
  stats$var_branch_length_int  <- try_stat(phylo,
                                           treestats::var_branch_length_int)
  stats$var_branch_length_ext  <- try_stat(phylo,
                                           treestats::var_branch_length_ext)

  stats$treeness <- try_stat(phylo, treestats::treeness)

  stats$root_imbalance <- try_stat(phylo, treestats::root_imbalance)

  return(stats)
}
