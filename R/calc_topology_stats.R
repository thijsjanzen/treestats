#' Calculate all topology based statistics for a single tree
#' @param phylo phylo object
#' @param normalize if set to TRUE, results are normalized (if possible) under
#' either the  Yule expectation (if available), or the number of tips
#' @return list with statistics
#' @export
#' @description this function calculates all tree statistics based on topology
#' available in this package for a single tree, being:
#' \itemize{
#'   \item{area_per_pair}
#'   \item{average_leaf_depth}
#'   \item{avg_ladder}
#'   \item{avg_vert_depth}
#'   \item{b1}
#'   \item{b2}
#'   \item{beta}
#'   \item{blum}
#'   \item{cherries}
#'   \item{colless}
#'   \item{colless_corr}
#'   \item{colless_quad}
#'   \item{diameter}
#'   \item{double_cherries}
#'   \item{eigen_centrality}
#'   \item{ew_colless}
#'   \item{four_prong}
#'   \item{i_stat}
#'   \item{il_number}
#'   \item{imbalance_steps}
#'   \item{j_one}
#'   \item{max_betweenness}
#'   \item{max_closeness}
#'   \item{max_del_width}
#'   \item{max_depth}
#'   \item{max_ladder}
#'   \item{max_width}
#'   \item{mw_over_md}
#'   \item{pitchforks}
#'   \item{rogers}
#'   \item{root_imbalance}
#'   \item{rquartet}
#'   \item{sackin}
#'   \item{stairs}
#'   \item{stairs2}
#'   \item{symmetry_nodes}
#'   \item{tot_coph}
#'   \item{tot_internal_path}
#'   \item{tot_path_length}
#'   \item{var_depth}
#' }
#'
calc_topology_stats <- function(phylo, normalize = FALSE) {

  stats <- list()

  stats$sackin             <- try_stat(phylo, treestats::sackin,
                                       normalize, c("yule", "none"))
  stats$colless            <- try_stat(phylo, treestats::colless,
                                       normalize, c("yule", "none"))

  stats$colless_corr       <- try_stat(phylo, treestats::colless_corr)

  stats$colless_quad       <- try_stat(phylo, treestats::colless_quad)

  stats$beta               <- try_stat(phylo, treestats::beta_statistic)

  stats$blum               <- try_stat(phylo, treestats::blum, normalize)

  stats$avg_ladder         <- try_stat(phylo, treestats::avg_ladder)
  stats$max_ladder         <- try_stat(phylo, treestats::max_ladder)
  stats$cherries           <- try_stat(phylo, treestats::cherries,
                                       normalize, c("yule", "none"))

  stats$double_cherries    <- try_stat(phylo, treestats::double_cherries)

  stats$il_number          <- try_stat(phylo, treestats::ILnumber,
                                       normalize, c("tips", "none"))

  stats$pitchforks         <- try_stat(phylo, treestats::pitchforks,
                                       normalize, c("tips", "none"))

  stats$four_prong         <- try_stat(phylo, treestats::four_prong)

  stats$stairs             <- try_stat(phylo, treestats::stairs)

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

  stats$max_width          <- try_stat(phylo, treestats::max_width,
                                       normalize, c("tips", "none"))

  stats$mw_over_md         <- try_stat(phylo, treestats::mw_over_md)

  stats$rogers             <- try_stat(phylo, treestats::rogers,
                                       normalize, c("tips", "none"))


  stats$stairs2            <- try_stat(phylo, treestats::stairs2)

  stats$tot_coph           <- try_stat(phylo, treestats::tot_coph,
                                       normalize, c("yule", "none"))

  stats$var_depth          <- try_stat(phylo, treestats::var_leaf_depth,
                                       normalize, c("yule", "none"))

  stats$symmetry_nodes     <- try_stat(phylo, treestats::sym_nodes,
                                       normalize, c("tips", "none"))

  stats$rquartet           <- try_stat(phylo, treestats::rquartet,
                                       normalize, c("yule", "none"))

  stats$imbalance_steps    <- treestats::imbalance_steps(phylo,
                                                      normalization = normalize)

  stats$j_one              <- try_stat(phylo, treestats::j_one)

  stats$diameter           <- try_stat(phylo, treestats::diameter)

  stats$avg_vert_depth     <- try_stat(phylo, treestats::avg_vert_depth)

  calc_local_eigen_centrality <- function(x) {
    res <- treestats::eigen_centrality(x, weight = FALSE, scale = FALSE)
    return(max(res$eigenvector))
  }

  stats$eigen_centrality   <- try_stat(phylo, calc_local_eigen_centrality)

  stats$max_betweenness    <- try_stat(phylo, treestats::max_betweenness)

  calc_max_closeness <- function(x, n) {
    return(treestats::max_closeness(x, weight = FALSE, normalization = n))
  }

  stats$max_closeness      <- try_stat(phylo, calc_max_closeness,
                                     normalize, c("tips", "none"))

  stats$root_imbalance     <- try_stat(phylo, treestats::root_imbalance)

  stats$tot_internal_path  <- try_stat(phylo, treestats::tot_internal_path)

  stats$tot_path_length    <- try_stat(phylo, treestats::tot_path_length)

  stats <- unlist(stats)
  stats <- stats[order(names(stats))]

  return(stats)
}
