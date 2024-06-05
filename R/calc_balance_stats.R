#' Apply all balance statistics to a single tree
#' @param phylo phylo object
#' @param normalize if set to TRUE, results are normalized (if possible) under
#' either the  Yule expectation (if available), or the number of tips
#' @return list with statistics
#' @export
#' @description this function applies all tree statistics available in
#' this package to a single tree, being:
#' \itemize{
#'   \item{Sackin}
#'   \item{Colless}
#'   \item{Aldous' beta statistic}
#'   \item{Blum}
#'   \item{Average Ladder Size}
#'   \item{cherries}
#'   \item{ILnumber}
#'   \item{pitchforks}
#'   \item{stairs}
#'   \item{stairs2}
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
#'   \item{rquartet index}
#'   \item{j_one statistic}
#'   \item{diameter}
#' }
#'
calc_balance_stats <- function(phylo, normalize = FALSE) {

  stats <- list()

  stats$sackin             <- try_stat(phylo, treestats::sackin,
                                       normalize, c("yule", "none"))
  stats$colless            <- try_stat(phylo, treestats::colless,
                                       normalize, c("yule", "none"))

  stats$beta               <- try_stat(phylo, treestats::beta_statistic)

  stats$blum               <- try_stat(phylo, treestats::blum, normalize)

  stats$avg_ladder         <- try_stat(phylo, treestats::avg_ladder)
  stats$max_ladder         <- try_stat(phylo, treestats::max_ladder)
  stats$cherries           <- try_stat(phylo, treestats::cherries,
                                       normalize, c("yule", "none"))

  stats$il_number          <- try_stat(phylo, treestats::ILnumber,
                                       normalize, c("tips", "none"))

  stats$pitchforks         <- try_stat(phylo, treestats::pitchforks,
                                       normalize, c("tips", "none"))

  stats$stairs              <- try_stat(phylo, treestats::stairs)

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

  stats$rogers             <- try_stat(phylo, treestats::rogers,
                                       normalize, c("tips", "none"))


  stats$stairs2            <- try_stat(phylo, treestats::stairs2)

  stats$tot_coph           <- try_stat(phylo, treestats::tot_coph,
                                       normalize, c("yule", "none"))

  stats$var_depth          <- try_stat(phylo, treestats::var_depth,
                                       normalize, c("yule", "none"))

  stats$symmetry_nodes     <- try_stat(phylo, treestats::sym_nodes,
                                       normalize, c("tips", "none"))

  stats$rquartet           <- try_stat(phylo, treestats::rquartet,
                                       normalize, c("yule", "none"))

  stats$imbalance_steps    <- treestats::imbalance_steps(phylo,
                                                      normalization = normalize)

  stats$j_one              <- try_stat(phylo, treestats::j_one)

  stats$diameter           <- try_stat(phylo, treestats::diameter)

  return(stats)
}
