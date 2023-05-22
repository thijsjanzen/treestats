#' function to apply all balance statistics to a single tree
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
#' }
#'
calc_balance_stats <- function(phylo, normalize = FALSE) {

  stats <- list()

  stats$sackin             <- treestats::sackin(phylo,
                                                normalization =
                                                  ifelse(normalize,
                                                         "yule", "none"))
  stats$colless            <- treestats::colless(phylo,
                                                 normalization =
                                                   ifelse(normalize,
                                                          "yule", "none"))

  stats$beta          <- tryCatch(expr = {treestats::beta_statistic(phylo) }, #nolint
                                  error = function(e) {return(NA) }) #nolint

  stats$blum               <- treestats::blum(phylo,
                                              normalization = normalize)

  stats$avg_ladder          <- treestats::avg_ladder(phylo) #nolint
  stats$max_ladder          <- treestats::max_ladder(phylo) #nolint
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

  stats$average_leaf_depth  <-
        treestats::average_leaf_depth(phylo,
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
  stats$tot_coph      <- treestats::tot_coph(phylo,
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

  stats$rquartet     <- treestats::rquartet(phylo,
                                            normalization =
                                              ifelse(normalize,
                                                     "yule", "none"))

  stats$imbalance_steps <- treestats::imbalance_steps(phylo,
                                                      normalization = normalize)

  stats$j_one       <- treestats::j_one(phylo)

  stats$diameter    <- treestats::diameter(phylo)

  return(stats)
}
