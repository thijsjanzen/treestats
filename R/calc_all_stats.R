#' Apply all available tree statistics to a single tree
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
#'   \item corrected Colless
#'   \item quadratic Colless
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
#'   \item double cherries
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
#'   \item treeness statistic
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

  stats$crown_age          <- treestats::crown_age(phylo)
  stats$tree_height        <- treestats::tree_height(phylo)
  stats$pigot_rho          <- treestats::pigot_rho(phylo)

  stats$number_of_lineages <- treestats::number_of_lineages(phylo)
  stats$nltt_base          <- treestats::nLTT_base(phylo)
  stats$phylogenetic_div   <- treestats::phylogenetic_diversity(phylo)

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

  temp_stats <- try_stat(phylo,
                      function(x) {
                        return(treestats::minmax_laplace(x, TRUE))
                      })

  if (length(temp_stats) >= 2) {
    stats$min_laplace <- temp_stats$min
    stats$max_laplace <- temp_stats$max
  } else {
    stats$min_laplace <- NA
    stats$max_laplace <- NA
  }

  temp_stats <- try_stat(phylo,
                         function(x) {
                           return(treestats::minmax_adj(x, TRUE))
                          })

  if (length(temp_stats) >= 2) {
    stats$min_adj <- temp_stats$min
    stats$max_adj <- temp_stats$max
  } else {
    stats$min_adj <- NA
    stats$max_adj <- NA
  }

  stats$imbalance_steps  <-
    treestats::imbalance_steps(phylo,
                               normalization = normalize)

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

  stats$eigen_centrality   <- try_stat(phylo,
                                       function(x) {
                          return(max(treestats::eigen_centrality(x,
                                       weight = FALSE)$eigenvector))}) #nolint

  stats$eigen_centralityW  <- try_stat(phylo,
                                       function(x) {
                                      return(max(treestats::eigen_centrality(x,
                                           weight = TRUE)$eigenvector))}) #nolint


  stats$mean_branch_length <- treestats::mean_branch_length(phylo)
  stats$var_branch_length  <- treestats::var_branch_length(phylo)

  stats$mean_branch_length_int <- treestats::mean_branch_length_int(phylo)
  stats$mean_branch_length_ext <- treestats::mean_branch_length_ext(phylo)
  stats$var_branch_length_int <- treestats::var_branch_length_int(phylo)
  stats$var_branch_length_ext <- treestats::var_branch_length_ext(phylo)

  stats$treeness <- treestats::treeness(phylo)

  stats$root_imbalance <- try_stat(phylo, treestats::root_imbalance)

  stats <- unlist(stats)
  stats <- stats[order(names(stats))]

  return(stats)
}

#' @keywords internal
try_stat <- function(phylo,
                     func,
                     normalize = FALSE,
                     norm_res = c(TRUE, FALSE)) {
  res <- NA
  if (normalize) {
    res <- tryCatch(expr = {
                            func(phylo,
                                 normalization =
                                   ifelse(normalize,
                                          norm_res[1], norm_res[2]))
                           },
                    error = function(e) {return(NA)})     #nolint
  } else {
    res <- tryCatch(expr = {func(phylo)},  #nolint
                    error = function(e) {return(NA)})     #nolint
  }
  return(res)
}
