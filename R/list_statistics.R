#' Provides a list of all available statistics in the package
#' @param only_balance_stats only return those statistics associated with
#' measuring balance of a tree
#' @return vector with names of summary statistics
#' @export
list_statistics <- function(only_balance_stats = FALSE) {
  all_statistics <- c("gamma", "sackin",
                      "colless", "colless_corr", "colless_quad",
                      "beta", "blum", "crown_age",
                      "tree_height", "pigot_rho", "number_of_lineages",
                      "nltt_base", "phylogenetic_div", "avg_ladder",
                      "max_ladder", "cherries", "double_cherries",
                      "il_number", "pitchforks", "four_prong",
                      "stairs", "laplace_spectrum_a", "laplace_spectrum_p",
                      "laplace_spectrum_e", "laplace_spectrum_g",
                      "min_laplace", "max_laplace",
                      "min_adj", "max_adj",
                      "imbalance_steps", "j_one", "b1", "b2", "area_per_pair",
                      "average_leaf_depth", "avg_vert_depth", "i_stat",
                      "ew_colless",
                      "max_del_width", "max_depth", "max_width", "mw_over_md",
                      "rogers", "stairs2", "tot_coph", "var_depth",
                      "symmetry_nodes",
                      "mpd", "psv", "vpd", "mntd", "j_stat", "rquartet",
                      "wiener", "max_betweenness",
                      "max_closeness", "max_closenessW",
                      "diameter",
                      "eigen_centrality", "eigen_centralityW",
                      "tot_path", "tot_internal_path",
                      "mean_branch_length", "var_branch_length",
                      "mean_branch_length_int", "mean_branch_length_ext",
                      "var_branch_length_int", "var_branch_length_ext",
                      "treeness", "root_imbalance")

  all_statistics <- sort(all_statistics)

  bal_stats <- c("sackin", "colless", "beta", "blum", "avg_ladder",
                 "max_ladder", "cherries", "il_number", "pitchforks", "stairs",
                 "b1", "b2", "area_per_pair", "average_leaf_depth", "i_stat",
                 "ew_colless", "max_del_width", "max_depth", "max_width",
                 "rogers", "stairs2", "tot_coph", "var_depth", "symmetry_nodes",
                 "rquartet", "imbalance_steps", "j_one", "diameter")

  bal_stats <- sort(bal_stats)

  if (only_balance_stats) return(bal_stats)

  return(all_statistics)

}
