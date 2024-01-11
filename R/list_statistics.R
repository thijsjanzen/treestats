#' Provides a list of all available statistics in the package
#' @param only_balance_stats only return those statistics associated with
#' measuring balance of a tree
#' @return vector with names of summary statistics
#' @export
list_statistics <- function(only_balance_stats = FALSE) {
  all_statistics <- c("gamma", "sackin", "colless", "beta", "blum", "crown_age",
                      "tree_height", "pigot_rho", "number_of_lineages",
                      "nltt_base", "phylogenetic_div", "avg_ladder",
                      "max_ladder", "cherries", "il_number", "pitchforks",
                      "stairs", "laplace_spectrum_a", "laplace_spectrum_p",
                      "laplace_spectrum_e", "laplace_spectrum_g",
                      "imbalance_steps", "j_one", "b1", "b2", "area_per_pair",
                      "average_leaf_depth", "i_stat", "ew_colless",
                      "max_del_width", "max_depth", "max_width", "rogers",
                      "stairs2", "tot_coph", "var_depth", "symmetry_nodes",
                      "mpd", "psv", "vpd", "mntd", "j_stat", "rquartet",
                      "wiener", "max_betweenness", "max_closeness", "diameter",
                      "eigenvector", "mean_branch_length", "var_branch_length",
                      "mean_branch_length_int", "mean_branch_length_ext",
                      "var_branch_length_int", "var_branch_length_ext")

  bal_stats <- c("sackin", "colless", "beta", "blum", "avg_ladder",
                 "max_ladder", "cherries", "il_number", "pitchforks", "stairs",
                 "b1", "b2", "area_per_pair", "average_leaf_depth", "i_stat",
                 "ew_colless", "max_del_width", "max_depth", "max_width",
                 "rogers", "stairs2", "tot_coph", "var_depth", "symmetry_nodes",
                 "rquartet", "imbalance_steps", "j_one", "diameter")

  if (only_balance_stats) return(bal_stats)

  return(all_statistics)

}
