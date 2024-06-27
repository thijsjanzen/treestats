context("polytomies")

test_that("usage", {

  # Thanks to Fien Strijthaegen for providing a testing tree
  poly_tree <- ape::read.tree(text = "(151:34,(152:2,((((153:10,154:16):1,155:24):2,156:2):6,(157:5,((158:0,(((159:0,160:5):14,161:1):0,((((162:4,163:0):6,164:0):0,(165:22,166:4):3):1,167:13):4):9):3,(168:25,169:3):29):2):0):2):0,170:11,171:4):0;") #nolint

  all_stats <- treestats::calc_all_stats(poly_tree)

  test_na <- function(stat_name, local_stats) {
    testthat::expect_true(is.na(local_stats[names(local_stats) == stat_name]))
  }

  test_na("gamma", all_stats)
  test_na("sackin", all_stats)
  test_na("colless", all_stats)
  test_na("beta", all_stats)
  test_na("blum", all_stats)
  test_na("avg_ladder", all_stats)
  test_na("max_ladder", all_stats)
  test_na("cherries", all_stats)
  test_na("stairs", all_stats)
  test_na("j_one", all_stats)
  test_na("b1", all_stats)
  test_na("b2", all_stats)
  test_na("area_per_pair", all_stats)
  test_na("average_leaf_depth", all_stats)
  test_na("i_stat", all_stats)
  test_na("ew_colless", all_stats)
  test_na("rogers", all_stats)
  test_na("stairs2", all_stats)
  test_na("tot_coph", all_stats)
  test_na("symmetry_nodes", all_stats)
  test_na("rquartet", all_stats)
  test_na("wiener", all_stats)
  test_na("max_betweenness", all_stats)
  test_na("max_closeness", all_stats)
  test_na("max_closenessW", all_stats)
  test_na("diameter", all_stats)
  test_na("eigen_centrality", all_stats)
  test_na("eigen_centralityW", all_stats)
  test_na("min_laplace", all_stats)
  test_na("max_laplace", all_stats)
  test_na("min_adj", all_stats)
  test_na("max_adj", all_stats)

  brts_stats <- treestats::calc_brts_stats(poly_tree)
  test_na("gamma", brts_stats)

  bal_stats <- treestats::calc_topology_stats(poly_tree)


  test_na("sackin", bal_stats)
  test_na("colless", bal_stats)
  test_na("beta", bal_stats)
  test_na("blum", bal_stats)
  test_na("avg_ladder", bal_stats)
  test_na("max_ladder", bal_stats)
  test_na("cherries", bal_stats)
  test_na("stairs", bal_stats)
  test_na("j_one", bal_stats)
  test_na("b1", bal_stats)
  test_na("b2", bal_stats)
  test_na("area_per_pair", bal_stats)
  test_na("average_leaf_depth", bal_stats)
  test_na("i_stat", bal_stats)
  test_na("ew_colless", bal_stats)
  test_na("rogers", bal_stats)
  test_na("stairs2", bal_stats)
  test_na("tot_coph", bal_stats)
  test_na("symmetry_nodes", bal_stats)
  test_na("rquartet", bal_stats)
  test_na("diameter", bal_stats)

  # test phylo_to_l and rebase_ltable on polytomies:
  testthat::expect_silent(
    ltab  <- treestats::phylo_to_l(poly_tree))
  testthat::expect_silent(
    ltab2 <- treestats::rebase_ltable(ltab)
  )
})
