context("polytomies")

test_that("usage", {

  # Thanks to Fien Strijthaegen for providing a testing tree
  poly_tree <- ape::read.tree(text = "(151:34,(152:2,((((153:10,154:16):1,155:24):2,156:2):6,(157:5,((158:0,(((159:0,160:5):14,161:1):0,((((162:4,163:0):6,164:0):0,(165:22,166:4):3):1,167:13):4):9):3,(168:25,169:3):29):2):0):2):0,170:11,171:4):0;") #nolint

  all_stats <- treestats::calc_all_stats(poly_tree)

  test_na <- function(stat_name, local_stats) {
    testthat::expect_true(is.na(local_stats[names(local_stats) == stat_name]))
  }

  test_na("gamma", all_stats)
  test_na("colless", all_stats)
  test_na("beta", all_stats)
  test_na("avg_ladder", all_stats)
  test_na("max_ladder", all_stats)
  test_na("stairs", all_stats)
  test_na("j_one", all_stats)
  test_na("b2", all_stats)
  test_na("i_stat", all_stats)
  test_na("ew_colless", all_stats)
  test_na("rogers", all_stats)
  test_na("stairs2", all_stats)
  test_na("symmetry_nodes", all_stats)
  test_na("max_betweenness", all_stats)
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

  test_na("colless", bal_stats)
  test_na("beta", bal_stats)
  test_na("avg_ladder", bal_stats)
  test_na("max_ladder", bal_stats)
  test_na("cherries", bal_stats)
  test_na("stairs", bal_stats)
  test_na("j_one", bal_stats)
  test_na("b2", bal_stats)
  test_na("i_stat", bal_stats)
  test_na("ew_colless", bal_stats)
  test_na("rogers", bal_stats)
  test_na("stairs2", bal_stats)
  test_na("symmetry_nodes", bal_stats)
  test_na("diameter", bal_stats)

  # test phylo_to_l and rebase_ltable on polytomies:
  testthat::expect_silent(
    ltab  <- treestats::phylo_to_l(poly_tree))
  testthat::expect_silent(
    ltab2 <- treestats::rebase_ltable(ltab)
  )
})

test_that("manual trees", {
  # these trees were provided by Luise Haeuser.

  poly_tree <- ape::read.tree(text = "(((t9:1,t14:1)1:1,t4:1,t12:1)1:0.5,(((t6:1,t19:1,t10:1,t2:1)1:1,t8:1,(t5:1,t16:1)1:1)1:1,(((((t11:1,t17:1)1:1,(t7:1,t3:1)1:1)1:1,(t1:1,t18:1)1:1)1:1,t15:1)1:1,t13:1)1:1)1:0.5);") #nolint
  il_val <- treestats::ILnumber(poly_tree)
  testthat::expect_equal(il_val, 3) # reference value through solving by eye

  b1_val <- treestats::b1(poly_tree)
  b1_ref <- treebalance::B1I(poly_tree)
  testthat::expect_equal(b1_val, b1_ref)

  poly_tree <- ape::read.tree(text = "((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);") #nolint
  tip1 <- treestats::tot_internal_path(poly_tree)
  tip2 <- treebalance::totIntPathLen(poly_tree)
  testthat::expect_equal(tip1, tip2)

  b1_val <- treestats::b1(poly_tree)
  b1_ref <- treebalance::B1I(poly_tree)
  testthat::expect_equal(b1_val, b1_ref)
})
