context("average leaf depth")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  ald <- treestats::average_leaf_depth(focal_tree)
  ald_check <- treebalance::avgLeafDepI(focal_tree)
  testthat::expect_equal(ald, ald_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::average_leaf_depth(focal_tree),
                         treestats::average_leaf_depth(ltab))


  # with extinct species:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.2)[[1]]

  ald <- treestats::average_leaf_depth(focal_tree)
  ald_check <- treebalance::avgLeafDepI(focal_tree)
  testthat::expect_equal(ald, ald_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::average_leaf_depth(focal_tree),
                         treestats::average_leaf_depth(ltab))
})
