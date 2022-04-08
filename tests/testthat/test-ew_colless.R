context("ew_colless")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  colless <- treestats::ew_colless(focal_tree)
  colless_check <- treebalance::ewCollessI(focal_tree)
  testthat::expect_equal(colless, colless_check)

  focal_ltab <- treestats::phylo_to_l(focal_tree)

  colless <- treestats::ew_colless(focal_ltab)
  testthat::expect_equal(colless, colless_check)

  ## with extinct lineages:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.2)[[1]]

  colless <- treestats::ew_colless(focal_tree)
  colless_check <- treebalance::ewCollessI(focal_tree)
  testthat::expect_equal(colless, colless_check)

  focal_ltab <- treestats::phylo_to_l(focal_tree)

  colless <- treestats::ew_colless(focal_ltab)
  testthat::expect_equal(colless, colless_check)
})
