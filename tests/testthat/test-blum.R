context("blum")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  blum1 <- treestats::blum(focal_tree)
  blum_check <- castor::tree_imbalance(focal_tree, type = "Blum")
  testthat::expect_equal(blum1, blum_check)
})
