context("blum")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  blum1 <- treestats::blum(focal_tree)
  blum_check <- castor::tree_imbalance(focal_tree, type = "Blum")
  testthat::expect_equal(blum1, blum_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  blum2 <- treestats::blum(ltab)
  # testthat::expect_equal(blum1, blum2)  # FAIL TODO: fix!
})

test_that("abuse", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.3)[[1]]
  testthat::expect_error(
    treestats::blum(focal_tree),
    "can only calculate statistic for ultrametric tree"
  )
})
