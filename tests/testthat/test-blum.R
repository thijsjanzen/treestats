context("blum")

test_that("usage", {
  set.seed(42)
  castor_version <- packageVersion("castor")
  if (castor_version <= "1.7.2") {
    testthat::skip("old castor has wrong blum calculation")
    testthat::skip_on_cran("old castor has wrong blum calculation")
    testthat::skip_on_ci("old castor has wrong blum calculation")

  }
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  blum1 <- treestats::blum(focal_tree)
  blum_check <- castor::tree_imbalance(focal_tree, type = "Blum")
  testthat::expect_equal(blum1, blum_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  blum2 <- treestats::blum(ltab)
  testthat::expect_equal(blum1, blum2)
})

test_that("usage", {
  set.seed(5)
  focal_tree <- TreeSim::sim.bd.taxa(n = 4,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  blum1 <- treestats::blum(focal_tree)
  blum_check <- log(2 - 1) + log(3 - 1) + log(4 - 1)
  testthat::expect_equal(blum1, blum_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  blum2 <- treestats::blum(ltab)
  testthat::expect_equal(blum1, blum2)
})
