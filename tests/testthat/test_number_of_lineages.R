context("number_of_lineages")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  num_lin <- treestats::number_of_lineages(focal_tree)

  testthat::expect_equal(num_lin, 100)

  testthat::expect_error(
    treestats::number_of_lineages(num_lin),
    "object \"phy\" is not of class \"phylo\""
  )
})
