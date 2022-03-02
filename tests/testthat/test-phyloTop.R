context("phyloTop")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  avglad <- treestats::avgLadder(focal_tree)
  avglad_check <- phyloTop::avgLadder(focal_tree)
  testthat::expect_equal(avglad, avglad_check)

  strs <- treestats::stairs(focal_tree)
  strs_check <- phyloTop::stairs(focal_tree)[[1]]
  testthat::expect_equal(strs, strs_check)
})
