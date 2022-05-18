context("all_Stats")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  testthat::expect_invisible(
   all_stats <- treestats::calc_all_stats(focal_tree)
  )
  testthat::expect_equal(length(all_stats), 39)
})
