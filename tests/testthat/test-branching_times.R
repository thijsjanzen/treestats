context("branching_times")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 500,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  a1 <- ape::branching.times(focal_tree)
  a2 <- treestats::branching_times(focal_tree)

  diff_brts <- abs(sort(a1) - sort(a2))
  mean_diff <- mean(diff_brts)
  testthat::expect_lt(mean_diff, 1e-6)

  # again, but with extinct lineages:
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.3)[[1]]

  a1 <- ape::branching.times(focal_tree)
  a2 <- treestats::branching_times(focal_tree)

  diff_brts <- abs(sort(a1) - sort(a2))
  mean_diff <- mean(diff_brts)
  testthat::expect_lt(mean_diff, 1e-6)

  # using ltable
  l1 <- treestats::phylo_to_l(focal_tree)
  a3 <- treestats::branching_times(l1)
  testthat::expect_true(all.equal(a2, a3))
})
