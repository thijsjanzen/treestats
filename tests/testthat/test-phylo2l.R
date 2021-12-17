context("phylo2L")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 500,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  L1 <- DDD::phylo2L(focal_tree)
  L2 <- treestats::phylo_to_l(focal_tree)

  diff <- L1 - L2

  for (i in 1:4) {
    testthat::expect_lt(mean(diff[, i]), 1e-6)
  }

  # again, but with extinct lineages:
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.3)[[1]]

  L1 <- DDD::phylo2L(focal_tree)
  L2 <- treestats::phylo_to_l(focal_tree)

  diff <- L1 - L2
  for (i in 1:4) {
    testthat::expect_lt(mean(diff[, i]), 1e-6)
  }
})
