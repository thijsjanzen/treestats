context("beta")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]]

  beta_treestats <- treestats::beta(focal_tree)
  beta_ref       <- apTreeshape::maxlik.betasplit(apTreeshape::as.treeshape(focal_tree))

  testthat::expect_equal(beta_treestats, beta_ref$max_lik, tolerance = 0.05)

  beta_ltab <- treestats::beta(DDD::phylo2L(focal_tree))
  testthat::expect_equal(beta_treestats, beta_ltab)

  focal_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 10, lambda = 1, mu = 0)
  class(focal_tree) <- "multiPhylo"
  testthat::expect_warning (
   beta_treestats <- treestats::beta(focal_tree)
  )


  focal_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0.5)[[1]]
  testthat::expect_warning(
    beta_treestats <- treestats::beta(focal_tree)
  )
})
