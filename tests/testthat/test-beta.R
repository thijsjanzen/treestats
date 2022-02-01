context("beta")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  beta_treestats <- treestats::beta_statistic(focal_tree)
  beta_ref       <- apTreeshape::maxlik.betasplit(
    apTreeshape::as.treeshape(focal_tree))

  testthat::expect_equal(beta_treestats, beta_ref$max_lik, tolerance = 0.05)

  brts <- as.vector(ape::branching.times(focal_tree))

  bal_tree <- nodeSub::create_balanced_tree(brts)

  beta_treestats <- treestats::beta_statistic(bal_tree)
  beta_ref       <- apTreeshape::maxlik.betasplit(bal_tree)

  testthat::expect_equal(beta_treestats, beta_ref$max_lik, tolerance = 0.05)


  unbal_tree <- nodeSub::create_unbalanced_tree(brts)

  beta_treestats <- treestats::beta_statistic(unbal_tree)
  beta_ref       <- apTreeshape::maxlik.betasplit(unbal_tree)

  testthat::expect_equal(beta_ref$max_lik, -2, tolerance = 0.01)
  testthat::expect_equal(beta_treestats, beta_ref$max_lik, tolerance = 0.05)
  testthat::expect_equal(beta_treestats, -2, tolerance = 0.01)
})


# abuse
test_that("abuse", {
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
  numbsim = 1,
  lambda = 1, mu = 0.5)[[1]]

  testthat::expect_error(
    # throws error, because beta can only be calculated on extant tree
    treestats::beta_statistic(focal_tree)
  )

  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  testthat::expect_output(
    treestats::beta_statistic(focal_tree, algorithm = "none"),
    "no algorithm chosen")
})

