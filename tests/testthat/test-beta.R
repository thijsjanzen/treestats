context("beta")

test_that("usage", {
  if (requireNamespace("nodeSub") &&
      requireNamespace("apTreeshape")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

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


    # test methods
    comp_methods <- function(sim_tree) {
      beta1 <- treestats::beta_statistic(sim_tree, algorithm = "COBYLA")
      beta2 <- treestats::beta_statistic(sim_tree, algorithm = "subplex")
      beta3 <- treestats::beta_statistic(sim_tree, algorithm = "simplex")
      beta4 <- apTreeshape::maxlik.betasplit(sim_tree)$max_lik
      testthat::expect_equal(beta1, beta4, tolerance = 0.1)
      testthat::expect_equal(beta1, beta2, tolerance = 0.1)
      testthat::expect_equal(beta1, beta3, tolerance = 0.1)

      ltab <- treestats::phylo_to_l(sim_tree)
      beta5 <- treestats::beta_statistic(ltab)
      testthat::expect_equal(beta1, beta5, tolerance = 0.1)
    }

    comp_methods(bal_tree)
    comp_methods(focal_tree)
  }
})

# abuse
test_that("abuse", {
  focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0.5, fossils = TRUE)

  testthat::expect_error(
    # throws error, because beta can only be calculated on extant tree
    treestats::beta_statistic(focal_tree)
  )

  focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0)

  testthat::expect_output(
    treestats::beta_statistic(phy = focal_tree, algorithm = "none"),
    "no algorithm chosen")
})
