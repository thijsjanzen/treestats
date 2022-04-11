context("rogers")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  rog <- treestats::rogers(focal_tree)
  rog_check <- treebalance::rogersI(focal_tree)
  testthat::expect_equal(rog, rog_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::rogers(focal_tree),
                         treestats::rogers(ltab))


  # with extinct species:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.2)[[1]]

  rog <- treestats::rogers(focal_tree)
  rog_check <- treebalance::rogersI(focal_tree)
  testthat::expect_equal(rog, rog_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::rogers(focal_tree),
                         treestats::rogers(ltab))
})
