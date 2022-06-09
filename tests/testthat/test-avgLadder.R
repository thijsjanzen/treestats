context("avgLadder")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1, lambda = 1, mu = 0)[[1]]

  c1 <- treestats::avgLadder(focal_tree)
  c2 <- phyloTop::avgLadder(focal_tree)
  testthat::expect_equal(c1, c2)

  c3 <- treestats::avgLadder(treestats::phylo_to_l(focal_tree))
  testthat::expect_equal(c1, c3)


  focal_tree <- TreeSim::sim.bd.taxa(n = 30,
                                     numbsim = 1, lambda = 1, mu = 0.5)[[1]]

  c1 <- treestats::avgLadder(focal_tree)
  c2 <- phyloTop::avgLadder(focal_tree)
  testthat::expect_equal(c1, c2)

  c3 <- treestats::avgLadder(treestats::phylo_to_l(focal_tree))
  testthat::expect_equal(c1, c3)
})
