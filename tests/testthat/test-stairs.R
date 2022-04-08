context("stairs")

test_that("usage stairs1", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 30,
                                     numbsim = 1, lambda = 1, mu = 0)[[1]]

  c1 <- treestats::stairs(focal_tree)
  c2 <- phyloTop::stairs(focal_tree)[[1]]
  testthat::expect_equal(c1, c2)

  c3 <- treestats::stairs(treestats::phylo_to_l(focal_tree))
  testthat::expect_equal(c1, c3)


  focal_tree <- TreeSim::sim.bd.taxa(n = 30,
                                     numbsim = 1, lambda = 1, mu = 0.5)[[1]]

  c1 <- treestats::stairs(focal_tree)
  c2 <- phyloTop::stairs(focal_tree)[[1]]
  testthat::expect_equal(c1, c2)

  c3 <- treestats::stairs(treestats::phylo_to_l(focal_tree))
  testthat::expect_equal(c1, c3)
})

test_that("usage stairs2", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 30,
                                     numbsim = 1, lambda = 1, mu = 0)[[1]]

  c1 <- treestats::stairs2(focal_tree)
  c2 <- phyloTop::stairs(focal_tree)[[2]]
  testthat::expect_equal(c1, c2)

  c3 <- treestats::stairs2(treestats::phylo_to_l(focal_tree))
  testthat::expect_equal(c1, c3)


  focal_tree <- TreeSim::sim.bd.taxa(n = 30,
                                     numbsim = 1, lambda = 1, mu = 0.5)[[1]]

  c1 <- treestats::stairs2(focal_tree)
  c2 <- phyloTop::stairs(focal_tree)[[2]]
  testthat::expect_equal(c1, c2)

  c3 <- treestats::stairs2(treestats::phylo_to_l(focal_tree))
  testthat::expect_equal(c1, c3)
})
