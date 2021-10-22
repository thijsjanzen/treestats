context("nltt")

test_that("usage", {
  set.seed(42)
  tree1 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]]
  tree2 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 0.5, mu = 0)[[1]]

  nltt <- treestats::nLTT(tree1, tree2)
  nltt_check <- nLTT::nLTTstat(tree1, tree2)
  testthat::expect_equal(nltt, nltt_check)

  empty_tree <- TreeSim::sim.bd.taxa(n = 2, numbsim = 1, lambda = 1, mu = 0)[[1]]
  nltt_base <- treestats::nLTT_base(tree1)
  nltt_base2 <- treestats::nLTT(tree1, empty_tree)
  nltt_base3 <-  nLTT::nLTTstat(tree1, empty_tree)

  testthat::expect_equal(nltt_base, nltt_base2, tolerance = 1e-3)
  testthat::expect_equal(nltt_base, nltt_base3, tolerance = 1e-3)

  empty_tree <- ape::read.tree(text = "(1:4,2:4):0;")
  nltt_stat <- nLTT::nltt_diff(tree1, empty_tree)
  testthat::expect_equal(nltt_stat, nltt_base, tolerance = 1e-3)
  testthat::expect_equal(nltt_stat, nltt_base2, tolerance = 1e-3)
})
