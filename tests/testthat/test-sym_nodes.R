context("symNodes")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  a1 <- treestats::sym_nodes(focal_tree)
  a1_check <- treebalance::symNodesI(focal_tree)
  testthat::expect_equal(a1, a1_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::sym_nodes(focal_tree),
                         treestats::sym_nodes(ltab))


  # with extinct species:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.2)[[1]]

  a1 <- treestats::sym_nodes(focal_tree)
  a1_check <- treebalance::symNodesI(focal_tree)
  testthat::expect_equal(a1, a1_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::sym_nodes(focal_tree),
                         treestats::sym_nodes(ltab))
})
