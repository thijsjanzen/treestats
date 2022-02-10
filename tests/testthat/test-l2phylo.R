context("L2phylo")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  ltable_1 <- treestats::phylo_to_l(focal_tree)

  tree2 <- treestats::l_to_phylo(ltable_1)
  tree3 <- DDD::L2phylo(ltable_1)

  ax <- geiger::is.extinct(tree2)
  ay <- geiger::is.extinct(tree3)

  testthat::expect_true(all.equal(ax, ay))

  diff_edge_length <- sort(tree2$edge.length) - sort(tree3$edge.length)
  testthat::expect_equal(mean(diff_edge_length), 0.0, tolerance = 0.001)



  # and now with a tree with extinct tips
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.3)[[1]]

  ltable_1 <- treestats::phylo_to_l(focal_tree)

  tree2 <- treestats::l_to_phylo(ltable_1)
  tree3 <- DDD::L2phylo(ltable_1)

  ax <- geiger::is.extinct(tree2)
  ay <- geiger::is.extinct(tree3)

  testthat::expect_true(all.equal(ax, ay))
  diff_edge_length <- sort(tree2$edge.length) - sort(tree3$edge.length)
  testthat::expect_equal(mean(diff_edge_length), 0.0, tolerance = 0.001)

  # and again, but retain extinct lineages
  tree2 <- treestats::l_to_phylo(ltable_1, drop_extinct = FALSE)
  tree3 <- DDD::L2phylo(ltable_1, dropextinct = FALSE)

  ax <- geiger::is.extinct(tree2)
  ay <- geiger::is.extinct(tree3)

  testthat::expect_true(all.equal(ax, ay))
  diff_edge_length <- sort(tree2$edge.length) - sort(tree3$edge.length)
  testthat::expect_equal(mean(diff_edge_length), 0.0, tolerance = 0.001)
})