context("phylodiv")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 3,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  div1 <- treestats::phylogenetic_diversity(focal_tree)
  div2 <- sum(focal_tree$edge.length)
  testthat::expect_equal(div1, div2, tolerance = 1e-4)

  ca <- max(treestats::branching_times(focal_tree))
  pds <- treestats::phylogenetic_diversity(focal_tree,
                                           t = seq(0, ca, length.out =  100))
  testthat::expect_equal(length(pds), 100)
  for (i in 2:length(pds)) {
    testthat::expect_gt(pds[i], pds[i - 1])
  }

  # now check sub time
  brts <- ape::branching.times(focal_tree)
  tt <- (brts[1] - brts[2]) / 2
  div1 <- treestats::phylogenetic_diversity(focal_tree, tt)
  div2 <- tt * 2
  testthat::expect_equal(div1[[1]], div2[[1]])

  # now with extinct lineages
  focal_tree <- ape::read.tree(text =
                          "((t1:2.0, t2:2.0):1.0, (t3:1.0, t4:2.0):1.0):1.0;")

  div1 <- treestats::phylogenetic_diversity(focal_tree)
  testthat::expect_equal(div1, 8)

  div1 <- treestats::phylogenetic_diversity(focal_tree, 2.0)
  div2 <- treestats::phylogenetic_diversity(focal_tree, 2.02)
  testthat::expect_gt(div1, div2)

  # now with double extinct lineages
  focal_tree <- ape::read.tree(text =
                      "((:2.0, :2.0):1.0, ((:0.5, :1.0):0.5, :2.0):1.0):1.0;")
  div1 <- treestats::phylogenetic_diversity(focal_tree)
  testthat::expect_equal(div1, 8)
  div1 <- treestats::phylogenetic_diversity(focal_tree, 2)
  div2 <- treestats::phylogenetic_diversity(focal_tree, 2.01)
  testthat::expect_gt(div1, div2)




})
