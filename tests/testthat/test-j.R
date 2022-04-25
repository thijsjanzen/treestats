context("mean pair distance")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  a1 <- treestats::entropy_j(focal_tree)

  n <- length(focal_tree$tip.label)
  sample_mat <- matrix(data = 1, nrow = n, ncol = n)
  colnames(sample_mat) <- focal_tree$tip.label

  a2 <- picante::mpd(sample_mat, cophenetic(focal_tree),
                     abundance.weighted = FALSE)[[1]]
  testthat::expect_equal(a1, a2 / n)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::mean_pair_dist(focal_tree),
                         treestats::mean_pair_dist(ltab))
})