context("var pair distance")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  a1 <- treestats::var_pair_dist(focal_tree)

  n <- length(focal_tree$tip.label)
  sample_mat <- matrix(data = 1, nrow = n, ncol = n)
  colnames(sample_mat) <- focal_tree$tip.label

  a2 <- cophenetic(focal_tree)
  a2 <- a2[lower.tri(a2)]
  a3 <- var(a2, na.rm = TRUE, use = "everything")
  a3 <- a3 * (length(a2) - 1) / length(a2)  # var calculates / (n - 1), we use n

  testthat::expect_equal(a1, a3)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::var_pair_dist(focal_tree),
                         treestats::var_pair_dist(ltab))
})
