context("mean pair distance")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::mean_pair_dist(focal_tree)

  n <- length(focal_tree$tip.label)
  sample_mat <- matrix(data = 1, nrow = n, ncol = n)
  colnames(sample_mat) <- focal_tree$tip.label

  a2 <- picante::mpd(sample_mat, cophenetic(focal_tree),
                     abundance.weighted = FALSE)[[1]]
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::mean_pair_dist(focal_tree),
                         treestats::mean_pair_dist(ltab))
})
