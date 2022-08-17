context("betweenness")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1_1 <- treestats::max_betweenness(focal_tree)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 12795  # max(treeCentrality::computeBetweenness(focal_tree)) #nolint

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::diameter(focal_tree, weight = TRUE),
                         treestats::diameter(ltab, weight = TRUE))

  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1_1 <- treestats::max_betweenness(focal_tree)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 20315 # max(treeCentrality::computeBetweenness(focal_tree)) #nolint


  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::diameter(focal_tree, weight = TRUE),
                         treestats::diameter(ltab,       weight = TRUE))
})
