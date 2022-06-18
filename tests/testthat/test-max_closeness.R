context("max_closeness")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1_1 <- treestats::max_closeness(focal_tree, weight = TRUE)
  a1_2 <- treestats::max_closeness(focal_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 0.001400582  # max(treeCentrality::computeCloseness(focal_tree,
                       #                                      weight = TRUE))
  a2_2 <- 0.0007843137 # max(treeCentrality::computeCloseness(focal_tree,
                       #                                      weight = FALSE))

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::max_closeness(focal_tree, weight = TRUE),
                         treestats::max_closeness(ltab, weight = TRUE))

  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1_1 <- treestats::max_closeness(focal_tree, weight = TRUE)
  a1_2 <- treestats::max_closeness(focal_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 0.001018909  # max(treeCentrality::computeCloseness(focal_tree,
  #                                      weight = TRUE))
  a2_2 <- 0.0006035003 # max(treeCentrality::computeCloseness(focal_tree,
  #                                      weight = FALSE))


  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::max_closeness(focal_tree, weight = TRUE),
                         treestats::max_closeness(ltab,       weight = TRUE))
})
