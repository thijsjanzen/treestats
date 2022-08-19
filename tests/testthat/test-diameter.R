context("diameter")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1_1 <- treestats::diameter(focal_tree, weight = TRUE)
  a1_2 <- treestats::diameter(focal_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 8.750583  # treeCentrality::computeDiameter(focal_tree,
  #                                                   weight = TRUE))
  a2_2 <- 21        # treeCentrality::computeDiameter(focal_tree,
  #                                                      weight = FALSE))

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::diameter(focal_tree, weight = TRUE),
                         treestats::diameter(ltab, weight = TRUE))

  testthat::expect_equal(treestats::diameter(focal_tree, weight = FALSE),
                         treestats::diameter(ltab, weight = FALSE))

  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1_1 <- treestats::diameter(focal_tree, weight = TRUE)
  a1_2 <- treestats::diameter(focal_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 11.21469 # treeCentrality::computeDiameter(focal_tree,
  #                                                  weight = TRUE))
  a2_2 <- 21       # treeCentrality::computeDiameter(focal_tree,
  #                                                  weight = FALSE))


  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::diameter(focal_tree, weight = TRUE),
                         treestats::diameter(ltab,       weight = TRUE))

  testthat::expect_equal(treestats::diameter(focal_tree, weight = FALSE),
                         treestats::diameter(ltab,       weight = FALSE))
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::diameter(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::diameter(list()),
    "input object has to be phylo or ltable"
  )
})
