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


test_that("normalization", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::max_closeness(focal_tree, weight = TRUE)
  c2 <- treestats::max_closeness(focal_tree,  weight = TRUE,
                                 normalization = "tips")
  testthat::expect_lt(c1, c2)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::max_closeness(focal_tree))
    stats2 <- c(stats2, treestats::max_closeness(focal_tree,
                                                 normalization = "tips"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(abs(a2), abs(a1))
  testthat::expect_lt(abs(a2), 0.5)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::max_closeness(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::max_closeness(list()),
    "input object has to be phylo or ltable"
  )
})
