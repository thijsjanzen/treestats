context("wiener")

test_that("usage", {

  calc_using_ape <- function(phy) {
    av2 <- ape::dist.nodes(phy)
    return(sum(av2[lower.tri(av2)]))
  }

  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1_1 <- treestats::wiener(focal_tree, weight = TRUE)
  a1_2 <- treestats::wiener(focal_tree, weight = FALSE)
  a1_3 <- treestats::wiener(focal_tree, normalization = TRUE, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 120552.3 # treeCentrality::computeWienerIndex(focal_tree, #noLint
                   #                                    weight = TRUE)  #noLint
  a2_2 <- 217028   # treeCentrality::computeWienerIndex(focal_tree,  #noLint
                   #                                    weight = FALSE) #noLint
  a2_3 <- 11.01609 # treeCentrality::computeWienerIndex(focal_tree, #noLint
                   #                                    norm = TRUE, #noLint
                   #                                    weight = FALSE) #noLint

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)
  testthat::expect_equal(a1_3, a2_3, tolerance = 0.01)

  testthat::expect_equal(a1_1, calc_using_ape(focal_tree))

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::wiener(focal_tree, weight = TRUE),
                         treestats::wiener(ltab, weight = TRUE))

  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1_1 <- treestats::wiener(focal_tree, weight = TRUE)
  a1_2 <- treestats::wiener(focal_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 216717.3 # treeCentrality::computeWienerIndex(focal_tree,
                   #                                    weight = TRUE)
  a2_2 <- 362732   # treeCentrality::computeWienerIndex(focal_tree,
                   #                                    weight = FALSE)

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::wiener(focal_tree, weight = TRUE),
                         treestats::wiener(ltab,       weight = TRUE))
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::wiener(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::wiener(list()),
    "input object has to be phylo or ltable"
  )
})
