context("var pair distance")

test_that("usage", {

  test_func <- function(phy) {
    n <- length(phy$tip.label)

    a2 <- stats::cophenetic(phy)
    a2 <- a2[lower.tri(a2)]
    a3 <- var(a2, na.rm = TRUE, use = "everything")
    # var calculates / (n - 1), we use n
    a3 <- a3 * (length(a2) - 1) / length(a2)
    return(a3)
  }

  standard_test(treestats::var_pair_dist, test_func)

  test_polytomies(treestats::var_pair_dist, test_func)
})

test_that("unrooted", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
  a1 <- treestats::var_pair_dist(focal_tree)
  a2 <- treestats::var_pair_dist(ape::unroot(focal_tree))
  testthat::expect_equal(a1, a2)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::var_pair_dist(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::var_pair_dist(list()),
    "input object has to be phylo or ltable"
  )
})
