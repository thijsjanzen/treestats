context("inv_branch_dist")

test_that("usage", {
  # let's create an artifical tree for which we know the values
  set.seed(42)
  focal_tree <- ape::rphylo(n = 3, 1, 0)
  focal_tree <- stats::reorder(focal_tree)
  focal_tree$edge.length <- c(0.5, 1, 1, 1.5)
  expected_vals <- c(1 + 1 / 0.5, 1 + 1 / 0.5, 1 / 1.5)

  inv_dist <- treestats::inv_branch_dist(focal_tree, add_one = FALSE)
  testthat::expect_equal(as.vector(inv_dist), expected_vals)

  ltab <- treestats::phylo_to_l(focal_tree)
  inv_dist_ltab <- treestats::inv_branch_dist(ltab, add_one = FALSE)
  testthat::expect_equal(sum(as.vector(inv_dist_ltab)), sum(expected_vals))

  # now we add one
  expected_vals <- c(1 / (1 + 1) + 1 / 1.5, 1 / (1 + 1) + 1 / 1.5, 1 / 2.5)

  inv_dist <- treestats::inv_branch_dist(focal_tree, add_one = TRUE)
  testthat::expect_equal(as.vector(inv_dist), expected_vals)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::inv_branch_dist(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::inv_branch_dist(list()),
    "input object has to be phylo or ltable"
  )
})
