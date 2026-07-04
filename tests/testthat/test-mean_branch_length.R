context("mean_branch_length")

test_that("usage", {
  calc_mean_br_r <- function(focal_tree) {
    return(mean(focal_tree$edge.length, na.rm = TRUE))
  }

  focal_tree <- ape::read.tree(text = "(1:4,2:4):0;")

  a1 <- treestats::mean_branch_length(focal_tree)
  a2 <- calc_mean_br_r(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::mean_branch_length(focal_tree),
                         treestats::mean_branch_length(ltab))

  a1 <- treestats::mean_branch_length(extant_tree)
  a2 <- calc_mean_br_r(extant_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(extant_tree)
  testthat::expect_equal(a2,
                         treestats::mean_branch_length(ltab))


  # with extinct species:
  a1 <- treestats::mean_branch_length(complete_tree)
  a2 <- calc_mean_br_r(complete_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(complete_tree)
  testthat::expect_equal(a2,
                         treestats::mean_branch_length(ltab))
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::mean_branch_length(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mean_branch_length(list()),
    "input object has to be phylo or ltable"
  )
})
