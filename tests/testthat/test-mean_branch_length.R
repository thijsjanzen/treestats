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



  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::mean_branch_length(focal_tree)
  a2 <- calc_mean_br_r(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(a2,
                         treestats::mean_branch_length(ltab))


  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1 <- treestats::mean_branch_length(focal_tree)
  a2 <- calc_mean_br_r(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
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
