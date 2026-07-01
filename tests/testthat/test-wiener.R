context("wiener")

test_that("usage", {

  calc_using_ape <- function(phy) {
    av2 <- ape::dist.nodes(phy)
    return(sum(av2[lower.tri(av2)]))
  }

  standard_test(treestats::wiener, calc_using_ape)

  test_polytomies(treestats::wiener, calc_using_ape)
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
