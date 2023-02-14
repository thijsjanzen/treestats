context("brts_stats")

test_that("usage", {
  testthat::skip_on_cran() # these tests take very long
  set.seed(42)
  focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0)

  testthat::expect_invisible(
    all_stats <- treestats::calc_brts_stats(focal_tree)
  )
  testthat::expect_equal(length(all_stats), 9)
})
