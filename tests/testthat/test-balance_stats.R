context("balance_stats")

test_that("usage", {
  testthat::skip_on_cran() # these tests take very long
  set.seed(42)
  focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0)

  testthat::expect_invisible(
    all_stats <- treestats::calc_balance_stats(focal_tree)
  )
  testthat::expect_equal(length(all_stats), 26)

  testthat::expect_invisible(
    all_stats <- treestats::calc_balance_stats(focal_tree,
                                               normalize = TRUE)
  )
  testthat::expect_equal(length(all_stats), 24)

  focal_tree <- ape::rphylo(n = 23171, birth = 1, death = 0)
  testthat::expect_warning(
    all_stats <- treestats::calc_balance_stats(focal_tree)
  )

  testthat::expect_true(is.na(all_stats$area_per_pair))
  testthat::expect_true(is.na(all_stats$tot_coph))
})
