context("all_statistics")

test_that("usage", {
  testthat::skip_on_cran() # these tests take very long
  testthat::skip_on_ci()
  testthat::skip_on_covr()

  set.seed(42)
  focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0)

  testthat::expect_invisible(
   all_stats <- treestats::calc_all_stats(focal_tree)
  )
  testthat::expect_equal(length(all_stats), 70)

  focal_tree <- ape::rphylo(n = 23175, birth = 1, death = 0)
  all_stats <- treestats::calc_all_stats(focal_tree)

  testthat::expect_true(
    is.na(all_stats[names(all_stats) == "laplace_spectrum_a"]))
  testthat::expect_true(
    is.na(all_stats[names(all_stats) == "laplace_spectrum_p"]))
  testthat::expect_true(
    is.na(all_stats[names(all_stats) == "laplace_spectrum_e"]))
  testthat::expect_true(
    is.na(all_stats[names(all_stats) == "laplace_spectrum_g"]))

  testthat::expect_true(is.na(all_stats[names(all_stats) == "max_laplace"]))
  testthat::expect_true(is.na(all_stats[names(all_stats) == "min_laplace"]))
  testthat::expect_true(is.na(all_stats[names(all_stats) == "max_adj"]))
  testthat::expect_true(is.na(all_stats[names(all_stats) == "min_adj"]))

  testthat::expect_true(is.na(all_stats[names(all_stats) == "vpd"]))
})
