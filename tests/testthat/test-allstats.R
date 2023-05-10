context("all_statistics")

test_that("usage", {
  testthat::skip_on_cran() # these tests take very long
   testthat::skip("fast checking")
   set.seed(42)
  focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0)

  testthat::expect_invisible(
   all_stats <- treestats::calc_all_stats(focal_tree)
  )
  testthat::expect_equal(length(all_stats), 54)

  focal_tree <- ape::rphylo(n = 23175, birth = 1, death = 0)
  testthat::expect_warning(
    all_stats <- treestats::calc_all_stats(focal_tree)
  )
  testthat::expect_true(is.na(all_stats$laplac_spectrum_a))
  testthat::expect_true(is.na(all_stats$laplac_spectrum_p))
  testthat::expect_true(is.na(all_stats$laplac_spectrum_e))
  testthat::expect_true(is.na(all_stats$laplac_spectrum_g))

  testthat::expect_true(is.na(all_stats$mntd))
  testthat::expect_true(is.na(all_stats$psv))
  testthat::expect_true(is.na(all_stats$vpd))
  testthat::expect_true(is.na(all_stats$tot_coph))
  testthat::expect_true(is.na(all_stats$area_per_pair))
})
