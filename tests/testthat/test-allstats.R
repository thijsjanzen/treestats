context("all_Stats")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  testthat::expect_invisible(
   all_stats <- treestats::calc_all_stats(focal_tree)
  )
  testthat::expect_equal(length(all_stats), 41)


  focal_tree <- ape::rphylo(n = 23171, birth = 1, death = 0)
  testthat::expect_invisible(
    all_stats <- treestats::calc_all_stats(focal_tree)
  )
  testthat::expect_true(is.na(all_stats$laplac_spectrum_a))
  testthat::expect_true(is.na(all_stats$laplac_spectrum_p))
  testthat::expect_true(is.na(all_stats$laplac_spectrum_e))
  testthat::expect_true(is.na(all_stats$laplac_spectrum_g))

  testthat::expect_true(is.na(all_stats$mntd))
  testthat::expect_true(is.na(all_stats$mpd))
  testthat::expect_true(is.na(all_stats$psv))
  testthat::expect_true(is.na(all_stats$vpd))
  testthat::expect_true(is.na(all_stats$j_stat))
})
