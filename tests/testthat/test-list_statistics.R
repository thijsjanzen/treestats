context("list_statistics")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 5, birth = 1, death = 0)
  all_stats <- treestats::calc_all_stats(focal_tree)
  names_all_stats <- names(all_stats)
  all_names <- treestats::list_statistics(only_balance_stats = FALSE)
  testthat::expect_equal(sum(all_names %in% names_all_stats),
                         length(all_names))

  testthat::expect_equal(length(names_all_stats),
                         length(all_names))

  all_stats <- treestats::calc_topology_stats(focal_tree)
  names_all_stats <- names(all_stats)
  all_names <- treestats::list_statistics(only_balance_stats = TRUE)
  testthat::expect_equal(sum(all_names %in% names_all_stats),
                         length(all_names))
})
