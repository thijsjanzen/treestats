context("branching_times")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 500, birth = 1, death = 0)

  a1 <- ape::branching.times(focal_tree)
  a2 <- treestats::branching_times(focal_tree)

  testthat::expect_equal(a1, a2)

  # again, but with extinct lineages:
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.3, fossils = TRUE)

  a1 <- ape::branching.times(focal_tree)

  # oddly enough, ape::branching.times does not rescale correctly:
  if (min(a1) < 0) {
    a1 <- a1 - min(a1)
  }

  a2 <- treestats::branching_times(focal_tree)

  diff_brts <- abs(sort(a1) - sort(a2))
  mean_diff <- mean(diff_brts)
  testthat::expect_lte(mean_diff, min(a2))

  # using ltable
  l1 <- treestats::phylo_to_l(focal_tree)
  a3 <- treestats::branching_times(l1)
  testthat::expect_true(sum(a2 - a3) == 0)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::branching_times(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::branching_times(list()),
    "input object has to be phylo or ltable"
  )
})
