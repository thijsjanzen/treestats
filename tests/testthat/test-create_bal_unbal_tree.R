context("create_bal_unbal_trees")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 16, birth = 1, death = 0)
  bal_tree <- treestats::create_fully_balanced_tree(focal_tree)
  unbal_tree <- treestats::create_fully_unbalanced_tree(focal_tree)

  a1 <- treestats::cherries(bal_tree)
  testthat::expect_equal(a1, 8)
  a2 <- treestats::pitchforks(unbal_tree)
  testthat::expect_equal(a2, 1)

  a1 <- treestats::colless(bal_tree)
  a2 <- treestats::colless(unbal_tree)
  testthat::expect_gt(a2, a1)

  ax <- treestats::imbalance_steps(bal_tree)
  ay <- treestats::imbalance_steps(unbal_tree)
  testthat::expect_gt(ax, ay)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::create_fully_balanced_tree(10),
    "This function requires a phylogeny as input"
  )

  testthat::expect_error(
    treestats::create_fully_balanced_tree(list()),
    "This function requires a phylogeny as input"
  )

  testthat::expect_error(
    treestats::create_fully_unbalanced_tree(10),
    "This function requires a phylogeny as input"
  )

  testthat::expect_error(
    treestats::create_fully_unbalanced_tree(list()),
    "This function requires a phylogeny as input"
  )
})
